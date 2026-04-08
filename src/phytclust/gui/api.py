# src/phytclust/gui/api.py

import logging
import tempfile
from pathlib import Path
from typing import Any, Optional

from fastapi import FastAPI, HTTPException
from fastapi.responses import HTMLResponse, PlainTextResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from phytclust.algo.core import PhytClust
from phytclust.config import OutlierConfig, PeakConfig

logger = logging.getLogger("phytclust.gui")

# ------------------------------------------------------
# Create the FastAPI app
# ------------------------------------------------------
app = FastAPI(title="PhytClust Web API")

STATIC_DIR = Path(__file__).parent / "static"
app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")

# Keep the last PhytClust instance and result in memory (simple session)
LAST_PC: Optional[PhytClust] = None
LAST_RESULT: Optional[dict[str, Any]] = None
LAST_NEWICK: Optional[str] = None


# ------------------------------------------------------
# Serve index.html at root URL
# ------------------------------------------------------
@app.get("/", response_class=HTMLResponse)
def serve_frontend():
    index_path = Path(__file__).parent / "templates" / "index.html"
    return index_path.read_text()


# ------------------------------------------------------
# Request model
# ------------------------------------------------------
class PhytclustRequest(BaseModel):
    newick: str
    mode: str = "global"
    k: Optional[int] = None
    top_n: int = 1
    by_resolution: bool = False
    num_bins: Optional[int] = None
    max_k: Optional[int] = None
    max_k_limit: Optional[float] = None
    outgroup: Optional[str] = None
    root_taxon: Optional[str] = None
    min_cluster_size: Optional[int] = None
    compute_all_clusters: bool = False

    # Branch-support options
    use_branch_support: bool = False
    min_support: Optional[float] = None
    support_weight: Optional[float] = None

    # Outlier config
    outlier_size_threshold: Optional[int] = None
    outlier_prefer_fewer: bool = False
    outlier_ratio_weight: Optional[float] = None
    outlier_ratio_mode: Optional[str] = None

    # Zero-length / polytomy handling
    no_split_zero_length: bool = False
    optimize_polytomies: bool = True

    # Peak config
    lambda_weight: float = 0.7
    ranking_mode: str = "adjusted"
    min_prominence: Optional[float] = None
    use_relative_prominence: bool = False
    boundary_window_size: Optional[int] = None
    boundary_ratio_threshold: Optional[float] = None


# ------------------------------------------------------
# Helpers
# ------------------------------------------------------
def _leaf_key_to_name(x: Any) -> str:
    name = getattr(x, "name", None)
    if name is not None:
        return str(name)
    return str(x)


def _serialize_clusters(raw_clusters: Any) -> list[dict[str, int]]:
    out: list[dict[str, int]] = []

    if raw_clusters is None:
        return out

    if isinstance(raw_clusters, dict):
        out.append({_leaf_key_to_name(k): int(v) for k, v in raw_clusters.items()})
        return out

    for cmap in raw_clusters:
        out.append({_leaf_key_to_name(k): int(v) for k, v in cmap.items()})

    return out


# ------------------------------------------------------
# API endpoint
# ------------------------------------------------------
@app.post("/api/run")
def run_phytclust(req: PhytclustRequest):
    global LAST_PC, LAST_RESULT, LAST_NEWICK
    notes: list[str] = []

    if not req.newick.strip():
        raise HTTPException(status_code=400, detail="Empty Newick string.")

    mode = req.mode.lower()
    if mode not in {"k", "global", "resolution"}:
        raise HTTPException(
            status_code=400, detail="mode must be 'k', 'global', or 'resolution'."
        )

    if (
        LAST_NEWICK is not None
        and LAST_PC is not None
        and str(LAST_NEWICK) == str(req.newick)
    ):
        pc = LAST_PC
        pc.compute_all_clusters = req.compute_all_clusters
        notes.append("Reusing previous PhytClust instance.")
    else:
        try:
            # Build outlier config
            outlier_kwargs: dict[str, Any] = {}
            if req.outlier_size_threshold is not None:
                outlier_kwargs["size_threshold"] = req.outlier_size_threshold
            if req.outlier_prefer_fewer:
                outlier_kwargs["prefer_fewer"] = True
            if req.outlier_ratio_weight is not None:
                outlier_kwargs["ratio_weight"] = req.outlier_ratio_weight
            if req.outlier_ratio_mode is not None:
                outlier_kwargs["ratio_mode"] = req.outlier_ratio_mode

            kwargs: dict[str, Any] = dict(
                tree=req.newick,
                outgroup=req.outgroup,
                compute_all_clusters=req.compute_all_clusters,
                use_branch_support=req.use_branch_support,
                no_split_zero_length=req.no_split_zero_length,
                optimize_polytomies=req.optimize_polytomies,
            )
            if req.root_taxon:
                kwargs["root_taxon"] = req.root_taxon
            if req.min_cluster_size is not None:
                kwargs["min_cluster_size"] = req.min_cluster_size
            if req.min_support is not None:
                kwargs["min_support"] = req.min_support
            if req.support_weight is not None:
                kwargs["support_weight"] = req.support_weight
            if outlier_kwargs:
                kwargs["outlier"] = OutlierConfig(**outlier_kwargs)

            pc = PhytClust(**kwargs)
        except Exception as e:
            raise HTTPException(status_code=400, detail=f"Could not parse Newick: {e}")

    try:
        # Build peak config
        peak_kwargs: dict[str, Any] = {"lambda_weight": req.lambda_weight}
        if req.ranking_mode:
            peak_kwargs["ranking_mode"] = req.ranking_mode
        if req.min_prominence is not None:
            peak_kwargs["min_prominence"] = req.min_prominence
        if req.use_relative_prominence:
            peak_kwargs["use_relative_prominence"] = True
        if req.boundary_window_size is not None:
            peak_kwargs["boundary_window_size"] = req.boundary_window_size
        if req.boundary_ratio_threshold is not None:
            peak_kwargs["boundary_ratio_threshold"] = req.boundary_ratio_threshold

        cfg = PeakConfig(**peak_kwargs)

        if mode == "k":
            if req.k is None:
                raise HTTPException(
                    status_code=400, detail="k must be provided in k-mode."
                )
            result = pc.run(
                k=req.k,
                top_n=1,
                by_resolution=False,
                plot_scores=False,
                peak_config=cfg,
            )

        elif mode == "resolution":
            result = pc.run(
                k=None,
                top_n=req.top_n,
                by_resolution=True,
                num_bins=req.num_bins,
                max_k=req.max_k,
                max_k_limit=req.max_k_limit,
                plot_scores=False,
                peak_config=cfg,
            )

        elif mode == "global":
            result = pc.run(
                k=None,
                top_n=req.top_n,
                by_resolution=False,
                max_k=req.max_k,
                max_k_limit=req.max_k_limit,
                plot_scores=False,
                peak_config=cfg,
            )
        else:
            raise HTTPException(status_code=400, detail="Invalid mode.")

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"PhytClust error: {e}")

    clusters_json = _serialize_clusters(result.get("clusters"))
    scores = result.get("scores")
    if scores is not None:
        try:
            scores = list(scores)
        except Exception:
            pass

    # When compute_all_clusters is on, include every cached k->cluster mapping
    all_clusters_json = None
    all_ks = None
    if req.compute_all_clusters and pc.clusters:
        sorted_ks = sorted(pc.clusters.keys())
        all_ks = sorted_ks
        all_clusters_json = _serialize_clusters([pc.clusters[kv] for kv in sorted_ks])

    payload = {
        "mode": result.get("mode"),
        "k": result.get("k"),
        "ks": result.get("ks"),
        "peaks": result.get("peaks"),
        "outgroup": result.get("outgroup"),
        "notes": [str(n) for n in notes if n is not None],
        "newick": req.newick,
        "scores": scores,
        "clusters": clusters_json,
        "all_clusters": all_clusters_json,
        "all_ks": all_ks,
    }

    LAST_PC = pc
    LAST_NEWICK = req.newick
    LAST_RESULT = payload
    return payload


# ------------------------------------------------------
# Save endpoint
# ------------------------------------------------------
class SaveRequest(BaseModel):
    results_dir: str
    filename: str = "phytclust_results.tsv"
    top_n: int = 1
    outlier: bool = True
    output_all: bool = False


@app.post("/api/save")
def save_results(req: SaveRequest):
    global LAST_PC, LAST_RESULT
    if LAST_PC is None or LAST_RESULT is None:
        raise HTTPException(
            status_code=400,
            detail="No PhytClust result available. Please run PhytClust before saving.",
        )

    try:
        out_dir = Path(req.results_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        raise HTTPException(
            status_code=400,
            detail=f"Cannot create output directory '{req.results_dir}': {e}",
        )

    try:
        LAST_PC.save(
            results_dir=req.results_dir,
            top_n=req.top_n,
            filename=req.filename,
            outlier=req.outlier,
            output_all=req.output_all,
        )
        return {
            "status": "ok",
            "results_dir": req.results_dir,
            "filename": req.filename,
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to save results: {e}")


class ExportTSVRequest(BaseModel):
    filename: str = "phytclust_results.tsv"
    top_n: int = 1
    outlier: bool = True
    output_all: bool = False


@app.post("/api/export_tsv")
def export_tsv(req: ExportTSVRequest):
    global LAST_PC, LAST_RESULT
    if LAST_PC is None or LAST_RESULT is None:
        raise HTTPException(
            status_code=400,
            detail="No PhytClust result available. Please run PhytClust before saving.",
        )

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            LAST_PC.save(
                results_dir=tmpdir,
                top_n=req.top_n,
                filename=req.filename,
                outlier=req.outlier,
                output_all=req.output_all,
            )
            tsv_path = Path(tmpdir) / req.filename
            if not tsv_path.exists():
                raise HTTPException(status_code=500, detail="TSV file not generated.")
            text = tsv_path.read_text(encoding="utf-8")
            return PlainTextResponse(
                content=text, media_type="text/tab-separated-values"
            )
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to export TSV: {e}")
