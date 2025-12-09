# src/phytclust/gui/api.py

from pathlib import Path
from typing import Any, Dict, List, Optional

from fastapi import FastAPI, HTTPException
from fastapi.responses import HTMLResponse, PlainTextResponse
from pydantic import BaseModel

from phytclust.algo.core import PhytClust
import tempfile
from pathlib import Path

# ------------------------------------------------------
# Create the FastAPI app
# ------------------------------------------------------
app = FastAPI(title="PhytClust Web API")

from fastapi.staticfiles import StaticFiles

STATIC_DIR = Path(__file__).parent / "static"
app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")

# Keep the last PhytClust instance and result in memory (simple session)
LAST_PC: Optional[PhytClust] = None
LAST_RESULT: Optional[Dict[str, Any]] = None


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
    lambda_weight: float = 0.7
    outgroup: Optional[str] = None


# ------------------------------------------------------
# Helpers
# ------------------------------------------------------
def _leaf_key_to_name(x: Any) -> str:
    name = getattr(x, "name", None)
    if name is not None:
        return str(name)
    return str(x)


def _serialize_clusters(raw_clusters: Any) -> List[Dict[str, int]]:
    out: List[Dict[str, int]] = []

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
    if not req.newick.strip():
        raise HTTPException(status_code=400, detail="Empty Newick string.")

    mode = req.mode.lower()
    if mode not in {"k", "global", "resolution"}:
        raise HTTPException(
            status_code=400, detail="mode must be 'k', 'global', or 'resolution'."
        )

    try:
        pc = PhytClust(tree=req.newick, outgroup=req.outgroup)
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Could not parse Newick: {e}")

    try:
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
                lambda_weight=req.lambda_weight,
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
                lambda_weight=req.lambda_weight,
            )

        elif mode == "global":
            result = pc.run(
                k=None,
                top_n=req.top_n,
                by_resolution=False,
                max_k=req.max_k,
                max_k_limit=req.max_k_limit,
                plot_scores=False,
                lambda_weight=req.lambda_weight,
            )
        else:
            raise HTTPException(status_code=400, detail="Invalid mode.")

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"PhytClust error: {e}")

    clusters_json = _serialize_clusters(result.get("clusters"))
    scores = result.get("scores")
    if scores is not None:
        try:
            scores = list(scores)
        except Exception:
            pass

    payload = {
        "mode": result.get("mode"),
        "k": result.get("k"),
        "ks": result.get("ks"),
        "peaks": result.get("peaks"),
        "outgroup": result.get("outgroup"),
        "scores": scores,
        "clusters": clusters_json,
        "newick": req.newick,
    }

    # store last instance and result for use by /api/save
    global LAST_PC, LAST_RESULT
    LAST_PC = pc
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
        raise HTTPException(status_code=400, detail="No PhytClust result available. Please run PhytClust before saving.")

    try:
        # Ensure directory exists
        out_dir = Path(req.results_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Cannot create output directory '{req.results_dir}': {e}")

    try:
        LAST_PC.save(
            results_dir=req.results_dir,
            top_n=req.top_n,
            filename=req.filename,
            outlier=req.outlier,
            output_all=req.output_all,
        )
        return {"status": "ok", "results_dir": req.results_dir, "filename": req.filename}
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
        raise HTTPException(status_code=400, detail="No PhytClust result available. Please run PhytClust before saving.")

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
            return PlainTextResponse(content=text, media_type="text/tab-separated-values")
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to export TSV: {e}")
