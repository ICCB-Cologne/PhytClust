# src/phytclust/gui/api.py

from pathlib import Path
from typing import Any, Dict, List, Optional

from fastapi import FastAPI, HTTPException
from fastapi.responses import HTMLResponse
from pydantic import BaseModel

from phytclust.algo.core import PhytClust

# ------------------------------------------------------
# Create the FastAPI app
# ------------------------------------------------------
app = FastAPI(title="PhytClust Web API")

from fastapi.staticfiles import StaticFiles

STATIC_DIR = Path(__file__).parent / "static"
app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")


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
    max_k: Optional[int] = None
    max_k_limit: Optional[float] = None
    lambda_weight: float = 0.7


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
        pc = PhytClust(tree=req.newick)
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
                num_bins=None,
                max_k=req.max_k,
                max_k_limit=req.max_k_limit,
                plot_scores=False,
                lambda_weight=req.lambda_weight,
            )

        else:  # global mode
            result = pc.run(
                k=None,
                top_n=req.top_n,
                by_resolution=False,
                max_k=req.max_k,
                max_k_limit=req.max_k_limit,
                plot_scores=False,
                lambda_weight=req.lambda_weight,
            )

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"PhytClust error: {e}")

    clusters_json = _serialize_clusters(result.get("clusters"))
    scores = result.get("scores")
    if scores is not None:
        try:
            scores = list(scores)
        except Exception:
            pass

    return {
        "mode": result.get("mode"),
        "k": result.get("k"),
        "ks": result.get("ks"),
        "peaks": result.get("peaks"),
        "scores": scores,
        "clusters": clusters_json,
        "newick": req.newick,
    }
