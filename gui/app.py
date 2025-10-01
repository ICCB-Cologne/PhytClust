import base64
import io
import plot_phylo

import matplotlib.pyplot as plt

from Bio import Phylo
from phytclust import PhytClust
from flask import Flask, render_template, request, jsonify

app = Flask(__name__)


def process_newick_string(newick_string: str) -> dict[str, str]:
    try:
        quoted_newick_string: str = f"\"{newick_string.rstrip(";")}\""
        tree = Phylo.read(
            io.StringIO(quoted_newick_string), "newick"
        )

        phytclust = PhytClust(tree, should_plot_scores=False)
        clusters = phytclust.best_global(plot_scores=False)

        fig, ax = plt.subplots(figsize=(8, 6))
        plot_phylo.plot_phylo(newick_string, ax)
        buf = io.BytesIO()
        plt.savefig(buf, format="png", bbox_inches="tight")
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode("utf-8")
        plt.close(fig)

        return {
            "input": newick_string,
            "clusters": str(clusters[0]),
            "tree_plot": img_base64,
        }
    except Exception as e:
        return {
            "error": str(e),
        }


@app.route("/")
def index():
    return render_template(
        "index.html"
    )


@app.route("/run-phytclust", methods=["POST"])
def run_phytclust():
    data: dict[str, object] = request.get_json()
    newick_string: str | None = data.get("newick_string", None)

    if not newick_string:
        return jsonify({"error": "no input string provided"}), 400

    result: dict[str, str] = process_newick_string(newick_string)
    return jsonify(result)


if __name__ == "__main__":
    app.run(debug=True, host="127.0.0.1", port=8008)
