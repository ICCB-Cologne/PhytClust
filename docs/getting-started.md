# Getting Started

This page gets you from zero to your first clustered tree in about five minutes.

## Requirements

- **Python 3.10** or newer
- A rooted tree in **Newick format** (PhytClust ships with example trees if you don't have one handy)

## Install

The simplest route:

```bash
pip install phytclust
```

If you want to work from source (e.g. to contribute or use the latest development version):

```bash
git clone https://github.com/ICCB-Cologne/PhytClust.git
cd PhytClust
pip install -e ".[dev]"
```

Check that it worked:

```bash
phytclust --version
phytclust --help
```

## Your first run

Let's cluster the bundled sample tree into 5 groups:

```bash
phytclust examples/sample_tree.nwk --k 5 --save-fig --out-dir results/quickstart
```

This tells PhytClust: *"Take this tree, split it into exactly 5 monophyletic clusters, save the figures, and put everything in `results/quickstart/`."*

You'll get:

- **`phytclust_results.tsv`** — which leaf belongs to which cluster
- **`tree_k5.png`** — the tree with clusters colored in
- **`scores.png`** — a score curve showing how good each possible *k* is (useful context even in exact-*k* mode)

Open the tree image to see your clusters. Each color is a clade.

## Let PhytClust pick *k* for you

Often you don't know how many clusters to expect. PhytClust can scan across all possible *k* values and find the ones where the tree has natural breakpoints:

```bash
phytclust examples/sample_tree.nwk --top-n 3 --save-fig --out-dir results/global
```

This returns the top 3 *k* values by peak prominence, along with cluster assignments and tree plots for each. Check `peaks_by_rank.txt` to see which *k* values were selected and in what order.

## Multi-resolution: a panel of granularities

If you want a bird's-eye view at multiple scales — coarse through fine — use resolution mode:

```bash
phytclust examples/sample_tree.nwk --resolution --bins 4 --save-fig --out-dir results/resolution
```

This divides the *k* range into 4 logarithmic bins and picks the best *k* in each bin. You get a compact panel instead of dozens of neighboring solutions.

## Rooting

PhytClust expects a rooted tree. If yours is unrooted, you can root it on the fly:

```bash
# root on a named taxon
phytclust tree.nwk --k 5 --root-taxon "species_A"

# midpoint rooting
phytclust tree.nwk --k 5 --root-taxon midpoint
```

## Web GUI (experimental)

If you prefer a visual interface, PhytClust includes an experimental web GUI. Install the extra dependencies and start it:

```bash
pip install "uvicorn[standard]" fastapi
uvicorn phytclust.gui.api:app --reload --host 127.0.0.1 --port 8000
```

Then open [http://127.0.0.1:8000](http://127.0.0.1:8000) in your browser. You can paste a Newick string directly and explore clusters interactively.

!!! note
    The GUI is a prototype — great for quick exploration, but for reproducible analysis use the CLI or Python API.

## What's next?

Now that you have PhytClust running, pick your path:

- **Want to understand what's happening under the hood?** Read [How It Works](concepts.md) — no math prerequisites, just intuition.
- **Prefer the command line?** The [CLI Tutorial](tutorials/cli-end-to-end.md) walks you through progressively more advanced runs.
- **Working in Python / notebooks?** The [Python Tutorial](tutorials/python-end-to-end.md) covers the same ground with full scripting control.
- **Just need to look up a flag?** Jump to the [CLI Reference](reference/cli.md) or [Configuration Reference](reference/configuration.md).
