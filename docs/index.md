# PhytClust

**Tree-aware clustering for rooted phylogenies.**

PhytClust takes a rooted phylogenetic tree and splits its leaves into clusters — with the guarantee that every cluster is a true clade. No post-hoc filtering, no heuristic merging: monophyly is baked into the optimization from the start.

The core algorithm is dynamic programming directly on the tree topology. It finds the partition of leaves into *k* groups that minimizes within-cluster branch length, while never breaking monophyly. Because the DP explores the full solution space for a given *k*, the result is exact — not an approximation.

---

## Who is this for?

PhytClust was built for researchers who work with phylogenies and need to group taxa into meaningful clusters — whether for downstream analysis, visualization, or summarizing large trees. Typical use cases:

- **Cancer genomics**: grouping clones from tumor phylogenies into subclonal populations
- **Comparative genomics**: partitioning gene or species trees into coherent groups
- **Ecology / metagenomics**: summarizing diversity at multiple resolutions
- **Any tree**: if you have a rooted Newick tree, PhytClust can cluster it

## What can it do?

### Three clustering modes

You can tell PhytClust exactly how many clusters you want, or let it figure that out for you:

- **Exact *k***: "Give me exactly 5 clusters." Useful when biology or a downstream pipeline dictates the number.
- **Global peak search**: "Find the best *k* values." PhytClust computes a score curve over all possible *k* and picks the most prominent peaks — the *k* values where the tree has natural breakpoints.
- **Multi-resolution**: "Show me the tree at different granularities." Returns one representative *k* per logarithmic bin, so you get a coarse-to-fine panel in a single run.

### Practical controls

Real trees are messy. PhytClust handles the common headaches:

- **Minimum cluster size** — hard constraint so you don't get singleton noise
- **Outlier handling** — small clusters can be flagged as outliers, with optional preference for fewer outlier groups
- **Polytomies** — native support for multifurcations (hard or soft mode), no need to resolve them beforehand
- **Zero-length edges** — optional guard against splitting clades joined by zero-length branches
- **Outgroups and rooting** — exclude an outgroup taxon or re-root on the fly

### Three interfaces

- **CLI** — one command, immediate results. Good for quick runs and shell pipelines.
- **Python API** — full programmatic control for notebooks, batch pipelines, and custom analysis.
- **Web GUI** (experimental) — paste a Newick string and explore clusters interactively in the browser.

### Publication-ready output

Every run can produce:

| File | What it is |
|------|-----------|
| `phytclust_results.tsv` | Leaf-to-cluster assignments (one or multiple *k* values) |
| `scores.png` | Score curve with annotated peaks |
| `tree_k{K}.png` | Colored tree for each selected *k* |
| `peaks_by_rank.txt` | Selected *k* values in rank order |

Plots use a colorblind-safe palette and are configurable down to font sizes, marker sizes, and axis scales.

---

## Where to start

If you're new here, go through these pages in order:

1. **[Getting Started](getting-started.md)** — install, verify, run your first clustering
2. **[How It Works](concepts.md)** — what the algorithm actually does (no math prerequisites)
3. **[CLI Tutorial](tutorials/cli-end-to-end.md)** — practical walkthrough from first run to tuned output
4. **[Python Tutorial](tutorials/python-end-to-end.md)** — same journey, but in scripts and notebooks

When you need to look something up:

- [CLI Flags Reference](reference/cli.md)
- [Configuration Reference](reference/configuration.md)
- [Python API Reference](reference/api.md)
