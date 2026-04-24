# CLI Tutorial: From First Run to Tuned Output

This tutorial walks through a real workflow with the `phytclust` command line tool. Each section builds on the last, so by the end you'll know how to handle most situations you'll encounter in practice.

We'll use the bundled example trees throughout. If you want to follow along with your own data, just swap the file path.

---

## 1. Start simple: exact-*k* clustering

When you already know how many clusters you need — maybe from prior biology, or because a downstream tool expects a fixed number — use `--k`:

```bash
phytclust examples/sample_tree.nwk --k 5 --save-fig --out-dir results/cli_baseline
```

Open the output folder. You'll find:

- **`phytclust_results.tsv`** — one row per leaf, showing which cluster it belongs to
- **`tree_k5.png`** — the tree with each cluster colored differently
- **`scores.png`** — the score curve across all *k* values

Even in exact-*k* mode, the score plot is worth a glance. If *k*=5 sits at a peak, that's reassuring — it means the tree has a natural breakpoint there. If it's in a flat region, you might want to explore nearby values.

## 2. Let PhytClust find good *k* values

Most of the time you don't know the right *k* in advance. The global peak search scans the full score curve and picks the *k* values where clustering quality jumps:

```bash
phytclust examples/sample_tree.nwk \
  --top-n 3 \
  --max-k 120 \
  --save-fig \
  --out-dir results/cli_global
```

What to look at:

- **`peaks_by_rank.txt`** — the selected *k* values, ranked by prominence. The top one is usually the strongest signal.
- **`scores.png`** — are the peaks well-separated, or are they clustered together? Well-separated peaks suggest distinct structural scales in the tree. Peaks bunched in the same region suggest one broad optimum.
- **`phytclust_results.tsv`** — has cluster assignments for each selected *k*, so you can compare them side by side.

!!! tip
    If `--max-k` is omitted, PhytClust defaults to scanning up to 90% of the leaf count. For large trees, setting `--max-k` explicitly can speed things up if you know the plausible range.

## 3. Multi-resolution: one *k* per scale

Sometimes you want to see the tree at multiple granularities — a coarse grouping, a medium one, and a fine one — rather than three neighboring peaks. Resolution mode splits the *k* range into logarithmic bins and picks the best *k* in each:

```bash
phytclust examples/sample_tree.nwk \
  --resolution \
  --bins 4 \
  --save-fig \
  --out-dir results/cli_resolution
```

This gives you 4 solutions spread across the *k* range, from broad to fine-grained. It's particularly useful for exploratory analysis or when presenting results to collaborators who want to see "the big picture and the details."

## 4. Constrain the clustering

Real trees produce messy clusters — tiny singletons, noise clades, biologically implausible groupings. PhytClust gives you two levers.

### Hard constraint: minimum cluster size

If you never want clusters smaller than 3 leaves:

```bash
phytclust examples/sample_tree.nwk \
  --top-n 3 \
  --min-cluster-size 3 \
  --out-dir results/cli_min_size
```

The DP enforces this during optimization, not as a post-filter. If a particular *k* can't be achieved without violating the constraint, it simply won't appear in the results.

### Soft control: outlier handling

A gentler approach: let small clusters exist but flag them as outliers. You can also tell the optimizer to prefer solutions with fewer outlier groups:

```bash
phytclust examples/sample_tree.nwk \
  --top-n 3 \
  --outlier-size-threshold 3 \
  --prefer-fewer-outliers \
  --out-dir results/cli_outlier
```

In the output TSV, outlier clusters are marked with `-1`. This is handy when you want to keep the information ("these three leaves are weird") without letting noise dominate the clustering.

## 5. Polytomies and zero-length edges

### Polytomies (multifurcations)

Many real trees have polytomies — internal nodes with more than two children. By default, PhytClust handles these natively using an optimized DP that reasons about all children simultaneously:

```bash
# default behavior: optimized polytomy handling
phytclust examples/sample_polytomy.newick --k 5
```

If you want the older behavior (inserting dummy nodes to make the tree binary first), use:

```bash
phytclust examples/sample_polytomy.newick --k 5 --no-optimize-polytomies
```

For fine control, you can choose between two polytomy modes:

- **`hard`** (default) — each child of a polytomy goes entirely into one cluster. Fast and predictable.
- **`soft`** — allows subsets of children to be merged across cluster boundaries. More flexible, but computationally expensive for high-degree nodes.

```bash
phytclust tree.nwk --k 5 --polytomy-mode soft --soft-polytomy-max-degree 15
```

!!! warning
    Soft mode is exponential in node degree. The `--soft-polytomy-max-degree` flag (default 18) is a safety guardrail — if a polytomy exceeds this degree, it falls back to hard mode for that node.

### Zero-length edges

Some trees have internal edges of length zero (often from collapsed uncertain nodes). By default PhytClust will split at these edges if the DP says so. If you'd rather keep zero-length clades together:

```bash
phytclust tree.nwk --k 5 --no-split-zero-length
```

This is useful for trees with many collapsed or near-identical leaves where splitting at zero-length branches would be arbitrary.

## 6. Use a config file for reproducibility

Once you've settled on good settings, save them in a YAML file so every run is reproducible:

```yaml
# my_config.yaml
peak:
  lambda_weight: 0.6
  ranking_mode: adjusted
  min_prominence: 5.0
  resolution_fallback_mode: max_score

runtime:
  plot:
    scores:
      log_scale_y: false
      fig_width: 20
    cluster:
      height_scale: 0.25
      marker_size: 55
  save:
    csv_name: phytclust_results.tsv
```

Then pass it with `--config`:

```bash
phytclust examples/sample_tree.nwk \
  --top-n 3 \
  --max-k 120 \
  --config my_config.yaml \
  --save-fig \
  --out-dir results/cli_config
```

CLI flags override config file values, so you can use the file as a baseline and tweak individual settings per run. See the [Configuration Reference](../reference/configuration.md) for all available fields.

## 7. Before you move on: a sanity checklist

Before feeding PhytClust results into downstream analysis, it's worth a quick check:

- [ ] **Score plot**: do the selected peaks look like real peaks, or are they in noisy/flat regions?
- [ ] **Peak positions**: are the selected *k* values biologically plausible? All peaks near `max-k` can indicate the tree has no strong structure at coarser scales.
- [ ] **Outlier clusters**: do the `-1` marked clusters match your expectation of "small weird clades"?
- [ ] **Tree plots**: are they legible? For large trees, adjust `height_scale` and `marker_size` in a config file.
- [ ] **Reproducibility**: is your run captured in a config file or script?
