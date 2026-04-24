# Configuration Reference

PhytClust settings are split into two groups: **PeakConfig** (how peaks are detected and ranked) and **RuntimeConfig** (how plots look and where output goes). You can set these in Python or in a YAML config file.

---

## Quick example

=== "Python"

    ```python
    from phytclust import (
        PeakConfig, RuntimeConfig, PlotConfig,
        ScorePlotConfig, ClusterPlotConfig,
    )

    peak_cfg = PeakConfig(lambda_weight=0.6, min_prominence=5.0)

    runtime_cfg = RuntimeConfig(
        plot=PlotConfig(
            scores=ScorePlotConfig(log_scale_y=False, fig_width=22),
            cluster=ClusterPlotConfig(height_scale=0.25, marker_size=55),
        )
    )
    ```

=== "YAML"

    ```yaml
    peak:
      lambda_weight: 0.6
      min_prominence: 5.0
      ranking_mode: adjusted

    runtime:
      plot:
        scores:
          log_scale_y: false
          fig_width: 22
        cluster:
          height_scale: 0.25
          marker_size: 55
      save:
        csv_name: phytclust_results.tsv
        outlier: true
    ```

    Pass to the CLI with `--config`:

    ```bash
    phytclust tree.nwk --top-n 3 --config my_config.yaml
    ```

CLI flags always override config file values.

---

## PeakConfig fields

These control how PhytClust finds and ranks peaks in the score-vs-*k* curve.

| Field | Type | Default | What it does |
|-------|------|---------|-------------|
| `lambda_weight` | float | `0.7` | Blends raw and adjusted ranking. 0 = pure raw score, 1 = pure adjusted score. Higher values give more weight to outlier-corrected rankings. |
| `ranking_mode` | str | `"adjusted"` | `"raw"` uses the score directly. `"adjusted"` applies outlier correction before ranking. |
| `boundary_window_size` | int | `5` | Window size for evaluating whether *k*=2 is a valid peak (compared against right-side neighbors). |
| `boundary_ratio_threshold` | float | `1.5` | Minimum score ratio for *k*=2 to qualify as a peak vs. its right-window median. |
| `min_prominence` | float or None | `None` | Minimum peak prominence. `None` means auto-derived from the score range — usually a good default. |
| `use_log_peak_input` | bool | `False` | Detect peaks on log-transformed scores instead of raw values. |
| `use_relative_prominence` | bool | `False` | Rank peaks by relative prominence (fold-change) instead of absolute. |
| `prominence_k_power` | float | `0.0` | Scale the prominence threshold by *k*^power. At 0 the threshold is constant. |
| `min_k` | int | `2` | Ignore any peak below this *k*. |
| `resolution_fallback_mode` | str | `"none"` | What to do when a resolution bin has no peak: `"none"` skips it, `"max_score"` uses the *k* with the highest score in that bin. |

## RuntimeConfig structure

RuntimeConfig nests plotting and save settings:

```
RuntimeConfig
  plot: PlotConfig
    scores: ScorePlotConfig    ← score curve appearance
    cluster: ClusterPlotConfig ← tree plot appearance
  save: SaveConfig             ← output file settings
```

### ScorePlotConfig

Controls the score-vs-*k* curve figure.

| Field | Default | What it does |
|-------|---------|-------------|
| `title_fontsize` | `50` | Title font size |
| `axis_label_fontsize` | `35` | Axis label font size |
| `tick_labelsize` | `30` | Tick label font size |
| `peak_labelsize` | `50` | Font size for peak annotations |
| `bin_labelsize` | `30` | Font size for resolution bin labels |
| `fig_width` | `18` | Figure width in inches |
| `fig_height` | `10` | Figure height in inches |
| `clamp_negative_to_zero` | `True` | Replace negative scores with zero in the plot |
| `log_scale_y` | `True` | Use log scale for the Y axis |
| `x_axis_mode` | `"log"` | `"log"` or `"linear"` for the X axis |
| `log_base` | `None` | Base for log scale (None = matplotlib default) |
| `prefer_unsmoothed_primary` | `True` | Show raw scores as the primary line (smoothed as secondary) |
| `show_secondary_score_plot` | `False` | Add a second panel with the alternative score view |
| `colorblind_palette` | 8 colors | Palette used for peak and bin markers (colorblind-safe by default) |

### ClusterPlotConfig

Controls the colored tree figures.

| Field | Default | What it does |
|-------|---------|-------------|
| `cmap` | `"phytclust"` | Matplotlib colormap name. `"phytclust"` uses the project's brand palette. |
| `width_scale` | `2.5` | Horizontal stretch factor for the tree |
| `height_scale` | `0.30` | Vertical stretch factor (increase for large trees) |
| `marker_size` | `40` | Size of leaf markers |
| `show_branch_lengths` | `False` | Annotate branches with their lengths |
| `hide_internal_nodes` | `True` | Hide internal node labels for a cleaner look |

### SaveConfig

Controls output file behavior.

| Field | Default | What it does |
|-------|---------|-------------|
| `csv_name` | `"phytclust_results.csv"` | Default filename for the cluster assignment table |
| `outlier` | `True` | Mark outlier clusters as `-1` in output |

---

## Backward-compatible config blocks

Older config files may use flat top-level keys. These still work but the nested `runtime` structure is preferred:

| Old key | Maps to |
|---------|---------|
| `plot` | `runtime.plot.cluster` (subset) |
| `scores_plot` | `runtime.plot.scores` |
| `save` | `runtime.save` |

---

## Recipes

Here are a few ready-made configurations for common scenarios.

### Conservative (trust strong signals only)

```yaml
peak:
  lambda_weight: 0.8
  min_prominence: 10.0
  min_k: 3
```

### Exploratory (catch weaker peaks too)

```yaml
peak:
  lambda_weight: 0.4
  min_prominence: 1.0
  use_relative_prominence: true
```

### Publication-ready figures

```yaml
runtime:
  plot:
    scores:
      fig_width: 24
      fig_height: 14
      title_fontsize: 60
      axis_label_fontsize: 40
      peak_labelsize: 55
    cluster:
      height_scale: 0.20
      marker_size: 35
      hide_internal_nodes: true
      show_branch_lengths: false
```

### Large tree (500+ leaves)

```yaml
runtime:
  plot:
    cluster:
      height_scale: 0.10
      marker_size: 15
      width_scale: 3.0
    scores:
      x_axis_mode: log
      log_scale_y: true
```
