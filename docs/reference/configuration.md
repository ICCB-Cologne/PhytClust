# Configuration Reference

PhytClust has two configuration layers:

1. `PeakConfig`: controls score-peak detection and ranking
2. `RuntimeConfig`: controls plotting and save defaults

## Python Usage

```python
from phytclust import (
    PeakConfig,
    RuntimeConfig,
    PlotConfig,
    ScorePlotConfig,
    ClusterPlotConfig,
)

peak_cfg = PeakConfig(lambda_weight=0.6, min_prominence=5.0)
runtime_cfg = RuntimeConfig(
    plot=PlotConfig(
        scores=ScorePlotConfig(log_scale_y=False, fig_width=22),
        cluster=ClusterPlotConfig(height_scale=0.25, marker_size=55),
    )
)
```

## YAML Usage (CLI)

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

Run with:

```bash
phytclust tree.nwk --top-n 3 --config examples/phytclust.config.yaml --out-dir out
```

## PeakConfig Fields

| Field | Type | Default | Notes |
| --- | --- | --- | --- |
| `lambda_weight` | `float` | `0.7` | Blend between raw and adjusted ranking terms. |
| `ranking_mode` | `str` | `"adjusted"` | `raw` or `adjusted`. |
| `boundary_window_size` | `int` | `5` | Right-side window for k=2 boundary candidate checks. |
| `boundary_ratio_threshold` | `float` | `1.5` | Minimum ratio for accepting k=2 boundary candidate. |
| `min_prominence` | `float | None` | `None` | If `None`, auto-derived from score range. |
| `min_k` | `int` | `2` | Ignore peaks below this k. |
| `resolution_fallback_mode` | `str` | `"none"` | `none` or `max_score`. |

## RuntimeConfig Shape

```text
RuntimeConfig
  plot: PlotConfig
    cluster: ClusterPlotConfig
    scores: ScorePlotConfig
  save: SaveConfig
```

### ScorePlotConfig

| Field | Default |
| --- | --- |
| `title_fontsize` | `50` |
| `axis_label_fontsize` | `35` |
| `tick_labelsize` | `30` |
| `peak_labelsize` | `50` |
| `bin_labelsize` | `30` |
| `fig_width` | `18` |
| `fig_height` | `10` |
| `clamp_negative_to_zero` | `True` |
| `log_scale_y` | `True` |
| `x_axis_mode` | `"log"` |
| `log_base` | `None` |
| `prefer_unsmoothed_primary` | `True` |
| `show_secondary_score_plot` | `False` |
| `colorblind_palette` | 8-color palette |

### ClusterPlotConfig

| Field | Default |
| --- | --- |
| `cmap` | `"tab20"` |
| `width_scale` | `2.5` |
| `height_scale` | `0.30` |
| `marker_size` | `40` |
| `show_branch_lengths` | `False` |
| `hide_internal_nodes` | `True` |

### SaveConfig

| Field | Default |
| --- | --- |
| `csv_name` | `"phytclust_results.csv"` |
| `outlier` | `True` |

## Backward-Compatible Config Blocks

CLI also accepts older top-level keys and maps them into runtime fields:

- `plot` -> `runtime.plot.cluster` subset
- `scores_plot` -> `runtime.plot.scores`
- `save` -> `runtime.save`

Prefer new `runtime` nesting for forward compatibility.
