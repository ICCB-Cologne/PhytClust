#!/usr/bin/env python3
from __future__ import annotations

import argparse
import pathlib
import sys
import textwrap
import logging
import time
from typing import Any

from Bio import Phylo
from phytclust import PhytClust

try:
    from importlib.metadata import version, PackageNotFoundError
except Exception:
    try:
        from importlib_metadata import version, PackageNotFoundError  # type: ignore
    except Exception:  # pragma: no cover
        version = None
        PackageNotFoundError = Exception  # type: ignore
# optional rich UI
try:
    from rich.console import Console
    from rich.status import Status
except Exception:
    Console = None
    Status = None

console = Console() if Console else None

LOG = logging.getLogger("phytclust.cli")


# ───────────────────────── helpers ─────────────────────────


def _positive_int(value: str) -> int:
    ivalue = int(value)
    if ivalue < 1:
        raise argparse.ArgumentTypeError("value must be ≥ 1")
    return ivalue


def _existing_path_or_stdin(p: str) -> pathlib.Path | str:
    if p == "-":
        return p
    path = pathlib.Path(p)
    if not path.exists():
        raise argparse.ArgumentTypeError(f"file not found: {p}")
    return path


def _package_version() -> str:
    if version is None:
        return "unknown"
    try:
        return version("phytclust")
    except PackageNotFoundError:  # pragma: no cover
        return "unknown"


def _load_config(path: pathlib.Path | None) -> dict:
    if not path:
        return {}
    try:
        text = path.read_text()
        # try YAML first if available; fallback to JSON
        try:
            import yaml  # type: ignore

            data = yaml.safe_load(text)
            return data or {}
        except Exception:
            import json

            return json.loads(text or "{}")
    except Exception as exc:
        LOG.warning("Could not read config %s: %s", path, exc)
        return {}


def _add_common_run_flags(sp: argparse.ArgumentParser) -> None:
    sp.add_argument(
        "--plot/--no-plot",
        dest="plot",
        default=True,
        action=argparse.BooleanOptionalAction,
        help="Display interactive plots after running.",
    )
    sp.add_argument(
        "--save-fig",
        action="store_true",
        help="Write scores.png and coloured tree PNG(s) to --out-dir.",
    )
    sp.add_argument(
        "--save-tree",
        action="store_true",
        help="Explicitly write coloured tree PNG(s) (in addition to --save-fig).",
    )
    sp.add_argument(
        "--save-all-k",
        action="store_true",
        help="Write CSV rows for *every* k from 1..max_k (can be large!).",
    )
    sp.add_argument("--csv-name", default="phytclust_results.csv")
    sp.add_argument(
        "--no-csv", action="store_true", help="Skip writing the results CSV."
    )
    sp.add_argument(
        "--dpi", type=_positive_int, default=150, help="PNG resolution (DPI)."
    )
    sp.add_argument(
        "-o", "--out-dir",
        type=pathlib.Path,
        help="Where to write PNGs / CSV. If omitted here, top-level --out-dir is used.",
    )


def _configure_logging(
    verbosity: int, quiet: int, no_color: bool, log_format: str
) -> None:
    """
    Configure a single handler for the 'phytclust' logger hierarchy so that
    all submodules (phytclust.*, using logging.getLogger(__name__)) inherit it.
    """
    # Decide level from -v / -q
    level = logging.WARNING
    if verbosity > 0:
        level = logging.INFO
    if verbosity > 1:
        level = logging.DEBUG
    if quiet > 0:
        level = logging.ERROR
    if quiet > 1:
        level = logging.CRITICAL

    # Remove any existing handlers (avoid duplicate logs during tests)
    root = logging.getLogger("phytclust")
    root.handlers.clear()
    root.setLevel(level)
    root.propagate = False

    # Try rich (nicer formatting/colors) unless disabled
    handler = None
    if not no_color:
        try:
            from rich.logging import RichHandler  # type: ignore

            handler = RichHandler(
                level=level,
                markup=True,
                rich_tracebacks=True,
                show_time=True,
                show_level=True,
                show_path=False,
            )
            formatter = logging.Formatter(log_format)
            handler.setFormatter(formatter)
        except Exception:
            handler = None

    if handler is None:
        # Plain stdlib handler
        handler = logging.StreamHandler()
        formatter = logging.Formatter(log_format)
        handler.setFormatter(formatter)

    root.addHandler(handler)


class phase:
    """
    Spinner + timing in one place. Shows a 'rich' spinner if --progress is set
    and rich is installed. Always logs elapsed time at INFO.
    """

    def __init__(self, enabled: bool, label: str):
        self.enabled = bool(enabled and console is not None)
        self.label = label
        self.t0 = None
        self.status = None

    def __enter__(self):
        import time

        self.t0 = time.perf_counter()
        if self.enabled:
            self.status = console.status(f"[bold]{self.label}…", spinner="dots")
            self.status.start()
        return self

    def __exit__(self, exc_type, exc, tb):
        import time

        dt = time.perf_counter() - self.t0
        if self.enabled and self.status:
            self.status.stop()
        # single source of truth for timing lines:
        LOG.info("⏱ %s: %.3fs", self.label, dt)


# def phase(enabled: bool, label: str):
#     """
#     Context manager to time phases when --time is set.
#     Usage: with phase(args.time, "load tree"):
#     """

#     class _Ctx:
#         def __enter__(self):
#             self.t0 = time.perf_counter()
#             return self

#         def __exit__(self, exc_type, exc, tb):
#             if enabled:
#                 dt = time.perf_counter() - self.t0
#                 LOG.info("⏱ %s: %.3fs", label, dt)

#     return _Ctx()


def build_parser() -> argparse.ArgumentParser:
    epilog = textwrap.dedent(
        """\
        examples:
          # exact k=5
          phytclust data/tree.nwk k --k 5 --save-fig

          # global CalBow peaks: top 3, cap k at 200
          phytclust data/tree.nwk global --top-n 3 --max-k 200 --save-fig --out-dir out

          # one peak per 4 log-bins, headless save
          phytclust - resolution --bins 4 --no-plot --save-all-k --save-fig
        """
    )

    p = argparse.ArgumentParser(
        prog="phytclust",
        description="Monophyletic clustering of phylogenetic trees (dynamic programming).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=epilog,
    )

    p.add_argument(
        "tree", type=_existing_path_or_stdin, help="Newick file or '-' for stdin."
    )
    p.add_argument(
        "-o", "--out-dir", type=pathlib.Path, default=pathlib.Path("results")
    )
    p.add_argument("--outgroup", help="Taxon to exclude from all clusters.")
    p.add_argument(
        "--no-outlier",
        action="store_true",
        help="Do NOT mark singleton clusters as –1 in CSV.",
    )
    p.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase logging verbosity (-vv for DEBUG).",
    )
    p.add_argument(
        "-q",
        "--quiet",
        action="count",
        default=0,
        help="Reduce logging verbosity (-qq for CRITICAL).",
    )
    p.add_argument("--no-color", action="store_true", help="Disable colorized logging.")
    p.add_argument(
        "--log-format",
        default="%(message)s",
        help="Logging format (plain logging.Formatter style).",
    )
    p.add_argument(
        "--time", action="store_true", help="Report elapsed time per phase and total."
    )
    p.add_argument("--progress", action="store_true",
               help="Show progress spinners for major phases (requires rich).")
    p.add_argument("--config", type=pathlib.Path,
               help="Optional JSON/YAML config with plotting/saving options.")

    p.add_argument(
        "--version", action="version", version=f"phytclust { _package_version() }"
    )

    sub = p.add_subparsers(dest="mode", required=True, help="Choose a run-mode.")

    s_k = sub.add_parser("k", help="Exact k-way partition")
    s_k.add_argument(
        "--k", type=_positive_int, required=True, help="Number of clusters."
    )
    _add_common_run_flags(s_k)

    s_g = sub.add_parser("global", help="Global CalBow peak search")
    s_g.add_argument(
        "--top-n", type=_positive_int, default=1, help="How many top peaks to return."
    )
    s_g.add_argument(
        "--max-k", type=_positive_int, help="Upper bound for k search space."
    )
    _add_common_run_flags(s_g)

    s_r = sub.add_parser("resolution", help="One peak per log-spaced bin")
    s_r.add_argument(
        "--bins", type=_positive_int, default=3, help="Number of log bins."
    )
    s_r.add_argument(
        "--max-k", type=_positive_int, help="Upper bound for k search space."
    )
    _add_common_run_flags(s_r)

    return p


# ───────────────────────── main ─────────────────────────


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    _configure_logging(args.verbose, args.quiet, args.no_color, args.log_format)
    if getattr(args, "out_dir", None) is None:
        pass
    # Load config (YAML or JSON)
    cfg = _load_config(args.config)  # expects args.config to be a path string
    plot_cfg = dict(cfg.get("plot") or {})
    save_cfg = dict(cfg.get("save") or {})

    t0_total = time.perf_counter()

    # Prepare output dir early if we might write
    if args.save_fig or args.save_tree or not args.no_csv:
        try:
            args.out_dir.mkdir(parents=True, exist_ok=True)
        except Exception as exc:
            LOG.error("Cannot create output directory '%s': %s", args.out_dir, exc)
            sys.exit(2)

    # Load tree
    try:
        with phase(args.time, "load tree"):
            handle: Any = sys.stdin if args.tree == "-" else args.tree
            tree = Phylo.read(handle, "newick")
    except Exception as exc:
        LOG.error("Cannot read tree: %s", exc)
        sys.exit(2)

    try:
        with phase(args.time, "initialization"):
            pc = PhytClust(tree=tree, outgroup=args.outgroup)
        LOG.info("Tree terminals: %d", pc.num_terminals)
    except Exception as exc:
        LOG.error("Failed to initialize PhytClust: %s", exc)
        sys.exit(1)

    # Run selected mode
    try:
        if args.mode == "k":
            with phase(args.time, "compute k-clusters"):
                LOG.info("Exact partition (k=%d)…", args.k)
                try:
                    pc.get_clusters(args.k)
                except ValueError as e:
                    # Common case: "Cannot form k clusters on this tree."
                    LOG.warning("No valid k=%d partition: %s", args.k, e)

        elif args.mode == "global":
            with phase(args.time, "global peak search"):
                LOG.info(
                    "Global peak search (top_n=%d, max_k=%s)…",
                    args.top_n,
                    args.max_k or "auto",
                )
                pc.best_global(
                    top_n=args.top_n, max_k=args.max_k, plot_scores=args.plot
                )

        elif args.mode == "resolution":
            with phase(args.time, "resolution peak search"):
                LOG.info(
                    "Resolution mode (bins=%d, max_k=%s)…",
                    args.bins,
                    args.max_k or "auto",
                )
                pc.best_by_resolution(
                    num_bins=args.bins, max_k=args.max_k, plot_scores=args.plot
                )

    except Exception as exc:
        LOG.error("Clustering failed: %s", exc)
        sys.exit(1)

    # Plot
    try:
        if args.plot or args.save_fig or args.save_tree:
            with phase(args.time, "render/plot"):
                pc.plot(
                    results_dir=(
                        args.out_dir if (args.save_fig or args.save_tree) else None
                    ),
                    save=(args.save_fig or args.save_tree),
                    **plot_cfg,  # load defaults from config
                    outlier=not args.no_outlier,
                )
    except Exception as exc:
        LOG.warning("Plotting encountered a problem: %s", exc)

    # CSV
    if not args.no_csv:
        try:
            with phase(args.time, "write CSV"):
                pc.save(
                    results_dir=args.out_dir,
                    **save_cfg,  # load defaults from config
                    outlier=not args.no_outlier,
                    output_all=args.save_all_k,
                )
        except Exception as exc:
            LOG.error("Failed to write CSV: %s", exc)
            sys.exit(1)

    if args.time:
        LOG.info("Total runtime: %.3fs", time.perf_counter() - t0_total)

    LOG.info("Done.")
