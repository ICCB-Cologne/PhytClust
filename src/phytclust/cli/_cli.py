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
from ..algo.core import PhytClust
from pathlib import Path
from importlib import resources

try:
    from importlib.metadata import version, PackageNotFoundError
except Exception:
    try:
        from importlib_metadata import version, PackageNotFoundError
    except Exception:
        version = None
        PackageNotFoundError = Exception

try:
    from rich.console import Console
    from rich.status import Status
    from rich.text import Text
    from rich.panel import Panel
except Exception:
    Console = None
    Status = None
    Text = None
    Panel = None

console = Console() if Console else None

LOG = logging.getLogger("phytclust.cli")


def print_banner():
    try:
        text = (
            resources.files("phytclust")
            .joinpath("ascii_logo.txt")
            .read_text(encoding="utf-8")
        )
        print(text)
    except Exception:
        LOG.debug("ascii_logo.txt not found; skipping banner.")


def _positive_int(value: str) -> int:
    ivalue = int(value)
    if ivalue < 1:
        raise argparse.ArgumentTypeError("value must be ≥ 1")
    return ivalue


def _min_int(min_value: int):
    def _check(value: str) -> int:
        ivalue = int(value)
        if ivalue < min_value:
            raise argparse.ArgumentTypeError(f"value must be ≥ {min_value}")
        return ivalue

    return _check


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
    except PackageNotFoundError:
        return "unknown"


def _load_config(path: pathlib.Path | None) -> dict:
    if not path:
        return {}
    try:
        text = path.read_text()
        try:
            import yaml

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
        help="Write tsv rows for *every* k from 1..max_k (can be large!).",
    )
    sp.add_argument("--tsv-name", default="phytclust_results.tsv")
    sp.add_argument(
        "--no-tsv", action="store_true", help="Skip writing the results tsv."
    )
    sp.add_argument(
        "--dpi", type=_positive_int, default=150, help="PNG resolution (DPI)."
    )


def _configure_logging(
    verbosity: int, quiet: int, no_color: bool, log_format: str
) -> None:
    """
    Configure a single handler for the 'phytclust' logger hierarchy so that
    all submodules (phytclust.*, using logging.getLogger(__name__)) inherit it.
    """
    level = logging.WARNING
    if verbosity > 0:
        level = logging.INFO
    if verbosity > 1:
        level = logging.DEBUG
    if quiet > 0:
        level = logging.ERROR
    if quiet > 1:
        level = logging.CRITICAL

    root = logging.getLogger("phytclust")
    root.handlers.clear()
    root.setLevel(level)
    root.propagate = False

    handler = None
    if not no_color:
        try:
            from rich.logging import RichHandler

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
        LOG.info("⏱ %s: %.3fs", self.label, dt)


def build_parser() -> argparse.ArgumentParser:
    epilog = textwrap.dedent(
        """\
        examples:
          # exact k=5
          phytclust data/tree.nwk -k 5 --save-fig

          # global CalBow peaks: top 3, cap k at 200
          phytclust data/tree.nwk --top-n 3 --max-k 200 --save-fig --out-dir out

          # one peak per 4 log-bins, headless save
          phytclust data/tree.nwk --resolution --bins 4 --no-plot --save-all-k --save-fig
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
        "-o",
        "--out-dir",
        type=pathlib.Path,
        default=pathlib.Path("results"),
        help="Directory for PNGs / tsv (created if needed).",
    )
    p.add_argument("--outgroup", help="Taxon to exclude from all clusters.")
    p.add_argument(
        "--no-outlier",
        action="store_true",
        help="Do NOT mark singleton clusters as -1 in tsv.",
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
        "--time",
        action="store_true",
        help="Report elapsed time per phase and total.",
    )
    p.add_argument(
        "--progress",
        action="store_true",
        help="Show progress spinners for major phases (requires rich).",
    )
    p.add_argument(
        "--config",
        type=pathlib.Path,
        help="Optional JSON/YAML config with plotting/saving options.",
    )
    p.add_argument(
        "--version",
        action="version",
        version=f"phytclust { _package_version() }",
    )

    p.add_argument(
        "-k",
        "--k",
        type=_min_int(2),
        help="Exact k-way partition. If given, overrides --top-n/--resolution.",
    )
    p.add_argument(
        "--top-n",
        type=_positive_int,
        default=1,
        help="Number of global CalBow peaks to return (ignored if -k is set or --resolution is used).",
    )
    p.add_argument(
        "--resolution",
        action="store_true",
        help="Multi-resolution mode: one peak per log-spaced bin.",
    )
    p.add_argument(
        "--bins",
        type=_positive_int,
        default=3,
        help="Number of log bins when using --resolution.",
    )
    p.add_argument(
        "--max-k",
        type=_min_int(4),
        help="Upper bound for k in peak search (global or resolution modes).",
    )
    p.add_argument(
        "--max-k-limit",
        type=float,
        help=(
            "If --max-k is not given, set max_k = ceil(max_k_limit * num_leaves). "
            "Ignored when -k is specified."
        ),
    )
    p.add_argument(
        "--lambda-weight",
        type=float,
        dest="lambda_weight",
        default=0.5,
        help="Peak prominence parameter for score peak selection.",
    )

    _add_common_run_flags(p)

    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    _configure_logging(args.verbose, args.quiet, args.no_color, args.log_format)

    show_phase_ui = bool(console is not None and sys.stderr.isatty())

    print_banner()
    _emit_header(args)

    cfg = _load_config(args.config)
    plot_cfg = dict(cfg.get("plot") or {})
    save_cfg = dict(cfg.get("save") or {})

    t0_total = time.perf_counter()

    will_write_anything = bool(args.save_fig or args.save_tree or (not args.no_tsv))
    if will_write_anything:
        out_dir = Path(args.out_dir) if args.out_dir is not None else Path("results")
        out_dir.mkdir(parents=True, exist_ok=True)
        args.out_dir = str(out_dir)

    try:
        with phase(show_phase_ui, "load tree"):
            handle: Any = sys.stdin if args.tree == "-" else args.tree
            tree = Phylo.read(handle, "newick")
    except Exception as exc:
        LOG.error("Cannot read tree: %s", exc)
        sys.exit(2)

    try:
        with phase(show_phase_ui, "initialization"):
            pc = PhytClust(tree=tree, outgroup=args.outgroup)
        LOG.info("Tree terminals (after outgroup handling): %d", pc.num_terminals)
    except Exception as exc:
        LOG.error("Failed to initialize PhytClust: %s", exc)
        sys.exit(1)

    try:
        with phase(show_phase_ui, "clustering"):
            LOG.info(
                "Running PhytClust (k=%s, top_n=%d, resolution=%s, max_k=%s)…",
                args.k or "auto",
                args.top_n,
                args.resolution,
                args.max_k or "auto",
            )

            run_kwargs = dict(
                k=args.k,
                top_n=args.top_n,
                by_resolution=args.resolution,
                num_bins=args.bins,
                max_k=args.max_k,
                max_k_limit=args.max_k_limit,
                plot_scores=args.plot,
                lambda_weight=args.lambda_weight,
            )

            t0_run = time.perf_counter()
            result = pc.run(**run_kwargs)
            dt_run = time.perf_counter() - t0_run

            LOG.info("⏱ phytclust.run(): %.3fs", dt_run)

            ks = []
            if isinstance(result, dict):
                if "k" in result and result["k"] is not None:
                    ks = [int(result["k"])]
                elif "ks" in result and result["ks"] is not None:
                    ks = [int(x) for x in result["ks"]]

            if ks:
                LOG.info("Selected k: %s", ks if len(ks) > 1 else ks[0])

    except ValueError as e:
        LOG.error("Clustering failed: %s", e)
        sys.exit(1)
    except Exception as exc:
        LOG.error("Unexpected error during clustering: %s", exc)
        sys.exit(1)

    try:
        if args.plot or args.save_fig or args.save_tree:
            with phase(show_phase_ui, "render/plot"):
                pc.plot(
                    results_dir=(
                        args.out_dir if (args.save_fig or args.save_tree) else None
                    ),
                    save=(args.save_fig or args.save_tree),
                    n=args.k,
                    **plot_cfg,
                )
    except Exception as exc:
        LOG.warning("Plotting encountered a problem: %s", exc)

    if not args.no_tsv:
        try:
            with phase(show_phase_ui, "write tsv"):
                pc.save(
                    results_dir=args.out_dir,
                    filename=args.tsv_name,
                    outlier=not args.no_outlier,
                    output_all=args.save_all_k,
                    **save_cfg,
                )
        except Exception as exc:
            LOG.error("Failed to write tsv: %s", exc)
            sys.exit(1)

    if will_write_anything:
        _emit_artifacts(args.out_dir)

    LOG.info("Total runtime: %.3fs", time.perf_counter() - t0_total)

    LOG.info("Done.")
    return 0


def _emit_header(args: argparse.Namespace) -> None:
    if console is None or Panel is None or Text is None:
        LOG.info(
            "phytclust | k=%s | top_n=%s | resolution=%s | out_dir=%s",
            args.k or "auto",
            args.top_n,
            args.resolution,
            args.out_dir,
        )
        return

    mode = (
        "k" if args.k is not None else ("resolution" if args.resolution else "global")
    )
    body = Text()
    body.append("Mode: ", style="bold")
    body.append(f"{mode}\n")
    body.append("Tree: ", style="bold")
    body.append(f"{args.tree}\n")
    body.append("Output: ", style="bold")
    body.append(f"{args.out_dir}\n")
    body.append("Flags: ", style="bold")
    body.append(
        f"plot={args.plot} save_fig={args.save_fig} save_tree={args.save_tree} tsv={not args.no_tsv}"
    )

    console.print(Panel(body, title="phytclust", expand=False))


def _emit_artifacts(out_dir: str) -> None:
    try:
        p = Path(out_dir)
        if not p.exists():
            return
        files = sorted([x.name for x in p.iterdir() if x.is_file()])
        if not files:
            return

        if console is None:
            LOG.info("Wrote %d file(s) to %s", len(files), out_dir)
            return

        from rich.table import Table  # type: ignore

        t = Table(title="Output files", show_lines=False)
        t.add_column("File", style="cyan")
        for f in files:
            t.add_row(f)
        console.print(t)
    except Exception:
        return
