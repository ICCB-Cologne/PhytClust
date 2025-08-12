import subprocess
import pathlib
import sys

TREE_FILE = pathlib.Path(__file__).parent / "test_tree.nwk"


def run_cli(args):
    """Run the CLI as if from terminal."""
    cmd = [sys.executable, "-m", "phytclust.cli"] + args
    return subprocess.run(cmd, capture_output=True, text=True)


def test_cli_k_mode(tmp_path):
    out_dir = tmp_path / "results"
    result = run_cli(
        [str(TREE_FILE), "k", "--k", "2", "--save-fig", "--out-dir", str(out_dir)]
    )
    assert result.returncode == 0
    assert out_dir.exists()


def test_cli_resolution_mode(tmp_path):
    out_dir = tmp_path / "results"
    result = run_cli(
        [
            str(TREE_FILE),
            "resolution",
            "--bins",
            "2",
            "--save-fig",
            "--out-dir",
            str(out_dir),
        ]
    )
    assert result.returncode == 0
