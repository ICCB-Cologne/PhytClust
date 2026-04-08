import subprocess
import pathlib
import sys

TREE_FILE = pathlib.Path(__file__).parent / "test_tree.nwk"


def run_cli(args):
    cmd = [sys.executable, "-m", "phytclust.cli"] + args
    return subprocess.run(cmd, capture_output=True, text=True)


def _assert_ok(result):
    if result.returncode != 0:
        raise AssertionError(
            "CLI failed\n"
            f"returncode={result.returncode}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}\n"
        )


def test_cli_k_mode(tmp_path):
    out_dir = tmp_path / "results"
    result = run_cli(
        [str(TREE_FILE), "--k", "2", "--save-fig", "--out-dir", str(out_dir)]
    )

    _assert_ok(result)
    assert out_dir.exists()
    assert any(out_dir.iterdir())


def test_cli_resolution_mode(tmp_path):
    out_dir = tmp_path / "results"
    result = run_cli(
        [
            str(TREE_FILE),
            "--bins",
            "2",
            "--resolution",
            "--save-fig",
            "--out-dir",
            str(out_dir),
        ]
    )

    _assert_ok(result)
    assert out_dir.exists()
    assert any(out_dir.iterdir())


def test_cli_soft_polytomy_mode(tmp_path):
    out_dir = tmp_path / "results"
    result = run_cli(
        [
            str(TREE_FILE),
            "--k",
            "2",
            "--polytomy-mode",
            "soft",
            "--soft-polytomy-max-degree",
            "8",
            "--save-fig",
            "--out-dir",
            str(out_dir),
        ]
    )

    _assert_ok(result)
    assert out_dir.exists()
    assert any(out_dir.iterdir())
