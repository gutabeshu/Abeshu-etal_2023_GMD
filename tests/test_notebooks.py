"""
Smoke tests for the Abeshu-etal_2023_GMD figure notebooks.

Run with:  pytest tests/
"""
import json
import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).parent.parent
FIGURES_DIR = REPO_ROOT / "figures"
NOTEBOOKS = sorted(FIGURES_DIR.glob("*.ipynb"))


def _config_keys(path):
    keys = set()
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or ":" not in line:
            continue
        keys.add(line.split(":", 1)[0].strip())
    return keys

# ---------------------------------------------------------------------------
# Config file tests
# ---------------------------------------------------------------------------

def test_config_yml_exists():
    assert (FIGURES_DIR / "config.yml").exists(), "figures/config.yml is missing"


def test_config_py_exists():
    assert (FIGURES_DIR / "config.py").exists(), "figures/config.py is missing"


def test_config_yml_required_keys():
    cfg = _config_keys(FIGURES_DIR / "config.yml")
    for key in ("dir_in", "dir_out"):
        assert key in cfg, f"config.yml is missing required key: '{key}'"


def test_config_py_valid_syntax():
    source = (FIGURES_DIR / "config.py").read_text()
    compile(source, "config.py", "exec")


# ---------------------------------------------------------------------------
# Notebook structural tests
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("nb_path", NOTEBOOKS, ids=lambda p: p.name)
def test_notebook_is_valid_json(nb_path):
    with open(nb_path) as f:
        nb = json.load(f)
    assert "cells" in nb
    assert "nbformat" in nb


@pytest.mark.parametrize("nb_path", NOTEBOOKS, ids=lambda p: p.name)
def test_notebook_has_no_embedded_outputs(nb_path):
    with open(nb_path) as f:
        nb = json.load(f)
    for cell in nb["cells"]:
        if cell["cell_type"] == "code":
            assert cell.get("outputs", []) == [], (
                f"{nb_path.name} has embedded cell outputs — run nbstripout"
            )


@pytest.mark.parametrize("nb_path", NOTEBOOKS, ids=lambda p: p.name)
def test_notebook_imports_config(nb_path):
    """Every notebook must import from config rather than hardcode paths."""
    with open(nb_path) as f:
        nb = json.load(f)
    full_source = "\n".join(
        "".join(cell["source"])
        for cell in nb["cells"]
        if cell["cell_type"] == "code"
    )
    assert "from config import" in full_source, (
        f"{nb_path.name} does not import from config.py"
    )


# ---------------------------------------------------------------------------
# Hardcoded-path guard
# ---------------------------------------------------------------------------

HARDCODED_PATH_PATTERNS = [
    r"D:/XanthosDev",
    r"D:\\\\XanthosDev",
    r"gwabeshu\.COUGARNET",
    r"D:/Xanthos-Repo",
    r"D:\\\\Xanthos-Repo",
    r"D:/RunningYalingsData",
    r"D:\\\\RunningYalingsData",
]


@pytest.mark.parametrize("nb_path", NOTEBOOKS, ids=lambda p: p.name)
def test_notebook_has_no_hardcoded_windows_paths(nb_path):
    with open(nb_path) as f:
        nb = json.load(f)
    violations = []
    for cell in nb["cells"]:
        if cell["cell_type"] != "code":
            continue
        for line in cell.get("source", []):
            if line.strip().startswith("#"):
                continue
            for pat in HARDCODED_PATH_PATTERNS:
                if re.search(pat, line):
                    violations.append(line.strip()[:100])
                    break
    assert not violations, (
        f"{nb_path.name} contains hardcoded Windows paths:\n"
        + "\n".join(violations)
    )
