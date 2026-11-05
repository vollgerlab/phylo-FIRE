"""Shared pytest setup: import the pipeline scripts as top-level modules.

Run with ``pixi run unit``. Tests that need the AGC archive request the
``agc_archive`` fixture and skip when it (or the ``agc`` binary) is missing.
Set ``PHYLO_FIRE_AGC`` to point at another archive, or at a missing path to
force those skips.
"""
import os
import shutil
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1]
REPO = SCRIPTS.parents[1]
sys.path[:0] = [str(SCRIPTS), str(SCRIPTS / "tools")]

DEFAULT_AGC = REPO / "resources" / "HPRC_r2_assemblies_0.6.1.agc"


@pytest.fixture(scope="session")
def agc_archive() -> str:
    """Path to the HPRC AGC archive; skips the test when it is unavailable."""
    agc = Path(os.environ.get("PHYLO_FIRE_AGC", DEFAULT_AGC))
    if not agc.is_file():
        pytest.skip(f"AGC archive not found: {agc}")
    if shutil.which("agc") is None:
        pytest.skip("agc binary not on PATH")
    return str(agc)
