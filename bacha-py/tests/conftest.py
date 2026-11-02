"""Shared pytest fixtures."""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pytest


@pytest.fixture
def small_scores_csv(tmp_path: Path) -> tuple[Path, np.ndarray, list[str], list[str]]:
    """A tiny gzipped bacha-style CSV: 3 motifs x 4 cells."""
    motifs = ["MOTIF_A", "MOTIF_B", "MOTIF_C"]
    cells = ["cellA", "cellB", "cellC", "cellD"]
    values = np.arange(12, dtype=np.float32).reshape(3, 4)
    path = tmp_path / "scores.csv.gz"
    with gzip.open(path, "wt") as fh:
        fh.write("," + ",".join(cells) + "\n")
        for m, row in zip(motifs, values):
            fh.write(m + "," + ",".join(f"{v:.6f}" for v in row) + "\n")
    return path, values, motifs, cells
