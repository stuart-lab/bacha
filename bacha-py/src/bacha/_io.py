"""Parse the CSV output of `bacha count` into numpy arrays."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def read_scores(path: str | Path) -> tuple[np.ndarray, list[str], list[str]]:
    """Read a `bacha count` CSV/CSV.GZ into a (motifs × cells) score matrix.

    Returns
    -------
    scores
        Float numpy array of shape (n_motifs, n_cells).
    motif_ids
        Row labels (motif / feature names).
    cell_ids
        Column labels (cell barcodes).
    """
    df = pd.read_csv(path, index_col=0)
    motif_ids = list(map(str, df.index))
    cell_ids = list(map(str, df.columns))
    scores = df.to_numpy(dtype=np.float32, copy=False)
    return scores, motif_ids, cell_ids


def align_cells(
    scores: np.ndarray,
    cell_ids: list[str],
    target_cells: list[str],
) -> np.ndarray:
    """Reorder columns of ``scores`` to match ``target_cells``.

    Raises if any target cell is missing from ``cell_ids``.
    """
    if cell_ids == target_cells:
        return scores
    pos = {c: i for i, c in enumerate(cell_ids)}
    missing = [c for c in target_cells if c not in pos]
    if missing:
        raise ValueError(
            f"{len(missing)} requested cells missing from bacha output "
            f"(first few: {missing[:5]})"
        )
    idx = np.fromiter((pos[c] for c in target_cells), dtype=np.int64,
                       count=len(target_cells))
    return scores[:, idx]
