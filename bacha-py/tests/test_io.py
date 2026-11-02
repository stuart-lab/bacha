"""Unit tests for bacha._io."""

from __future__ import annotations

import numpy as np
import pytest

from bacha import _io


def test_read_scores_round_trip(small_scores_csv):
    path, values, motifs, cells = small_scores_csv
    scores, motif_ids, cell_ids = _io.read_scores(path)
    assert motif_ids == motifs
    assert cell_ids == cells
    np.testing.assert_allclose(scores, values, atol=1e-5)


def test_align_cells_reorders(small_scores_csv):
    _, values, _, cells = small_scores_csv
    target = [cells[2], cells[0], cells[3], cells[1]]
    out = _io.align_cells(values, cells, target)
    np.testing.assert_array_equal(out[:, 0], values[:, 2])
    np.testing.assert_array_equal(out[:, 1], values[:, 0])
    np.testing.assert_array_equal(out[:, 3], values[:, 1])


def test_align_cells_missing(small_scores_csv):
    _, values, _, cells = small_scores_csv
    with pytest.raises(ValueError, match="missing"):
        _io.align_cells(values, cells, [*cells, "ghost_cell"])


def test_align_cells_identity_short_circuit(small_scores_csv):
    _, values, _, cells = small_scores_csv
    out = _io.align_cells(values, cells, list(cells))
    assert out is values
