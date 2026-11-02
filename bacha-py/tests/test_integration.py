"""End-to-end integration tests against the real ``bacha`` binary.

Unlike the rest of the suite (which mocks the CLI), these run the actual
``bacha`` binary on the small chr19 data bundled in the repo's ``test/``
directory. They are skipped automatically when the binary or that data is not
available, so a plain ``pip install`` checkout still collects green.

To run them, build the binary first (``cargo build --release`` in the repo
root) or install it (``cargo install --path .``).
"""

from __future__ import annotations

import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

import bacha

# Repo layout: this file is at bacha_code/bacha-py/tests/.
REPO = Path(__file__).resolve().parents[2]
TEST_DATA = REPO / "test"
INDEX = TEST_DATA / "hg38_chr19_test.bi"
FRAGMENTS = TEST_DATA / "GSM7673557_chr19_fragments.tsv.gz"
CELLS = TEST_DATA / "cells_subset.txt"

# Prefer an installed `bacha`; fall back to the local release build.
_release_build = REPO / "target" / "release" / "bacha"
_on_path = shutil.which("bacha")
BACHA_BIN: Path | None = (
    Path(_on_path) if _on_path else (_release_build if _release_build.is_file() else None)
)

_data_present = INDEX.is_file() and FRAGMENTS.is_file() and CELLS.is_file()

pytestmark = [
    pytest.mark.skipif(BACHA_BIN is None, reason="bacha binary not found (build with `cargo build --release`)"),
    pytest.mark.skipif(not _data_present, reason=f"bundled test data not found under {TEST_DATA}"),
]


@pytest.fixture(scope="module")
def barcodes() -> list[str]:
    return [ln.strip() for ln in CELLS.read_text().splitlines() if ln.strip()]


@pytest.fixture
def adata(barcodes) -> ad.AnnData:
    return ad.AnnData(
        X=np.zeros((len(barcodes), 1), dtype=np.float32),
        obs=pd.DataFrame(index=barcodes),
    )


def test_summary_runs(barcodes):
    info = bacha.summary(INDEX, bacha_path=BACHA_BIN)
    assert info["format"] == 1
    assert info["genome"] == "BSgenome.Hsapiens.UCSC.hg38"
    assert info["motif_set"] == "JASPAR2024"
    assert info["window"] == 50
    assert info["chromosomes"] == [("chr19", 58617616)]
    assert len(info["sites"]) == 5
    assert len(info["controls"]) == 5
    assert [s[0] for s in info["sites"]] == [c[0] for c in info["controls"]]


def test_count_not_inplace(adata, barcodes):
    result = bacha.count(
        adata, fragments=FRAGMENTS, index=INDEX,
        bacha_path=BACHA_BIN, inplace=False, verbose=False,
    )
    assert isinstance(result, ad.AnnData)
    assert result.shape == (len(barcodes), 5)
    assert list(result.obs_names) == barcodes
    assert not np.isnan(result.X).any()


def test_count_inplace_matches_not_inplace(adata):
    ref = bacha.count(
        adata, fragments=FRAGMENTS, index=INDEX,
        bacha_path=BACHA_BIN, inplace=False, verbose=False,
    )
    out = bacha.count(
        adata, fragments=FRAGMENTS, index=INDEX,
        bacha_path=BACHA_BIN, inplace=True, verbose=False,
    )
    assert out is None
    assert adata.obsm["bacha"].shape == ref.shape
    assert len(adata.uns["bacha_motifs"]) == ref.n_vars
    assert list(adata.uns["bacha_motifs"]) == list(ref.var_names)
    np.testing.assert_allclose(adata.obsm["bacha"], ref.X, equal_nan=True)
