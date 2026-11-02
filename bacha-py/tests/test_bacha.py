"""Unit tests for bacha.tl.bacha (CLI invocation mocked)."""

from __future__ import annotations

import stat
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

import bacha


def _make_fake_bacha(tmp_path: Path, motifs: list[str]) -> Path:
    """Create a fake ``bacha`` binary that emits a deterministic scores file.

    The shell script ignores all flags and writes a CSV whose values are
    derived from the row/column position so we can verify the parser.
    """
    motifs_space = " ".join(motifs)
    script = tmp_path / "bacha"
    script.write_text(
        f"""#!/bin/bash
set -e
while [ "$#" -gt 0 ]; do
  case "$1" in
    -f) frags="$2"; shift 2;;
    -i) index="$2"; shift 2;;
    -c) cells="$2"; shift 2;;
    -o) out="$2"; shift 2;;
    -l|-p|-t) shift 2;;
    count) shift;;
    *) shift;;
  esac
done
case "$out" in
  *.gz) writer="gzip -c";;
  *)    writer="cat";;
esac
{{
  printf ','
  awk 'NF' "$cells" | paste -sd, -
  i=0
  for m in {motifs_space}; do
    printf '%s' "$m"
    n=$(awk 'NF' "$cells" | wc -l)
    j=0
    while [ "$j" -lt "$n" ]; do
      printf ',%d' "$((i * 10 + j + 1))"
      j=$((j + 1))
    done
    printf '\\n'
    i=$((i + 1))
  done
}} | $writer > "$out"
"""
    )
    script.chmod(script.stat().st_mode | stat.S_IEXEC | stat.S_IXUSR | stat.S_IXGRP)
    return script


@pytest.fixture
def fake_inputs(tmp_path: Path):
    motifs = ["MOTIF_A", "MOTIF_B"]
    binary = _make_fake_bacha(tmp_path, motifs)
    frags = tmp_path / "frags.tsv.gz"
    frags.write_bytes(b"")
    index = tmp_path / "fake.bi"
    index.write_bytes(b"")
    adata = ad.AnnData(
        X=np.zeros((3, 5), dtype=np.float32),
        obs=pd.DataFrame(index=["cell1", "cell2", "cell3"]),
    )
    return {
        "bacha_path": binary,
        "fragments": frags,
        "index": index,
        "adata": adata,
        "motifs": motifs,
    }


def test_bacha_returns_anndata(fake_inputs):
    result = bacha.count(
        fake_inputs["adata"],
        fragments=fake_inputs["fragments"],
        index=fake_inputs["index"],
        bacha_path=fake_inputs["bacha_path"],
        inplace=False,
        verbose=False,
    )
    assert isinstance(result, ad.AnnData)
    assert result.shape == (3, 2)
    assert list(result.obs_names) == ["cell1", "cell2", "cell3"]
    assert list(result.var_names) == fake_inputs["motifs"]


def test_bacha_inplace_writes_obsm(fake_inputs):
    adata = fake_inputs["adata"]
    out = bacha.count(
        adata,
        fragments=fake_inputs["fragments"],
        index=fake_inputs["index"],
        bacha_path=fake_inputs["bacha_path"],
        inplace=True,
        verbose=False,
    )
    assert out is None
    assert "bacha" in adata.obsm
    assert adata.obsm["bacha"].shape == (3, 2)
    assert list(adata.uns["bacha_motifs"]) == fake_inputs["motifs"]


def test_bacha_missing_fragments(fake_inputs, tmp_path):
    with pytest.raises(FileNotFoundError):
        bacha.count(
            fake_inputs["adata"],
            fragments=tmp_path / "nope.tsv.gz",
            index=fake_inputs["index"],
            bacha_path=fake_inputs["bacha_path"],
            inplace=False,
        )


def test_bacha_outfile_must_not_exist(fake_inputs, tmp_path):
    existing = tmp_path / "exists.csv.gz"
    existing.write_bytes(b"")
    with pytest.raises(FileExistsError):
        bacha.count(
            fake_inputs["adata"],
            fragments=fake_inputs["fragments"],
            index=fake_inputs["index"],
            bacha_path=fake_inputs["bacha_path"],
            outfile=existing,
            inplace=False,
        )


def test_bacha_inplace_requires_matching_obs(fake_inputs):
    with pytest.raises(ValueError, match="adata.obs_names"):
        bacha.count(
            fake_inputs["adata"],
            fragments=fake_inputs["fragments"],
            index=fake_inputs["index"],
            bacha_path=fake_inputs["bacha_path"],
            cells=["cell1", "cell2"],
            inplace=True,
            verbose=False,
        )
