"""Unit tests for bacha._cli."""

from __future__ import annotations

import os
import stat
from pathlib import Path

import pytest

from bacha import _cli


def test_build_count_args_matches_r_package():
    args = _cli.build_count_args(
        fragments="/data/frags.tsv.gz",
        index="/idx/hg38.bi",
        cells_file="/tmp/cells.txt",
        outfile="/tmp/out.csv.gz",
        lambda_=10.0,
        pseudocount=1.0,
        threads=4,
    )
    assert args == [
        "count",
        "-f", "/data/frags.tsv.gz",
        "-i", "/idx/hg38.bi",
        "-c", "/tmp/cells.txt",
        "-o", "/tmp/out.csv.gz",
        "-l", "10.0",
        "-p", "1.0",
        "-t", "4",
    ]


def test_locate_explicit_path(tmp_path: Path):
    fake = tmp_path / "bacha"
    fake.write_text("#!/bin/sh\nexit 0\n")
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
    assert _cli.locate(fake) == fake.resolve()


def test_locate_missing_explicit_path(tmp_path: Path):
    with pytest.raises(_cli.BachaNotFoundError):
        _cli.locate(tmp_path / "does-not-exist")


def test_locate_not_on_path(monkeypatch):
    monkeypatch.setenv("PATH", "/nonexistent-path-xyz")
    with pytest.raises(_cli.BachaNotFoundError):
        _cli.locate(None)


def test_run_propagates_nonzero(tmp_path: Path):
    failing = tmp_path / "fail.sh"
    failing.write_text("#!/bin/sh\necho oops 1>&2\nexit 7\n")
    failing.chmod(failing.stat().st_mode | stat.S_IEXEC | stat.S_IXUSR)
    with pytest.raises(_cli.BachaCLIError) as exc:
        _cli.run(failing, ["whatever"], verbose=False)
    assert "code 7" in str(exc.value)


def test_build_count_args_manifest():
    args = _cli.build_count_args(
        manifest="/tmp/manifest.tsv",
        index="/idx/hg38.bi",
        outfile="/tmp/out.csv.gz",
        lambda_=10.0,
        pseudocount=1.0,
        threads=4,
    )
    assert args == [
        "count",
        "-m", "/tmp/manifest.tsv",
        "-i", "/idx/hg38.bi",
        "-o", "/tmp/out.csv.gz",
        "-l", "10.0",
        "-p", "1.0",
        "-t", "4",
    ]


def test_build_count_args_rejects_mixed_sources():
    with pytest.raises(ValueError, match="not both"):
        _cli.build_count_args(
            manifest="/tmp/m.tsv",
            fragments="/data/f.tsv.gz",
            index="/idx/hg38.bi",
            outfile="/tmp/out.csv.gz",
        )


def test_build_count_args_requires_a_source():
    with pytest.raises(ValueError, match="both required"):
        _cli.build_count_args(index="/idx/hg38.bi", outfile="/tmp/out.csv.gz")


def test_write_manifest(tmp_path: Path):
    path = tmp_path / "manifest.tsv"
    _cli.write_manifest(
        path,
        [("S1", "/data/a.tsv.gz", "/tmp/a.txt"), ("S2", "/data/b.tsv.gz", "/tmp/b.txt")],
    )
    assert path.read_text() == (
        "sample\tfragments\tcells\n"
        "S1\t/data/a.tsv.gz\t/tmp/a.txt\n"
        "S2\t/data/b.tsv.gz\t/tmp/b.txt\n"
    )


def test_write_manifest_rejects_tab_in_sample(tmp_path: Path):
    with pytest.raises(ValueError, match="may not contain a tab"):
        _cli.write_manifest(tmp_path / "m.tsv", [("bad\tname", "f", "c")])


# --- cells-file normalisation (bacha._count) ---------------------------------

def test_normalize_cells_accepts_every_form():
    from bacha import _count

    assert _count._normalize_cells(["a", "b"], what="x") == [("a", "a"), ("b", "b")]
    assert _count._normalize_cells([("a", "s_a")], what="x") == [("a", "s_a")]
    assert _count._normalize_cells([["a", "s_a"]], what="x") == [("a", "s_a")]
    # A mapping preserves insertion order, which is the output column order.
    assert _count._normalize_cells({"b": "s_b", "a": "s_a"}, what="x") == [
        ("b", "s_b"),
        ("a", "s_a"),
    ]


def test_normalize_cells_rejects_bad_entries():
    from bacha import _count

    with pytest.raises(ValueError, match="barcode, or a"):
        _count._normalize_cells([("a", "b", "c")], what="sample 's1'")


def test_write_cells_file_uses_one_column_when_nothing_is_renamed(tmp_path: Path):
    from bacha import _count

    p = tmp_path / "cells.txt"
    _count._write_cells_file(p, [("a", "a"), ("b", "b")])
    assert p.read_text() == "a\nb\n"

    _count._write_cells_file(p, [("a", "s_a"), ("b", "s_b")])
    assert p.read_text() == "a\ts_a\nb\ts_b\n"
