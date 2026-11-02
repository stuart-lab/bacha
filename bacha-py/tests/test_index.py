"""Unit tests for bacha.index._summary._parse."""

from __future__ import annotations

from bacha._summary import _parse


SAMPLE = """## format: binary v3
## genome: hg38
## motif set: H14CORE
## window: 200
##[chromosome sizes]
chr1\t248956422
chr2\t242193529
##[sites]
# id\tindex\tcount
MOTIF_A\t1\t12345
MOTIF_B\t2\t6789
##[controls]
# id\tindex\tcount
MOTIF_A\t3\t99999
MOTIF_B\t4\t88888
"""


def test_parse_extracts_header():
    out = _parse(SAMPLE)
    assert out["format"] == 3
    assert out["genome"] == "hg38"
    assert out["motif_set"] == "H14CORE"
    assert out["window"] == 200


def test_parse_extracts_chromosomes():
    out = _parse(SAMPLE)
    assert out["chromosomes"] == [("chr1", 248956422), ("chr2", 242193529)]


def test_parse_extracts_sites_and_controls():
    out = _parse(SAMPLE)
    assert out["sites"] == [("MOTIF_A", 1, 12345), ("MOTIF_B", 2, 6789)]
    assert out["controls"] == [("MOTIF_A", 3, 99999), ("MOTIF_B", 4, 88888)]
