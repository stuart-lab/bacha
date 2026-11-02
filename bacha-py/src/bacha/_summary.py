"""Parse the output of `bacha summary <index>` into a structured dict."""

from __future__ import annotations

import subprocess
from pathlib import Path

from . import _cli


def summary(
    index: str | Path,
    *,
    bacha_path: str | Path | None = None,
) -> dict:
    """Return header info for a `.bi` index file.

    Returned dict has the keys:

    - ``format`` (int): binary format version.
    - ``genome`` (str): genome assembly name.
    - ``motif_set`` (str): motif collection identifier.
    - ``window`` (int): site window size.
    - ``chromosomes`` (list[tuple[str, int]]): (name, size) pairs.
    - ``sites`` (list[tuple[str, int, int]]): (motif_id, column_index, count).
    - ``controls`` (list[tuple[str, int, int]]): same shape as ``sites``.
    """
    binary = _cli.locate(bacha_path)
    completed = subprocess.run(
        [str(binary), "summary", str(Path(index).expanduser())],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise _cli.BachaCLIError(
            f"bacha summary failed (exit {completed.returncode}):\n{completed.stderr}"
        )
    return _parse(completed.stdout)


def _parse(text: str) -> dict:
    result: dict = {
        "format": None,
        "genome": None,
        "motif_set": None,
        "window": None,
        "chromosomes": [],
        "sites": [],
        "controls": [],
    }
    section: str | None = None
    for raw in text.splitlines():
        line = raw.rstrip()
        if not line:
            continue
        if line.startswith("## format:"):
            v = line.split(":", 1)[1].strip()
            result["format"] = int(v.lstrip("binary v")) if "v" in v else int(v)
            continue
        if line.startswith("## genome:"):
            result["genome"] = line.split(":", 1)[1].strip()
            continue
        if line.startswith("## motif set:"):
            result["motif_set"] = line.split(":", 1)[1].strip()
            continue
        if line.startswith("## window:"):
            result["window"] = int(line.split(":", 1)[1].strip())
            continue
        if line.startswith("##[chromosome sizes]"):
            section = "chromosomes"
            continue
        if line.startswith("##[sites]"):
            section = "sites"
            continue
        if line.startswith("##[controls]"):
            section = "controls"
            continue
        if line.startswith("#"):
            continue
        parts = line.split("\t")
        if section == "chromosomes" and len(parts) == 2:
            result["chromosomes"].append((parts[0], int(parts[1])))
        elif section in ("sites", "controls") and len(parts) == 3:
            result[section].append((parts[0], int(parts[1]), int(parts[2])))
    return result
