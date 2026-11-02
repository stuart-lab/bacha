"""Locate and invoke the `bacha` Rust CLI."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Sequence


class BachaNotFoundError(RuntimeError):
    """Raised when the `bacha` binary cannot be located on PATH."""


class BachaCLIError(RuntimeError):
    """Raised when the `bacha` binary exits with a non-zero status."""


def locate(bacha_path: str | Path | None = None) -> Path:
    """Return an absolute path to the `bacha` binary.

    Parameters
    ----------
    bacha_path
        If given, validated and returned. Otherwise the binary is located via
        :func:`shutil.which`.
    """
    if bacha_path is not None:
        p = Path(bacha_path).expanduser().resolve()
        if not p.is_file():
            raise BachaNotFoundError(f"bacha binary not found at {p}")
        return p
    found = shutil.which("bacha")
    if not found:
        raise BachaNotFoundError(
            "bacha binary not found on PATH. Install it with "
            "`cargo install bacha` or pass `bacha_path=...`."
        )
    return Path(found).resolve()


def build_count_args(
    *,
    index: str | Path,
    outfile: str | Path,
    fragments: str | Path | None = None,
    cells_file: str | Path | None = None,
    manifest: str | Path | None = None,
    lambda_: float = 10.0,
    pseudocount: float = 1.0,
    threads: int = 6,
    check_bounds: bool = True,
) -> list[str]:
    """Build the argument list for `bacha count`.

    Pass either ``fragments`` + ``cells_file`` for a single fragment file, or
    ``manifest`` to score several files into one matrix. Mirrors
    ``bachaR::Bacha.default`` / ``bachaR:::bacha_multi`` exactly so output is
    bit-identical.
    """
    # Flag order mirrors the R package exactly, so the two front ends invoke
    # the CLI identically.
    if manifest is not None:
        if fragments is not None or cells_file is not None:
            raise ValueError("pass either manifest, or fragments + cells_file, not both")
        source = ["-m", str(manifest), "-i", str(index)]
    else:
        if fragments is None or cells_file is None:
            raise ValueError("fragments and cells_file are both required without a manifest")
        source = ["-f", str(fragments), "-i", str(index), "-c", str(cells_file)]

    args = [
        "count",
        *source,
        "-o", str(outfile),
        "-l", str(lambda_),
        "-p", str(pseudocount),
        "-t", str(threads),
    ]
    if not check_bounds:
        args.append("--no-bounds-check")
    return args


def write_manifest(
    path: str | Path,
    rows: Sequence[tuple[str, str | Path, str | Path]],
) -> None:
    """Write a `bacha count --manifest` TSV.

    ``rows`` is a sequence of ``(sample, fragments_path, cells_path)``. Paths
    are written as given, so absolute paths are the safe choice (relative ones
    are resolved by bacha against the manifest's own directory).
    """
    with open(path, "w") as handle:
        handle.write("sample\tfragments\tcells\n")
        for sample, fragments, cells in rows:
            if "\t" in sample or "\n" in sample:
                raise ValueError(f"sample name may not contain a tab or newline: {sample!r}")
            handle.write(f"{sample}\t{fragments}\t{cells}\n")


def run(
    bacha_path: str | Path,
    args: Sequence[str],
    *,
    verbose: bool = True,
) -> None:
    """Invoke `bacha` with the given arguments; raise on non-zero exit."""
    stdout = None if verbose else subprocess.DEVNULL
    stderr = None if verbose else subprocess.PIPE
    completed = subprocess.run(
        [str(bacha_path), *args],
        stdout=stdout,
        stderr=stderr,
        check=False,
    )
    if completed.returncode != 0:
        captured = completed.stderr.decode() if completed.stderr else ""
        raise BachaCLIError(
            f"bacha exited with code {completed.returncode}"
            + (f":\n{captured}" if captured else "")
        )
