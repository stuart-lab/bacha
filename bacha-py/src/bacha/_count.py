"""bacha.count — run Bacha scoring on an AnnData."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from anndata import AnnData

from . import _cli, _io

# A cell list, as accepted for one fragment file. Either bare barcodes, or
# (barcode, output cell name) pairs, or a {barcode: output cell name} mapping.
CellSpec = "Sequence[str] | Sequence[tuple[str, str]] | Mapping[str, str]"


def _normalize_cells(spec, *, what: str) -> list[tuple[str, str]]:
    """Return (barcode, output cell name) pairs from any accepted cells form.

    Bare barcodes get the barcode as their output name, which is what a
    single-sample run wants. The paired forms exist because barcodes repeat
    across 10x samples that share a whitelist: those are different cells, and
    they need distinct output names to be scored side by side.
    """
    if isinstance(spec, Mapping):
        pairs = [(str(b), str(n)) for b, n in spec.items()]
    else:
        pairs = []
        for item in spec:
            if isinstance(item, str):
                pairs.append((item, item))
            elif isinstance(item, (tuple, list)) and len(item) == 2:
                pairs.append((str(item[0]), str(item[1])))
            else:
                raise ValueError(
                    f"{what}: each cell must be a barcode, or a (barcode, name) pair; "
                    f"got {item!r}"
                )
    return pairs


def _write_cells_file(path: str | Path, pairs: Sequence[tuple[str, str]]) -> None:
    """Write a bacha cells file for `pairs`.

    Emits one column when every output name equals its barcode, so a run that
    does not rename anything produces exactly the file earlier versions wrote.
    """
    renamed = any(barcode != name for barcode, name in pairs)
    if renamed:
        body = "".join(f"{barcode}\t{name}\n" for barcode, name in pairs)
    else:
        body = "".join(f"{barcode}\n" for barcode, _ in pairs)
    Path(path).write_text(body)


def count(
    adata: AnnData,
    *,
    fragments: str | Path | Mapping[str, str | Path],
    index: str | Path,
    cells: Sequence[str]
    | Sequence[tuple[str, str]]
    | Mapping[str, str]
    | Mapping[str, Sequence[str] | Sequence[tuple[str, str]] | Mapping[str, str]]
    | None = None,
    lambda_: float = 10.0,
    pseudocount: float = 1.0,
    threads: int = 6,
    check_bounds: bool = True,
    bacha_path: str | Path | None = None,
    outfile: str | Path | None = None,
    cleanup: bool = True,
    obsm_key: str = "bacha",
    inplace: bool = True,
    verbose: bool = True,
) -> AnnData | None:
    """Run Bacha scoring on an AnnData and attach the result.

    Shells out to the ``bacha`` Rust CLI to score per-cell motif accessibility
    against matched control regions, then either writes the resulting
    (cell × motif) matrix into ``adata.obsm[obsm_key]`` or returns it as a new
    ``AnnData``.

    Parameters
    ----------
    adata
        AnnData whose ``.obs_names`` provide the default cell list and (when
        ``inplace=True``) the target for the result.
    fragments
        Path to a tabix-indexed fragments TSV (``.tsv.gz``), or a mapping of
        ``{sample: path}`` to score several fragment files into one matrix. In
        the multi-file case ``cells`` must be a matching ``{sample: barcodes}``
        mapping, and every cell must belong to exactly one file: a Bacha score
        is a non-linear function of a cell's insertion counts, so counts for
        one cell split across two files cannot be combined. bacha loads the
        index once for the whole run rather than once per file.
    index
        Path to a Bacha index (``.bi``) for the matching genome/motif set.
    cells
        Cells to score. Defaults to ``list(adata.obs_names)``. Each cell list
        may be bare barcodes, ``(barcode, output_name)`` pairs, or a
        ``{barcode: output_name}`` mapping; the paired forms are what let two
        samples that share a 10x barcode whitelist be scored side by side,
        since those are different cells and need distinct output names. Required
        as a ``{sample: cells}`` mapping when ``fragments`` names several files.
    lambda_
        Shrinkage strength passed to ``bacha count -l``. Default 10.0.
    pseudocount
        Log-zero guard passed to ``bacha count -p``. Default 1.0.
    threads
        Total threads for the CLI (``-t``). Against a bgzip-compressed fragment
        file 2 of them inflate it and the rest count it, so the two pools share
        this budget rather than adding to it. Default 6.
    check_bounds
        When ``True`` (default), abort if a fragment coordinate lies beyond the
        end of its chromosome in the index, which signals that the fragment
        file and index were built on different genome assemblies. Set to
        ``False`` to skip the check (passes ``--no-bounds-check``).
    bacha_path
        Path to the ``bacha`` binary. Defaults to ``shutil.which("bacha")``.
    outfile
        Where to write the CSV scores. Defaults to a temporary ``.csv.gz``.
        Must not already exist if specified.
    cleanup
        Remove the temporary cells file (and the outfile, if we created it)
        on success. Default ``True``.
    obsm_key
        Key under ``adata.obsm`` to store the result when ``inplace=True``.
        Motif IDs are written to ``adata.uns[f"{obsm_key}_motifs"]``.
    inplace
        If ``True``, modify ``adata`` and return ``None``. If ``False``,
        return a new ``AnnData`` of shape (n_cells, n_motifs).
    verbose
        Pass through CLI stdout/stderr when ``True``.

    Returns
    -------
    AnnData or None
        ``None`` when ``inplace=True``; a new AnnData otherwise.
    """
    multi = isinstance(fragments, Mapping)
    if multi:
        if not isinstance(cells, Mapping):
            raise ValueError(
                "`cells` must be a {sample: barcodes} mapping when `fragments` names "
                "several fragment files."
            )
        if set(cells) != set(fragments):
            raise ValueError(
                f"`fragments` and `cells` must name the same samples; got "
                f"{sorted(fragments)} and {sorted(cells)}."
            )
        samples = list(fragments)
        frag_paths = {}
        for sample in samples:
            p = Path(fragments[sample]).expanduser()
            if not p.is_file():
                raise FileNotFoundError(f"sample {sample!r}: fragments file not found: {p}")
            frag_paths[sample] = p.resolve()
    else:
        # A mapping here is unambiguous: with one fragment file there is no
        # sample dimension, so it can only be {barcode: output name}.
        single = Path(fragments).expanduser()
        if not single.is_file():
            raise FileNotFoundError(f"fragments file not found: {single}")
        samples = ["sample1"]
        frag_paths = {"sample1": single}

    index_path = Path(index).expanduser().resolve(strict=True)

    binary = _cli.locate(bacha_path)

    # Per-sample (barcode, output name) pairs, and the flat output-cell order
    # they imply. The flat order is what bacha's output columns follow, and what
    # the result is aligned against.
    if multi:
        per_sample = {
            s: _normalize_cells(cells[s], what=f"sample {s!r}") for s in samples
        }
        for sample, pairs in per_sample.items():
            if not pairs:
                raise ValueError(f"sample {sample!r}: no cells to score.")
    else:
        spec = cells if cells is not None else list(map(str, adata.obs_names))
        per_sample = {"sample1": _normalize_cells(spec, what="cells")}

    target_cells = [name for s in samples for _, name in per_sample[s]]
    if not target_cells:
        raise ValueError("No cells to score (adata has 0 obs and `cells` not provided).")
    duplicated = sorted({c for c in target_cells if target_cells.count(c) > 1})
    if duplicated:
        raise ValueError(
            f"{len(duplicated)} output cell name(s) are used more than once "
            f"(first few: {duplicated[:5]}). Barcodes shared between samples are fine, "
            "but each cell needs a distinct output name: pass `cells` as "
            "{sample: {barcode: name}} or as (barcode, name) pairs."
        )

    # Validated before running bacha, so an incompatible request fails in
    # milliseconds rather than after the whole count.
    if inplace and target_cells != list(map(str, adata.obs_names)):
        raise ValueError(
            "`inplace=True` requires the scored cells to match `adata.obs_names`. "
            "Pass `inplace=False` (or omit `cells`) instead."
        )

    # Temp files are registered before anything that can raise, so nothing is
    # left behind if the outfile check below rejects the request.
    tempfiles: list[str] = []
    outfile_path: str | None = None
    created_outfile = False
    try:
        cellfiles: dict[str, str] = {}
        for sample in samples:
            fd, path = tempfile.mkstemp(prefix="bacha_cells_", suffix=".txt")
            os.close(fd)
            tempfiles.append(path)
            _write_cells_file(path, per_sample[sample])
            cellfiles[sample] = path

        manifest_path: str | None = None
        if multi:
            fd, manifest_path = tempfile.mkstemp(prefix="bacha_manifest_", suffix=".tsv")
            os.close(fd)
            tempfiles.append(manifest_path)
            _cli.write_manifest(
                manifest_path,
                [(s, frag_paths[s], cellfiles[s]) for s in samples],
            )

        if outfile is None:
            out_fd, outfile_path = tempfile.mkstemp(prefix="bacha_scores_", suffix=".csv.gz")
            os.close(out_fd)
            # mkstemp creates the file; bacha count refuses an existing path, so remove it.
            os.unlink(outfile_path)
            created_outfile = True
        else:
            outfile_path = str(Path(outfile).expanduser())
            if os.path.exists(outfile_path):
                raise FileExistsError(f"outfile already exists: {outfile_path}")

        if multi:
            args = _cli.build_count_args(
                manifest=manifest_path,
                index=index_path,
                outfile=outfile_path,
                lambda_=lambda_,
                pseudocount=pseudocount,
                threads=threads,
                check_bounds=check_bounds,
            )
        else:
            args = _cli.build_count_args(
                fragments=frag_paths["sample1"],
                index=index_path,
                cells_file=cellfiles["sample1"],
                outfile=outfile_path,
                lambda_=lambda_,
                pseudocount=pseudocount,
                threads=threads,
                check_bounds=check_bounds,
            )
        _cli.run(binary, args, verbose=verbose)

        scores, motif_ids, cell_ids = _io.read_scores(outfile_path)
        scores = _io.align_cells(scores, cell_ids, target_cells)
    finally:
        if cleanup:
            for p in (*tempfiles, outfile_path if created_outfile else None):
                if p and os.path.exists(p):
                    try:
                        os.unlink(p)
                    except OSError:
                        pass

    cell_x_motif = scores.T  # (n_cells, n_motifs)

    if inplace:
        # Compatibility with `adata.obs_names` was checked before bacha ran.
        adata.obsm[obsm_key] = cell_x_motif
        adata.uns[f"{obsm_key}_motifs"] = np.asarray(motif_ids, dtype=object)
        return None

    return AnnData(
        X=cell_x_motif,
        obs=pd.DataFrame(index=pd.Index(target_cells, name=adata.obs_names.name)),
        var=pd.DataFrame(index=pd.Index(motif_ids, name="motif")),
    )
