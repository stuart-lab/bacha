# Bacha

Python interface to **Bacha**: bias-adjusted aggregatino of chromatin accessibility data.

This package mirrors the R package `bachaR`: it wraps the `bacha` Rust CLI and
returns per-cell-per-motif scores as an `AnnData` object.

## Install

```bash
# 1) Install the bacha CLI (Rust)
cargo install bacha

# 2) Install the Python package
pip install bacha
```

The Python package shells out to the `bacha` binary, so the CLI must be on
your `PATH` (or you can pass `bacha_path=...`).

## Usage

```python
import anndata as ad
import bacha

adata = ad.read_h5ad("pbmc.h5ad")

# Returns a new AnnData of shape (n_cells, n_motifs)
result = bacha.count(
    adata,
    fragments="pbmc/fragments.tsv.gz",
    index="hg38.h14.v1.bi",
    inplace=False,
)
```

Or write the scores back into the input object:

```python
bacha.count(
    adata,
    fragments="pbmc/fragments.tsv.gz",
    index="hg38.h14.v1.bi",
    inplace=True,  # scores in adata.obsm["bacha"], motif IDs in adata.uns["bacha_motifs"]
)
```

## Relationship to `bachaR`

`bacha.count()` is feature-equivalent to `bachaR::Bacha()` for scoring (same
CLI arguments, same output). Index *construction* (`BachaIndex()`) is
currently R-only.
