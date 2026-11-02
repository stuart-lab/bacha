"""Run Bacha on a single-cell ATAC dataset and cluster cells on the scores.

This mirrors the scverse workflow: score per-cell motif accessibility with
``bacha.count`` (which shells out to the ``bacha`` Rust CLI), store the result
in ``adata.obsm``, then run the standard neighbors -> UMAP -> Leiden pipeline
on the Bacha embedding.

Prerequisites
-------------
1. The ``bacha`` CLI on your PATH (``cargo install bacha``), or pass
   ``bacha_path=...`` to ``bacha.count``.
2. ``pip install bacha scanpy`` (scanpy is only needed for the downstream
   clustering at the bottom; the scoring itself needs only anndata + numpy).
3. Inputs:
   - an AnnData with one row per cell barcode (``.obs_names`` are the barcodes),
   - a tabix-indexed fragments file (``fragments.tsv.gz`` + ``.tsv.gz.tbi``),
   - a Bacha index (``.bi``) for the matching genome / motif set, from
     https://github.com/stuart-lab/bacha-index

"""

from __future__ import annotations
from pathlib import Path
import anndata as ad
import bacha

# --- paths -----------------------------------------------------------------
FRAGMENTS = "fragments.tsv.gz" 
INDEX = "hg38.h14.v1.bi"
H5AD = "pbmc.h5ad"


def main() -> None:
    # 0) Peek at the index header (genome, motif set, window, #motifs) -------
    info = bacha.summary(INDEX)
    print(f"index: genome={info['genome']} motif_set={info['motif_set']} "
          f"window={info['window']}bp n_motifs={len(info['sites'])}")

    # 1) Load cells -----------------------------------------------------------
    adata = ad.read_h5ad(H5AD)
    print(f"loaded {adata.n_obs} cells")

    # 2) Score in place: per-cell motif scores land in adata.obsm["bacha"] ---
    #    Motif IDs (the columns of that matrix) go to adata.uns["bacha_motifs"].
    #    obs_names are used as the cell list, so the rows stay aligned.
    bacha.count(
        adata,
        fragments=FRAGMENTS,
        index=INDEX,
        threads=8,
        inplace=True,
    )

    scores = adata.obsm["bacha"]              # (n_cells, n_motifs) ndarray
    motifs = adata.uns["bacha_motifs"]        # length n_motifs
    print(f"bacha scores: {scores.shape[0]} cells x {scores.shape[1]} motifs")

    # 3) (Alternative) get the scores as a standalone (cells x motifs) AnnData
    #    instead of writing them back. Its .X *is* the score matrix, so the
    #    standard scanpy pipeline runs directly on it.
    tfa = bacha.count(
        adata,
        fragments=FRAGMENTS,
        index=INDEX,
        inplace=False,
    )

    # 4) Downstream: cluster cells on the Bacha scores -----------------------
    #    Operate on `tfa`, whose .X is the (cells x motifs) score matrix.
    #    Scaling each motif to unit variance before PCA keeps a few
    #    high-variance motifs from dominating the embedding.
    try:
        import scanpy as sc
    except ImportError:
        print("scanpy not installed; skipping clustering")
        return

    sc.pp.scale(tfa)
    sc.pp.pca(tfa, n_comps=30)
    sc.pp.neighbors(tfa, n_pcs=20)
    sc.tl.umap(tfa)
    sc.tl.leiden(tfa, resolution=0.5, key_added="bacha_clusters")
    print(tfa.obs["bacha_clusters"].value_counts())


if __name__ == "__main__":
    main()
