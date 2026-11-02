#' Bias-adjusted aggregation of chromatin accessibility data
#' 
#' Bacha computes a per-cell aggregated accessibility score for a set of genomic
#' regions (for example, transcription factor motif sites).
#' 
#' Counts use paired insertion counting. A fragment marks two Tn5 insertions,
#' at its start and end coordinates, and each site or control interval
#' containing either insertion is counted once.
#'
#' The raw Bacha score is the log\eqn{_2} fold change in site insertion counts
#' at the query positions \eqn{x_i} relative to a set of control regions
#' \eqn{x_j}. Indexes built with `BachaIndex()` match controls for GC content
#' and sequence entropy by default, with Tn5 bias matching available as an
#' optional index-building feature:
#'
#' \deqn{
#' \text{LFC} = \log_2 \left( \frac{x_i + p}{x_j / r + p} \right)
#' }
#'
#' where \eqn{p} is a pseudocount and \eqn{r} is the ratio of control to
#' site region footprint. The LFC estimate is shrunk by a precision-based
#' weight:
#'
#' \deqn{
#' S = w \cdot \text{LFC}, \quad w = \frac{\text{prec}}{\text{prec} + \lambda},
#' \quad \text{prec} = \frac{(s + p)(c + p r)}{s + c + p + p r}
#' }
#'
#' The pseudocount \eqn{p} corresponds to a Gamma\eqn{(p, \cdot)} prior on each
#' per-bp rate, giving \eqn{p} pseudo-observations on the site footprint and
#' \eqn{p r} pseudo-observations on the (larger) control footprint. `prec` is
#' the resulting posterior precision of the LFC, so the weight is the standard
#' Bayesian posterior shrinkage toward zero under a Normal prior on
#' log-enrichment with prior precision \eqn{\lambda}. Larger values of
#' `lambda` give more aggressive shrinkage. Under this parameterisation, one-
#' sided evidence (e.g. many site insertions and no control insertions) still
#' receives non-zero weight scaling with the control window size, since the
#' control exposure itself is informative about the absence of insertions.
#'
#' The Bacha function requires that the command line Rust tool is available.
#' This can be installed from <https://github.com/stuart-lab/bacha>.
#'
#' When an assay holds several fragment files, each is scored separately and the
#' results are placed side by side. This requires that every cell's fragments
#' live in exactly one file, which is the usual arrangement for multi-sample
#' objects (sharing barcodes between samples is fine). If a cell's fragments are
#' split across files, `Bacha()` errors instead of returning a score: the score
#' is a non-linear function of the insertion counts, so per-file scores for one
#' cell cannot be recombined into a whole-cell score. Merge such fragment files
#' into a single sorted, bgzip-compressed file before running.
#'
#' @param object A Seurat object, ChromatinAssay object, Fragment object, or
#' a path to a fragment file.
#' @param ... Additional arguments passed to other methods
#' 
#' @export Bacha
#' @rdname Bacha
#' @return When running on a `Seurat` object, returns a modified version of
#' the input object with a new assay containing the bacha results. When running
#' on a `ChromatinAssay` object, returns a new Seurat `Assay` object
#' containing the bacha results. In the returned assay, underscores and pipes
#' in motif names are replaced with dashes, and name collisions are made
#' unique. Original index names are retained in the feature metadata column
#' `motif_id`. When running on a fragment file or
#' `Fragment` object, returns a feature x cell matrix.
Bacha <- function(object, ...) {
  UseMethod(generic = "Bacha", object = object)
}
