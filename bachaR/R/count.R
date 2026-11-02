
#' @param cells Vector of cell names to include.
#' @rdname Bacha
#' @export
#' @importFrom data.table fread
Bacha.default <- function(
    object,
    cells,
    index,
    lambda = 10,
    pseudocount = 1,
    threads = 6,
    check.bounds = TRUE,
    bacha.path = NULL,
    outfile = tempfile(fileext = ".csv.gz"),
    cleanup = TRUE,
    verbose = TRUE,
    ...
  ) {

  # check fragment file exists
  if (!file.exists(object)) {
    stop("Requested fragment file does not exist")
  }
  if (file.exists(outfile)) {
    stop("Requested output file already exists")
  }

  # convert index to absolute path
  index <- normalizePath(path = index, mustWork = TRUE)

  # find bacha
  bacha.path <- bacha.path %||% unname(obj = Sys.which(names = "bacha"))
  if (nchar(x = bacha.path) == 0) {
    stop("Bacha not found. Please install Bacha:",
         "https://crates.io/crates/bacha")
  }
  bacha.path <- normalizePath(path = bacha.path, mustWork = TRUE)

  # write cells file. Temporary files are removed on exit rather than only on
  # the success path, so an error in bacha or in reading its output does not
  # leave them behind.
  cellfile <- tempfile()
  if (cleanup) {
    on.exit(expr = unlink(x = cellfile), add = TRUE)
  }
  writeLines(cells, cellfile)

  # run bacha (use system2 so paths with spaces or shell metacharacters
  # are passed as a single argument rather than re-parsed by the shell)
  # When verbose, stream stderr live to the console; otherwise capture it so
  # that a failure can report the underlying error message rather than just an
  # exit code.
  errfile <- if (verbose) "" else tempfile()
  if (!verbose) {
    on.exit(expr = unlink(x = errfile), add = TRUE)
  }
  if (cleanup) {
    on.exit(expr = unlink(x = outfile), add = TRUE)
  }
  exit_code <- system2(
    command = bacha.path,
    args = c(
      "count",
      "-f", object,
      "-i", as.character(x = index),
      "-c", cellfile,
      "-o", outfile,
      "-l", as.character(x = lambda),
      "-p", as.character(x = pseudocount),
      "-t", as.character(x = threads),
      if (!check.bounds) "--no-bounds-check"
    ),
    stdout = if (verbose) "" else FALSE,
    stderr = errfile,
    wait = TRUE
  )
  if (exit_code != 0) {
    err.msg <- if (!verbose && file.exists(errfile)) {
      paste(readLines(con = errfile, warn = FALSE), collapse = "\n")
    } else {
      ""
    }
    stop(
      "bacha command failed with exit code ", exit_code,
      if (nzchar(err.msg)) paste0(":\n", err.msg) else "",
      call. = FALSE
    )
  }
  
  # read output matrix (temporary files are removed by the on.exit handlers
  # registered above)
  bacha.matrix <- as.matrix(
    x = fread(outfile, sep = ",", header = TRUE, nThread = threads),
    rownames = 1
  )

  return(bacha.matrix)
}


#' @rdname Bacha
#' @method Bacha Fragment2
#' @importFrom Signac GetFragmentData
#' @export
Bacha.Fragment2 <- function(
    object,
    index,
    lambda = 10,
    threads = 6,
    bacha.path = NULL,
    outfile = tempfile(fileext = ".csv.gz"),
    cleanup = TRUE,
    verbose = TRUE,
    ...
) {
  # get fragment path
  fpath <- GetFragmentData(object = object, slot = "file.path")

  # get cell vector
  # element is the name of the cell in the fragment file
  # name of each element is the cell name in the assay
  cells.use <- GetFragmentData(object = object, slot = "cells")

  mat <- Bacha(
    object = fpath,
    index = index,
    cells = unname(obj = cells.use),
    lambda = lambda,
    threads = threads,
    bacha.path = bacha.path,
    outfile = outfile,
    cleanup = cleanup,
    verbose = verbose,
    ...
  )

  # rename cells
  cellconvert <- names(x = cells.use)
  names(x = cellconvert) <- cells.use
  colnames(x = mat) <- cellconvert[colnames(x = mat)]
  return(mat)
}

#' @rdname Bacha
#' @method Bacha Fragment
#' @importFrom Signac GetFragmentData
#' @export
Bacha.Fragment <- function(
    object,
    index,
    lambda = 10,
    threads = 6,
    bacha.path = NULL,
    outfile = tempfile(fileext = ".csv.gz"),
    cleanup = TRUE,
    verbose = TRUE,
    ...
  ) {
  # get fragment path
  fpath <- GetFragmentData(object = object, slot = "path")

  # get cell vector
  # element is the name of the cell in the fragment file
  # name of each element is the cell name in the assay
  cells.use <- GetFragmentData(object = object, slot = "cells")

  mat <- Bacha(
    object = fpath,
    index = index,
    cells = unname(obj = cells.use),
    lambda = lambda,
    threads = threads,
    bacha.path = bacha.path,
    outfile = outfile,
    cleanup = cleanup,
    verbose = verbose,
    ...
  )

  # rename cells
  cellconvert <- names(x = cells.use)
  names(x = cellconvert) <- cells.use
  colnames(x = mat) <- cellconvert[colnames(x = mat)]
  return(mat)
}

#' @rdname Bacha
#' @method Bacha ChromatinAssay5
#' @export
Bacha.ChromatinAssay5 <- function(
    object,
    index,
    lambda = 10,
    threads = 6,
    bacha.path = NULL,
    outfile = tempfile(fileext = ".csv.gz"),
    cleanup = TRUE,
    verbose = TRUE,
    ...
) {
  bacha.assay <- Bacha.ChromatinAssay(
    object = object,
    index = index,
    lambda = lambda,
    threads = threads,
    bacha.path = bacha.path,
    outfile = outfile,
    cleanup = cleanup,
    verbose = verbose,
    ...
  )
  return(bacha.assay)
}

#' @rdname Bacha
#' @method Bacha ChromatinAssay
#' @importFrom Signac Fragments
#' @importFrom SeuratObject CreateAssay5Object LayerData<-
#' @importFrom Matrix sparseMatrix
#' @export
Bacha.ChromatinAssay <- function(
    object,
    index,
    lambda = 10,
    threads = 6,
    bacha.path = NULL,
    outfile = tempfile(fileext = ".csv.gz"),
    cleanup = TRUE,
    verbose = TRUE,
    ...
  ) {

  frags <- Fragments(object = object)

  if (length(x = frags) == 0) {
    stop("No fragment files present in the requested assay")
  }

  # All fragment files scored in a single bacha call
  bm.all <- bacha_multi(
    frags = frags,
    index = index,
    lambda = lambda,
    threads = threads,
    bacha.path = bacha.path,
    outfile = outfile,
    cleanup = cleanup,
    verbose = verbose,
    ...
  )

  # Cells in the assay that no fragment file scored are reported as zero. Build
  # the full matrix in one allocation and fill the scored columns into it
  assay.cells <- colnames(x = object)
  if (!identical(x = colnames(x = bm.all), y = assay.cells)) {
    scored <- intersect(x = assay.cells, y = colnames(x = bm.all))
    n.missing <- length(x = assay.cells) - length(x = scored)
    if (n.missing > 0 && verbose) {
      message("Filling ", n.missing,
              " cell(s) absent from the fragment file(s) with zero")
    }
    full <- matrix(
      data = 0,
      nrow = nrow(x = bm.all),
      ncol = length(x = assay.cells),
      dimnames = list(rownames(x = bm.all), assay.cells)
    )
    full[, scored] <- bm.all[, scored, drop = FALSE]
    bm.all <- full
  }

  # Normalize names before both the placeholder and the dense layer are
  # created: Seurat replaces underscores and pipes with dashes. Resolve
  # collisions (e.g. TF_A and TF-A) and retain the original index labels.
  motif.ids <- rownames(x = bm.all)
  rownames(x = bm.all) <- make.unique(names = gsub(
    pattern = "[_|]", replacement = "-", x = motif.ids
  ))

  # Create the assay, keeping the scores as a dense matrix.
  #
  # Bacha scores are a log fold change per cell per feature, so essentially
  # every entry is non-zero. `CreateAssay5Object()` coerces whatever it is
  # given to dgCMatrix, which for a matrix with no zeros costs 12 bytes per
  # entry (value + row index, plus the column pointers) against 8 for dense
  # 1.5x the memory, and slower row-wise access.
  #
  # Assay5 layers themselves accept any matrix-like object (this is how dense
  # scale.data layers are stored), so seed the assay with an empty sparse
  # matrix of the right shape and attach the scores with `LayerData<-`, which
  # stores them as given. Passing `bm.all` to CreateAssay5Object() and then
  # overwriting would work but builds the full sparse copy first and
  # doubles peak memory.
  placeholder <- sparseMatrix(
    i = integer(0), j = integer(0),
    dims = dim(x = bm.all), dimnames = dimnames(x = bm.all)
  )
  bacha.assay <- CreateAssay5Object(data = placeholder)
  LayerData(object = bacha.assay, layer = "data") <- bm.all
  bacha.assay[["motif_id"]] <- setNames(
    object = motif.ids, nm = rownames(x = bm.all)
  )
  return(bacha.assay)
}

# Path to the fragment file backing a Fragment (Signac v1) or Fragment2
# (Signac v2) object. The two classes name the slot differently.
#' @importFrom Signac GetFragmentData
fragment_file_path <- function(object) {
  slot <- if (inherits(x = object, what = "Fragment2")) "file.path" else "path"
  GetFragmentData(object = object, slot = slot)
}

# Score several fragment files in one bacha call.
#
# `frags` is a list of Fragment/Fragment2 objects. Each one carries a `cells`
# vector whose values are the barcodes as they appear in that fragment file and
# whose names are the corresponding cell names in the assay. That mapping is
# written straight into bacha's two-column cells file, so bacha labels its
# output columns with the assay's cell names and no renaming is needed here.
#
# Returns the feature x cell score matrix, with columns in manifest order.
#' @importFrom data.table fread
bacha_multi <- function(
    frags,
    index,
    lambda = 10,
    pseudocount = 1,
    threads = 6,
    check.bounds = TRUE,
    bacha.path = NULL,
    outfile = tempfile(fileext = ".csv.gz"),
    cleanup = TRUE,
    verbose = TRUE,
    ...
  ) {

  if (file.exists(outfile)) {
    stop("Requested output file already exists")
  }
  index <- normalizePath(path = index, mustWork = TRUE)

  bacha.path <- bacha.path %||% unname(obj = Sys.which(names = "bacha"))
  if (nchar(x = bacha.path) == 0) {
    stop("Bacha not found. Please install Bacha:",
         "https://crates.io/crates/bacha")
  }
  bacha.path <- normalizePath(path = bacha.path, mustWork = TRUE)

  # Temporary files are removed on exit rather than only on the success path,
  # so an error in bacha or in reading its output does not leave them behind.
  manifest <- tempfile(fileext = ".tsv")
  cellfiles <- replicate(n = length(x = frags), expr = tempfile())
  if (cleanup) {
    on.exit(expr = unlink(x = c(manifest, cellfiles)), add = TRUE)
  }

  rows <- vector(mode = "list", length = length(x = frags))
  for (i in seq_along(along.with = frags)) {
    fpath <- fragment_file_path(object = frags[[i]])
    if (!file.exists(fpath)) {
      stop("Fragment file does not exist: ", fpath, call. = FALSE)
    }
    cells.use <- GetFragmentData(object = frags[[i]], slot = "cells")
    if (length(x = cells.use) == 0) {
      stop("Fragment file ", i, " has no cells associated with it",
           call. = FALSE)
    }
    # fall back to the barcodes themselves so the output columns are still labelled.
    out.names <- names(x = cells.use) %||% unname(obj = cells.use)
    if (anyNA(x = out.names) || any(!nzchar(x = out.names))) {
      stop("Fragment file ", i, " has cells with a missing or empty name",
           call. = FALSE)
    }
    # Column 1: barcode in the fragment file. Column 2: name to use in the output
    utils::write.table(
      x = data.frame(
        barcode = unname(obj = cells.use),
        cell = out.names,
        stringsAsFactors = FALSE
      ),
      file = cellfiles[[i]],
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
    rows[[i]] <- data.frame(
      sample = paste0("file", i, "_",
                      tools::file_path_sans_ext(x = basename(path = fpath))),
      fragments = normalizePath(path = fpath, mustWork = TRUE),
      cells = cellfiles[[i]],
      stringsAsFactors = FALSE
    )
  }
  utils::write.table(
    x = do.call(what = rbind, args = rows),
    file = manifest,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )

  errfile <- if (verbose) "" else tempfile()
  if (!verbose) {
    on.exit(expr = unlink(x = errfile), add = TRUE)
  }
  if (cleanup) {
    on.exit(expr = unlink(x = outfile), add = TRUE)
  }
  exit_code <- system2(
    command = bacha.path,
    args = c(
      "count",
      "-m", manifest,
      "-i", as.character(x = index),
      "-o", outfile,
      "-l", as.character(x = lambda),
      "-p", as.character(x = pseudocount),
      "-t", as.character(x = threads),
      if (!check.bounds) "--no-bounds-check"
    ),
    stdout = if (verbose) "" else FALSE,
    stderr = errfile,
    wait = TRUE
  )
  if (exit_code != 0) {
    err.msg <- if (!verbose && file.exists(errfile)) {
      paste(readLines(con = errfile, warn = FALSE), collapse = "\n")
    } else {
      ""
    }
    stop(
      "bacha command failed with exit code ", exit_code,
      if (nzchar(err.msg)) paste0(":\n", err.msg) else "",
      call. = FALSE
    )
  }

  as.matrix(x = fread(outfile, sep = ",", header = TRUE, nThread = threads), rownames = 1)
}

#' @rdname Bacha
#' @method Bacha Assay
#' @export
Bacha.Assay <- function(object, ...) {
  stop("Bacha requires a ChromatinAssay with associated fragment files, ",
       "not a standard Assay")
}

#' @rdname Bacha
#' @method Bacha Assay5
#' @export
Bacha.Assay5 <- function(object, ...) {
  stop("Bacha requires a ChromatinAssay with associated fragment files, ",
       "not a standard Assay5")
}

#' @param index Path to the Bacha index file.
#' @param lambda Numeric. Prior precision on log-enrichment, used in the
#' precision-based shrinkage weight `prec / (prec + lambda)`, where
#' `prec = (s + p)(c + p*r) / (s + c + p + p*r)` is the exposure-consistent
#' posterior precision of the LFC under a Gamma(p, .) prior on each per-bp
#' rate. Larger values give stronger shrinkage toward zero. Default: 10.
#' @param pseudocount Numeric. Small constant added to counts before log fold
#' change calculation to avoid \code{log(0)}. Shrinkage is handled by the
#' precision weight, so this only needs to be a numerical guard. Default: 1.
#' @param threads Integer. Total number of threads to use, passed to bacha as
#' \code{-t}. Against a bgzip-compressed fragment file, 2 of them inflate the
#' file and the rest count it, so the two pools share \code{threads} rather
#' than adding to it. Default: 6.
#' @param check.bounds Logical. If TRUE (the default), bacha aborts as soon as
#' it sees a fragment coordinate that lies beyond the end of its chromosome in
#' the index, which indicates the fragment file and index were built on
#' different genome assemblies. Set to FALSE to skip this check (passes
#' \code{--no-bounds-check} to bacha).
#' @param assay Name of assay to use. If NULL, the default assay will be used.
#' @param bacha.path Path to bacha executable. If NULL, try to find bacha automatically.
#' @param outfile Path to store temporary output file.
#' @param cleanup Remove temporary file created by bacha.
#' @param new.assay.name Name to store new assay in the Seurat object containing
#' bacha scores per cell.
#' @param verbose Display messages
#' @rdname Bacha
#' @importFrom SeuratObject DefaultAssay
#' @method Bacha Seurat
#' @export
Bacha.Seurat <- function(
    object,
    index,
    lambda = 10,
    threads = 6,
    assay = NULL,
    bacha.path = NULL,
    outfile = tempfile(fileext = ".csv.gz"),
    cleanup = TRUE,
    new.assay.name = "bacha",
    verbose = TRUE,
    ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  bacha.assay <- Bacha(
    object = object[[assay]],
    index = index,
    lambda = lambda,
    threads = threads,
    bacha.path = bacha.path,
    outfile = outfile,
    cleanup = cleanup,
    verbose = verbose,
    ...
  )
  object[[new.assay.name]] <- bacha.assay
  return(object)
}
