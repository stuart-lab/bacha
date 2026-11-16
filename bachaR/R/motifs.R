#' Extract motif matrices from a Bacha index
#'
#' Returns the motif matrices stored in a Bacha index built from PWMs with
#' [BachaIndex()], so that the motifs used to build the index travel with it.
#' Indexes built from `positions` store no matrices.
#'
#' Requires the `bacha` command-line tool, since the index file is read via
#' `bacha motifs`.
#'
#' @param index Path to a Bacha index file (`.bi`).
#' @param bacha.path Path to the `bacha` executable. If `NULL`, looks it up
#'   on the system PATH.
#'
#' @return A list of numeric matrices named by feature ID, each with rows
#'   `A`, `C`, `G`, `T` and one column per motif position. Each matrix has a
#'   `name` attribute (the feature's display name) and a `type` attribute
#'   giving the values stored: `"counts"`, `"probabilities"`, `"logodds"` or
#'   `"other"`. The list can be passed to [Signac::CreateMotifObject()].
#' @export
IndexMotifs <- function(index, bacha.path = NULL) {
  if (!file.exists(index)) {
    stop("Index file does not exist: ", index)
  }
  index <- normalizePath(index, mustWork = TRUE)
  bacha.path <- bacha.path %||% unname(Sys.which("bacha"))
  if (nchar(bacha.path) == 0) {
    stop("bacha executable not found. Install from ",
         "https://github.com/stuart-lab/bacha")
  }
  bacha.path <- normalizePath(bacha.path, mustWork = TRUE)

  err <- tempfile()
  on.exit(expr = unlink(x = err), add = TRUE)
  lines <- suppressWarnings(
    system2(bacha.path, c("motifs", shQuote(index)), stdout = TRUE, stderr = err)
  )
  status <- attr(x = lines, which = "status")
  if (!is.null(x = status) && status != 0) {
    stop("bacha motifs failed: ", paste(readLines(con = err), collapse = "\n"),
         call. = FALSE)
  }
  parse_jaspar_motifs(lines = lines)
}

# Parse `bacha motifs` output: records of a ">ID<TAB>NAME<TAB>TYPE" header
# followed by one "B  [ v1 v2 ... ]" line per base.
parse_jaspar_motifs <- function(lines) {
  starts <- which(startsWith(x = lines, prefix = ">"))
  out <- lapply(X = starts, FUN = function(i) {
    hdr <- strsplit(x = substring(text = lines[i], first = 2L), split = "\t",
                    fixed = TRUE)[[1]]
    rows <- lines[i + seq_len(length.out = 4L)]
    vals <- lapply(X = rows, FUN = function(r) {
      body <- sub(pattern = "^[ACGT]\\s*\\[\\s*(.*?)\\s*\\]\\s*$", replacement = "\\1", x = r)
      as.numeric(x = strsplit(x = body, split = "\\s+")[[1]])
    })
    m <- do.call(what = rbind, args = vals)
    dimnames(x = m) <- list(c("A", "C", "G", "T"), NULL)
    attr(x = m, which = "name") <- hdr[2]
    attr(x = m, which = "type") <- hdr[3]
    m
  })
  names(x = out) <- vapply(
    X = starts,
    FUN = function(i) strsplit(x = substring(text = lines[i], first = 2L), split = "\t",
                               fixed = TRUE)[[1]][1],
    FUN.VALUE = ""
  )
  out
}
