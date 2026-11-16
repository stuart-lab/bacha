#' Summarise a Bacha index file
#'
#' Prints a summary of a Bacha index file: format version, genome and motif
#' set, window size, how the index was built (date and [BachaIndex()]
#' parameters), number of chromosomes covered, number of TFs, and a summary
#' of the per-TF number of motif sites. Optionally returns a data
#' frame of TF site counts.
#'
#' Requires the `bacha` command-line tool to be installed and available on
#' the system PATH (or passed via `bacha.path`), since the index file is
#' parsed via `bacha view`.
#'
#' @param index Path to a Bacha index file (`.bi`).
#' @param bacha.path Path to the `bacha` executable. If `NULL`, looks it up
#'   on the system PATH.
#' @param return.table Logical. If `TRUE`, invisibly return a data frame
#'   with columns `tf` (feature ID), `name` (display name) and `n_sites` in
#'   addition to printing the summary.
#'   Default: `FALSE`.
#'
#' @return Invisibly, `NULL` (or a data frame if `return.table = TRUE`).
#' @importFrom stats quantile
#' @importFrom utils read.table
#' @export
IndexSummary <- function(
    index,
    bacha.path = NULL,
    return.table = FALSE
) {
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

  lines <- system2(bacha.path, c("summary", index), stdout = TRUE, stderr = FALSE)

  # Header metadata: lines starting with "## "
  header_lines <- grep("^##[^\\[]", lines, value = TRUE)
  header_meta <- sub("^##\\s*", "", header_lines)

  # Build provenance: "key\tvalue" lines between "##[build]" and the next
  # "#" line. Absent when the index records none.
  build_block <- section_lines(lines = lines, header = "##[build]")
  build <- sub(pattern = "\t", replacement = ": ", x = build_block, fixed = TRUE)

  # Feature table: "id<TAB>name" lines after the "# id<TAB>name" header
  features <- section_lines(lines = lines, header = "##[features]")
  feature_name <- if (length(x = features)) {
    parts <- strsplit(x = features, split = "\t", fixed = TRUE)
    setNames(
      object = vapply(X = parts, FUN = function(p) p[2] %||% p[1], FUN.VALUE = ""),
      nm = vapply(X = parts, FUN = `[`, FUN.VALUE = "", 1L)
    )
  } else {
    character(0)
  }

  # Chromosome sizes: between "##[chromosome sizes]" and the next "#" line
  chrom_hdr <- grep("^##\\[chromosome sizes\\]", lines)
  n_chroms <- 0L
  if (length(chrom_hdr) == 1) {
    tail_lines <- lines[(chrom_hdr + 1):length(lines)]
    first_hash <- which(startsWith(tail_lines, "#"))
    chrom_block <- if (length(first_hash)) {
      tail_lines[seq_len(first_hash[1] - 1)]
    } else {
      tail_lines
    }
    n_chroms <- length(chrom_block)
  }

  # TF site table: block after the first "# id\tindex\tcount"
  # (a second identical header starts the controls block)
  tf_hdr <- grep("# id\tindex\tcount", lines, fixed = TRUE)
  if (length(tf_hdr) == 0) {
    stop("Could not locate TF table in index file")
  }
  after <- lines[(tf_hdr[1] + 1):length(lines)]
  first_hash <- which(startsWith(after, "#"))
  tf_block <- if (length(first_hash)) {
    after[seq_len(first_hash[1] - 1)]
  } else {
    after
  }
  tf_tab <- read.table(
    text = tf_block, sep = "\t", header = FALSE,
    col.names = c("tf", "idx", "n_sites"), stringsAsFactors = FALSE
  )
  tf_tab$name <- unname(feature_name[tf_tab$tf])
  tf_tab$name[is.na(x = tf_tab$name)] <- tf_tab$tf[is.na(x = tf_tab$name)]

  # Print summary
  cat("Bacha index:", index, "\n")
  for (m in header_meta) cat("  ", m, "\n", sep = "")
  if (length(x = build) > 0L) {
    cat("  build:\n")
    for (b in build) cat("    ", b, "\n", sep = "")
  }
  cat("  chromosomes: ", n_chroms, "\n", sep = "")
  cat("  TFs:         ", nrow(tf_tab), "\n", sep = "")
  cat("  motif sites per TF:\n")
  qs <- quantile(tf_tab$n_sites, probs = c(0, 0.25, 0.5, 0.75, 1))
  cat(sprintf(
    "    min = %d, Q1 = %.0f, median = %.0f, Q3 = %.0f, max = %d, mean = %.0f\n",
    qs[1], qs[2], qs[3], qs[4], qs[5], mean(tf_tab$n_sites)
  ))

  top <- tf_tab[order(-tf_tab$n_sites), ][seq_len(min(5, nrow(tf_tab))), ]
  bot <- tf_tab[order(tf_tab$n_sites), ][seq_len(min(5, nrow(tf_tab))), ]
  cat("  top 5 (most sites): ",
      paste(sprintf("%s (%d)", top$name, top$n_sites), collapse = ", "),
      "\n", sep = "")
  cat("  bottom 5 (fewest): ",
      paste(sprintf("%s (%d)", bot$name, bot$n_sites), collapse = ", "),
      "\n", sep = "")

  if (return.table) {
    return(invisible(tf_tab[, c("tf", "name", "n_sites")]))
  }
  invisible(NULL)
}

# Data lines of a "##[name]" section of `bacha view` output, up to the next
# line starting with "#". Empty if the section is absent.
section_lines <- function(lines, header) {
  hdr <- which(lines == header)
  if (length(x = hdr) != 1L || hdr == length(x = lines)) {
    return(character(0))
  }
  rest <- lines[(hdr + 1L):length(x = lines)]
  # skip a "# col1<TAB>col2" column header directly under the section header
  if (startsWith(x = rest[1], prefix = "# ")) {
    rest <- rest[-1L]
  }
  stop_at <- which(startsWith(x = rest, prefix = "#"))
  if (length(x = stop_at)) rest[seq_len(stop_at[1] - 1L)] else rest
}
