#' Build Bacha Index
#'
#' Create a bacha index file containing the position of motif match sites and
#' matched control sites. This file can be used to quantify per-cell activity
#' scores with the [Bacha()] function.
#'
#' By default, controls are matched for GC content and sequence entropy.
#'
#' # Index file format
#'
#' The output is a BGZF (Blocked GZIP) compressed binary file with the
#' following layout. BGZF is a sequence of independent gzip blocks (each up
#' to 64 KB), which bacha decompresses in parallel for fast index loading.
#' Writing the file requires the `bgzip` command from htslib to be on the
#' user's PATH.
#'
#' **Magic and version**
#'
#' - 4 bytes: magic string `"BACH"`
#' - 1 byte (u8): format version (currently 1)
#'
#' **Header**
#'
#' - Length-prefixed string (u16 LE length + UTF-8 bytes): genome name
#' - Length-prefixed string: motif set name
#' - 4 bytes (u32 LE): window size (bp)
#' - 4 bytes (u32 LE): number of TFs (`n_tfs`)
#' - 4 bytes (u32 LE): number of chromosomes (`n_chroms`)
#'
#' **Build provenance**
#'
#' - 2 bytes (u16 LE): number of entries (may be 0)
#' - For each entry: a length-prefixed key string then a length-prefixed value
#'   string. Values are stored as text. `BachaIndex()` writes `build_date`
#'   (UTC, ISO 8601), `bacha_version`, `r_version`, `input` (`"pwm"` or
#'   `"positions"`), `match_features`, `n_background`, `mindist`,
#'   `max_matches` and `p_cutoff` (motif scanning only), `gc_jitter_sd`,
#'   `seed`, `exclude`, `bias` (file name, when used), the number of
#'   `open_regions` and `mask_regions` supplied, and `genome_version` (the
#'   BSgenome package version, when installed), followed by any `metadata`
#'   entries. Readers ignore keys they do not recognise, so entries can be
#'   added without changing the format version. `bacha view` and
#'   `bacha summary` print them in a `##[build]` section.
#'
#' **Feature table**
#'
#' - For each of the `n_tfs` features: a length-prefixed ID string, then a
#'   length-prefixed display name string. IDs are unique. For PWMs from
#'   TFBSTools (including `Signac::ReadPWM()`) they are the matrix IDs and the
#'   names are the matrix names; otherwise both come from the list names.
#'   Names may repeat; `bacha count --feature-names` picks which is used to
#'   label the output.
#'
#' **Instance counts**
#'
#' - `[u32 LE; n_tfs]`: number of motif site instances per TF
#' - `[u32 LE; n_tfs]`: number of control instances per TF
#'
#' **Chromosome table**
#'
#' - For each chromosome: length-prefixed string (name) + u32 LE (length in bp)
#'
#' **Motif matrices**
#'
#' - 1 byte (u8): 1 if the index stores matrices, 0 if not (indexes built
#'   from `positions`). When 1, for each feature:
#'   - 1 byte (u8): value type, 0 = no matrix, 1 = counts, 2 = probabilities,
#'     3 = log-odds (e.g. a TFBSTools `PWMatrix`), 4 = other
#'   - 2 bytes (u16 LE): width (number of positions; 0 when there is no
#'     matrix)
#'   - `[f32 LE; 4 * width]`: values position by position, A, C, G, T within
#'     each position
#'
#' Export them with `bacha motifs` or [IndexMotifs()].
#'
#' **Varint columns**
#'
#' Most position columns are stored as unsigned LEB128 varints: each value
#' takes 1 to 5 bytes, 7 bits per byte, least significant group first, with
#' the high bit of a byte set when another byte follows. Values never exceed
#' 32 bits. A varint column starts with 8 bytes (u64 LE) giving its length in
#' bytes, followed by exactly that many bytes encoding the column's values;
#' the number of values comes from the block's row count.
#'
#' **Site positions** (sorted by chromosome then start)
#'
#' - 4 bytes (u32 LE): number of site rows (`n_sites`)
#' - varint column, `n_sites` values: chromosome IDs (0-based index into the
#'   chromosome table)
#' - varint column, `n_sites` values: start positions (1-based),
#'   delta-encoded within each chromosome run. The first start of each
#'   chromosome is absolute; subsequent starts within the same chromosome
#'   are differences from the previous start. Reader recovers absolute
#'   starts via cumulative sum reset at each chromosome boundary, then
#'   converts to 0-based half-open intervals for fragment endpoint counting.
#' - varint column, `n_sites` values: TF indices (1-based index into TF
#'   names)
#'
#' **Control positions** (sorted by chromosome then start)
#'
#' - 4 bytes (u32 LE): number of control rows (`n_ctrl`)
#' - varint column, `n_ctrl` values: chromosome IDs
#' - varint column, `n_ctrl` values: 1-based start positions, delta-encoded
#'   within each chromosome run (same scheme as site starts).
#' - varint column, `n_ctrl` values: number of TF associations per control
#'   (`n_per_row`)
#' - 4 bytes (u32 LE): total number of TF index entries (`total_indices`)
#' - varint column, `total_indices` values: flattened TF indices for each
#'   control row, delta-encoded within each row. Within a row the indices
#'   are sorted ascending; the first index of each row is absolute and
#'   subsequent indices are differences from the previous index. Each
#'   absolute index is encoded as `tf_index + n_tfs` to distinguish controls
#'   from sites in the shared ID space.
#'
#' **Trailing data**
#'
#' Readers stop after the control positions block and ignore anything that
#' follows it. Optional sections can therefore be appended in future without
#' changing the format version: files carrying them still load in older
#' readers, which skip them.
#'
#' @param pwm A list of position weight matrices or position frequency matrices.
#'   If provided, the `positions` argument should not be set.
#' @param genome A BSgenome object.
#' @param bias Optional path to a H5 file containing Tn5 insertion bias scores
#'   for each base in the genome. Required only when `match_features` includes
#'   `"meanBias"` or `"maxBias"`. See <https://zenodo.org/records/15224770>;
#'   <https://doi.org/10.1038/s41586-024-08443-4>.
#' @param open_regions A GRanges object containing open chromatin regions.
#'   Strand is ignored.
#' @param mask_regions A GRanges object containing regions to mask from the
#'   motif scan (e.g. blacklist regions or repeat regions).
#' @param output_path Path to write the output index file.
#' @param positions A [GenomicRanges::GRangesList] object containing
#'   pre-computed genomic positions of interest for a set of features. Each
#'   element should be a different feature (e.g. a TF motif or GWAS trait) and
#'   the genomic regions correspond to regions relevant to that feature to
#'   aggregate information across. Strand is ignored. If `positions` are
#'   supplied, `pwm` should
#'   not be provided and the motif scanning will be skipped.
#' @param tss_positions A GRanges object giving genomic TSS positions. Required
#'   if `"tssDist"` is included in `match_features`.
#' @param exclude Vector of chromosome names to exclude.
#' @param match_features DNA sequence features to control for. Any combination
#'   of `"GC"` (GC content), `"CpG"` (CpG dinucleotide frequency),
#'   `"entropy"` (Shannon entropy of A/C/G/T frequencies; captures
#'   residual-repeat / low-complexity context),
#'   `"meanBias"` (mean Tn5 insertion bias), `"maxBias"` (maximum Tn5 insertion
#'   bias), and `"tssDist"` (log-transformed distance to the nearest TSS), or
#'   `NULL` to randomly select control regions.
#' @param window Number of bases for each flanking region upstream and downstream
#'   of the motif.
#' @param max_matches Maximum number of motif matches to retain per motif.
#'   When a motif has more than this many matches, only the top matches by
#'   score are kept. Prevents GC-rich motifs from dominating with weak hits.
#' @param mindist Minimum distance between a control region and an instance of
#'   the query motif.
#' @param n_background Number of control regions to select per motif match.
#' @param motif_set_name Character string identifying the motif collection
#'   (e.g. `"JASPAR2026"`, `"HOCOMOCO14"`). Stored in the index file header for
#'   provenance tracking. If `NULL`, defaults to `"custom"`.
#' @param metadata Optional named character vector (or named list of single
#'   values) of extra provenance to store in the index's build block, such as
#'   `c(motif_source = "https://jaspar.elixir.no/...")`. Names must be unique
#'   and must not reuse a key `BachaIndex()` records itself; neither names nor
#'   values may contain tabs or newlines. Shown by `bacha view`,
#'   `bacha summary` and [IndexSummary()].
#' @param p_cutoff p-value cutoff for motif matches.
#' @param gc_jitter_sd Standard deviation of Gaussian jitter applied to GC
#'   values during matching only (not stored in the index). For an `N`-bp
#'   window, GC takes at most `N + 1` discrete values, so many candidates
#'   share each site's GC exactly. With unbroken ties the k-NN matcher
#'   returns the same handful of candidates for every site at a given GC
#'   value and the greedy selector exhausts them after picking a few
#'   thousand controls, regardless of how many were requested. A
#'   sub-integer jitter breaks ties while keeping the integer-rounded GC
#'   bin as the dominant matching signal. Set to `0` to disable.
#' @param seed Optional integer seed for index construction. Control selection
#'   draws on the random number generator (GC jitter, and random sampling
#'   when `match_features` is `NULL`), so two runs with the same inputs give
#'   different indexes unless a seed is fixed. The seed is applied with
#'   [withr::local_seed()], so the caller's random number state is restored
#'   when `BachaIndex()` returns.
#' @param verbose Display messages.
#' 
#' @importFrom GenomicRanges resize granges reduce slidingWindows distanceToNearest flank trim
#' @importFrom IRanges subsetByOverlaps Views viewMeans viewMaxs findOverlaps ranges overlapsAny
#' @importFrom hdf5r H5File
#' @importFrom motifmatchr matchMotifs
#' @importFrom utils txtProgressBar setTxtProgressBar head
#' @importFrom S4Vectors mcols mcols<- Rle endoapply queryHits
#' @importFrom GenomeInfoDb seqlevels seqlevelsInUse keepSeqlevels dropSeqlevels
#' @importFrom stats cov median rnorm setNames
#' @importFrom BiocGenerics start<- end unstrand
#' @importFrom withr local_seed
#' 
#' @export
#' @examples
#' \dontrun{
#' library(GenomicRanges)
#' library(AnnotationHub)
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' library(Signac)
#' library(TFBSTools)
#' 
#' # load open chromatin regions
#' ocr <- read.table("cre.bed")
#' ocr <- makeGRangesFromDataFrame(ocr, seqnames.field = "V1", start.field = "V2", end.field = "V3")
#' 
#' # load regions to exclude
#' ah <- AnnotationHub()
#' excludeGR.hg38 <- ah[["AH107305"]]
#' repeat_regions <- ah[['AH99003']]
#' mask_regions <- c(excludeGR.hg38, repeat_regions)
#' mask_regions <- reduce(mask_regions)
#' strand(mask_regions) <- "*"
#' 
#' # load PWMs
#' pwm <- ReadJASPAR("JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt",
#'                   mode = "jaspar", keep = "all")
#' 
#' # create Bacha index
#' BachaIndex(
#'   pwm = pwm,
#'   genome = BSgenome.Hsapiens.UCSC.hg38,
#'   open_regions = ocr,
#'   mask_regions = mask_regions,
#'   output_path = "hg38.j24.v1.bi"
#' )
#' }
#' 
BachaIndex <- function(
    pwm,
    genome,
    open_regions,
    mask_regions,
    output_path,
    positions,
    tss_positions,
    exclude = c("chrM", "chrY", "chrX"),
    match_features = c("GC", "entropy"),
    bias = NULL,
    window = 100,
    max_matches = 50000,
    mindist = 1000,
    n_background = 10,
    p_cutoff = 5e-05,
    motif_set_name = NULL,
    gc_jitter_sd = 0.5,
    seed = 1234,
    metadata = NULL,
    verbose = TRUE
) {
  # Seed for the duration of this call only; the caller's RNG state is
  # restored on exit rather than being clobbered.
  if (!is.null(x = seed)) {
    local_seed(seed = seed)
  }
  window <- as.integer(x = window)
  if (is.na(x = window) || window < 2L) {
    stop("'window' must be an integer >= 2")
  }
  valid_match_features <- c("GC", "CpG", "entropy", "meanBias", "maxBias", "tssDist")
  invalid_match_features <- setdiff(x = match_features, y = valid_match_features)
  if (length(x = invalid_match_features) > 0) {
    stop(
      "Invalid match feature(s): ",
      paste(invalid_match_features, collapse = ", "),
      ". Valid options are: ",
      paste(valid_match_features, collapse = ", ")
    )
  }
  bias_missing <- is.null(x = bias) ||
    length(x = bias) != 1L ||
    is.na(x = bias) ||
    !nzchar(x = bias)
  if (any(c("meanBias", "maxBias") %in% match_features) && bias_missing) {
    stop(
      "'bias' must be provided when 'match_features' includes ",
      "'meanBias' or 'maxBias'"
    )
  }
  if (!bias_missing && !any(c("meanBias", "maxBias") %in% match_features)) {
    warning(
      "'bias' was supplied but 'match_features' includes neither 'meanBias' ",
      "nor 'maxBias', so the bias file will be ignored. Add \"meanBias\" ",
      "(and/or \"maxBias\") to 'match_features' to match on Tn5 bias.",
      call. = FALSE
    )
  }
  # mutually exclusive arguments
  if (missing(x = pwm) && missing(x = positions)) {
    stop("Either 'pwm' or 'positions' must be provided")
  }
  if (!missing(x = pwm) && !missing(x = positions)) {
    stop("'pwm' and 'positions' are mutually exclusive. Provide only one.")
  }
  
  # first check that the output file path exists
  if (!dir.exists(paths = dirname(path = output_path))) {
    stop("Requested output path does not exist")
  }
  if (file.exists(output_path)) {
    stop("Requested output file already exists")
  }
  
  # Record how this index was built, before the inputs below are filtered.
  build_info <- index_build_info(
    input = if (missing(x = pwm)) "positions" else "pwm",
    match_features = match_features,
    n_background = n_background,
    mindist = mindist,
    max_matches = if (missing(x = pwm)) NULL else max_matches,
    p_cutoff = if (missing(x = pwm)) NULL else p_cutoff,
    gc_jitter_sd = gc_jitter_sd,
    seed = seed,
    exclude = exclude,
    # recorded only when matching uses it; an unused bias file is ignored
    bias = if (bias_missing || !any(c("meanBias", "maxBias") %in% match_features)) {
      NULL
    } else {
      bias
    },
    n_open_regions = length(x = open_regions),
    n_mask_regions = if (missing(x = mask_regions)) NULL else length(x = mask_regions),
    genome = genome
  )
  build_info <- add_build_metadata(build_info = build_info, metadata = metadata)

  # An index is unstranded: sites, flanks and controls are plain positions,
  # and the writer delta-encodes starts sorted by position alone. Stranded
  # input would otherwise sort as all `+` then all `-` within a chromosome.
  open_regions <- unstrand(x = open_regions)

  # Expand mask regions by `window` bp on each side. All motif sites, motif
  # flanks, and candidate controls are filtered against this padded mask.
  if (!missing(mask_regions)) {
    mask_regions <- resize(
      x = mask_regions,
      width = width(x = mask_regions) + 2L * window,
      fix = "center"
    )
    mask_regions <- trim(x = mask_regions)
  }

  # remove open regions overlapping mask regions
  if (!missing(mask_regions)) {
    open_regions <- subsetByOverlaps(
      x = open_regions, ranges = mask_regions, invert = TRUE
    )
  }
  
  # remove open regions on exclusion chromosomes
  open_regions <- open_regions[
    !(as.character(seqnames(x = open_regions)) %in% exclude)
  ]
  open_regions <- dropSeqlevels(open_regions, exclude, pruning.mode = "coarse")

  if ("tssDist" %in% match_features) {
    if (missing(tss_positions)) {
      stop("'tss_positions' must be provided when 'tssDist' is included in 'match_features'")
    }
    if (!inherits(tss_positions, "GRanges")) {
      stop("'tss_positions' must be a GRanges object")
    }
    tss_positions <- tss_positions[
      !(as.character(seqnames(x = tss_positions)) %in% exclude)
    ]
    tss_positions <- keepSeqlevels(
      x = tss_positions,
      value = intersect(seqlevelsInUse(x = tss_positions), seqlevels(open_regions)),
      pruning.mode = "coarse"
    )
  }
  
  if (!missing(pwm)) {
    # Features are keyed by a unique ID (the matrix ID where the PWM carries
    # one, else the list name) and carry a display name and the matrix
    # itself, which are stored in the index alongside the sites.
    motif_info <- motif_metadata(pwm = pwm)
    names(x = pwm) <- motif_info$id
    feature_names <- motif_info$name
    feature_matrices <- motif_info$matrix
    if (verbose) {
      message("Scanning genome for motif matches")
    }
    motif.matches <- matchMotifs(
      pwms = pwm,
      subject = open_regions,
      out = 'positions',
      genome = genome,
      p.cutoff = p_cutoff
    )

    # score-based cap: keep only top max_matches per motif
    motif.matches <- endoapply(motif.matches, function(gr) {
      if (length(gr) > max_matches) {
        gr <- gr[order(gr$score, decreasing = TRUE)[seq_len(max_matches)]]
      }
      gr
    })

    # per-motif match width, captured before reduce() merges adjacent matches
    motif_widths <- vapply(
      X = as.list(x = width(x = motif.matches)),
      FUN = function(w) {
        if (length(x = w) == 0L) NA_integer_ else as.integer(x = median(x = w))
      },
      FUN.VALUE = integer(1)
    )
    names(x = motif_widths) <- names(x = motif.matches)

    # merge nearby matches (within 50bp) to handle dimeric binding with spacers
    merge_gap <- 50L
    motif.matches <- reduce(
      x = motif.matches, min.gapwidth = merge_gap, ignore.strand = TRUE
    )

    # unlist and add motif name
    motif.names <- rep(
      x = names(x = motif.matches), elementNROWS(x = motif.matches)
    )
    motif.sites <- unlist(x = motif.matches, use.names = FALSE)
    motif.sites$motif <- motif.names

    # remove merged ranges that are too wide (dense clusters in
    # extreme-composition regions that produce biased flanking regions). The
    # cap is per motif, using that motif's OWN width W: keep a single match
    # (width W) or a dimeric merge of two instances within the gap
    # (width <= 2*W + merge_gap); anything wider implies >2 stacked instances.
    width_cap <- 2L * motif_widths[motif.sites$motif] + merge_gap
    width_cap[is.na(x = width_cap)] <- .Machine$integer.max
    motif.sites <- motif.sites[width(x = motif.sites) <= width_cap]

    # remove motif sites overlapping mask regions
    if (!missing(mask_regions)) {
      motif.sites <- subsetByOverlaps(
        x = motif.sites, ranges = mask_regions, invert = TRUE
      )
    }

    # create upstream and downstream flanking regions
    n_sites <- length(motif.sites)
    upstream <- flank(x = motif.sites, width = window, start = TRUE)
    downstream <- flank(x = motif.sites, width = window, start = FALSE)
    upstream$motif <- motif.sites$motif
    downstream$motif <- motif.sites$motif
    upstream$site_id <- seq_len(n_sites)
    downstream$site_id <- seq_len(n_sites)
    motif.positions <- c(upstream, downstream)
    # drop regions that overhang chromosome boundaries
    sl <- seqlengths(motif.positions)[as.character(seqnames(motif.positions))]
    sl[is.na(sl)] <- Inf
    motif.positions <- motif.positions[
      start(x = motif.positions) >= 1 & end(x = motif.positions) <= sl
    ]

    # remove flanks overlapping mask regions
    if (!missing(mask_regions)) {
      motif.positions <- subsetByOverlaps(
        x = motif.positions, ranges = mask_regions, invert = TRUE
      )
    }

    # remove flanks that overlap a motif site of the same TF
    if (verbose) {
      message("Removing flanks overlapping their own motif sites")
    }
    flank_motif_f <- factor(motif.positions$motif)
    site_motif_f <- factor(motif.sites$motif, levels = levels(flank_motif_f))
    flank_idx_by_motif <- split(seq_along(motif.positions), flank_motif_f)
    site_idx_by_motif <- split(seq_along(motif.sites), site_motif_f)
    has_own_motif <- logical(length(motif.positions))
    for (m in levels(flank_motif_f)) {
      f_idx <- flank_idx_by_motif[[m]]
      s_idx <- site_idx_by_motif[[m]]
      if (length(f_idx) == 0 || length(s_idx) == 0) next
      hit <- overlapsAny(motif.positions[f_idx], motif.sites[s_idx])
      has_own_motif[f_idx[hit]] <- TRUE
    }
    if (verbose && any(has_own_motif)) {
      pct <- round(100 * sum(has_own_motif) / length(has_own_motif), 1)
      message(
        "  Removed ", sum(has_own_motif),
        " flanks (", pct, "%) overlapping their own motif"
      )
    }
    motif.positions <- motif.positions[!has_own_motif]

    # keep only sites where both flanks survived filtering
    paired_ids <- names(which(table(motif.positions$site_id) == 2L))
    motif.positions <- motif.positions[motif.positions$site_id %in% paired_ids]
    motif.sites <- motif.sites[as.integer(paired_ids)]
    motif.positions <- interleave_flanks(x = motif.positions)
    motif.positions$site_id <- NULL

    motif_set <- unique(x = motif.positions$motif)
  } else {
    if (verbose) {
      message("Using provided genomic positions")
    }
    if (!inherits(positions, "GRangesList")) {
      stop("'positions' must be a GRangesList object")
    }
    names(x = positions) <- normalize_feature_names(
      x = names(x = positions),
      n = length(x = positions)
    )
    feature_names <- setNames(object = names(x = positions), nm = names(x = positions))
    feature_matrices <- NULL
    motif.names <- rep(
      x = names(x = positions), elementNROWS(x = positions)
    )
    positions <- unstrand(x = unlist(x = positions, use.names = FALSE))
    positions$motif <- motif.names
    # remove provided sites overlapping padded mask regions
    if (!missing(mask_regions)) {
      positions <- subsetByOverlaps(
        x = positions, ranges = mask_regions, invert = TRUE
      )
    }
    # create upstream and downstream flanking regions
    n_sites <- length(positions)
    upstream <- flank(x = positions, width = window, start = TRUE)
    downstream <- flank(x = positions, width = window, start = FALSE)
    upstream$motif <- positions$motif
    downstream$motif <- positions$motif
    upstream$site_id <- seq_len(n_sites)
    downstream$site_id <- seq_len(n_sites)
    motif.positions <- c(upstream, downstream)
    # drop regions that overhang chromosome boundaries
    sl <- seqlengths(motif.positions)[as.character(seqnames(motif.positions))]
    sl[is.na(sl)] <- Inf
    motif.positions <- motif.positions[start(x = motif.positions) >= 1 & end(x = motif.positions) <= sl]
    motif.positions <- subsetByOverlaps(
      x = motif.positions,
      ranges = open_regions,
    )
    # remove flanks overlapping mask regions
    if (!missing(mask_regions)) {
      motif.positions <- subsetByOverlaps(
        x = motif.positions, ranges = mask_regions, invert = TRUE
      )
    }
    # keep only sites where both flanks survived filtering
    paired_ids <- names(which(table(motif.positions$site_id) == 2L))
    motif.positions <- motif.positions[motif.positions$site_id %in% paired_ids]
    motif.sites <- positions[as.integer(paired_ids)]
    motif.positions <- interleave_flanks(x = motif.positions)
    motif.positions$site_id <- NULL
    motif_set <- unique(x = motif.positions$motif)
  }
  n_motif_total <- length(x = motif_set)

  if (verbose) {
    message("Generating control regions")
  }
  # generate background regions
  # place candidate points spaced `window` bp apart within each OCR to simulate
  # motif-sized regions, then take upstream and downstream flanks to match how
  # motif flank regions are determined
  extended_ocr <- resize(
    x = open_regions,
    width = width(x = open_regions) + 2L * window,
    fix = "center"
  )
  extended_ocr <- trim(x = extended_ocr)
  # An OCR trimmed at a chromosome edge can end up shorter than one window,
  # which has no valid full-width candidate inside it.
  extended_ocr <- extended_ocr[width(x = extended_ocr) >= window]
  candidate_bg <- unlist(x = slidingWindows(
    x = extended_ocr,
    width = window,
    step = window %/% 2L
  ))
  short <- width(x = candidate_bg) < window
  start(x = candidate_bg)[short] <- end(x = candidate_bg)[short] - window + 1L
  # drop regions that overhang chromosome boundaries
  sl <- seqlengths(candidate_bg)[as.character(seqnames(candidate_bg))]
  sl[is.na(sl)] <- Inf
  candidate_bg <- candidate_bg[start(x = candidate_bg) >= 1 & end(x = candidate_bg) <= sl]
  # remove controls overlapping padded mask regions; `extended_ocr` can
  # otherwise place controls back near masks after OCR filtering.
  if (!missing(mask_regions)) {
    candidate_bg <- subsetByOverlaps(
      x = candidate_bg, ranges = mask_regions, invert = TRUE
    )
  }
  candidate_bg <- keepSeqlevels(
    x = candidate_bg,
    value = seqlevelsInUse(x = candidate_bg),
    pruning.mode = "coarse"
  )
  
  motif.positions <- keepSeqlevels(
    x = motif.positions,
    value = seqlevelsInUse(x = motif.positions),
    pruning.mode = "coarse"
  )
  if (verbose) {
    message("Created ", format(x = length(x = candidate_bg), big.mark = ","),
            " control regions")
  }
  
  if ("GC" %in% match_features) {
    message("Computing genomic region GC content")
    candidate_bg <- add_gc(regions = candidate_bg, genome = genome)
    motif.positions <- add_gc(regions = motif.positions, genome = genome)
  }

  if ("CpG" %in% match_features) {
    message("Computing genomic region CpG dinucleotide frequency")
    candidate_bg <- add_cpg(regions = candidate_bg, genome = genome)
    motif.positions <- add_cpg(regions = motif.positions, genome = genome)
  }

  if ("entropy" %in% match_features) {
    message("Computing genomic region Shannon entropy")
    candidate_bg <- add_entropy(regions = candidate_bg, genome = genome)
    motif.positions <- add_entropy(regions = motif.positions, genome = genome)
  }

  if ("tssDist" %in% match_features) {
    message("Computing distance to nearest TSS")
    candidate_bg <- add_tss_distance(regions = candidate_bg, tss_positions = tss_positions)
    motif.positions <- add_tss_distance(regions = motif.positions, tss_positions = tss_positions)
  }
  
  if (("meanBias" %in% match_features) | ("maxBias" %in% match_features)) {
    message("Computing genomic region Tn5 bias")
    
    # load bias file
    # each element is a chromosome
    # position in vector corresponds to position in genome
    bias <- H5File$new(bias, mode = "r")
    on.exit(expr = bias$close(), add = TRUE)
    
    # Cover every chromosome carrying a control or a motif flank. Iterating
    # only the control seqlevels leaves motif regions on any other chromosome
    # with NA bias, which then propagates into the matching matrix.
    chromosomes <- union(
      x = seqlevelsInUse(x = candidate_bg),
      y = seqlevelsInUse(x = motif.positions)
    )
    has_bias <- vapply(
      X = chromosomes,
      FUN = function(ch) bias$exists(ch),
      FUN.VALUE = logical(1)
    )
    if (!all(has_bias)) {
      missing_bias <- chromosomes[!has_bias]
      stop(
        "The bias file has no data for chromosome(s): ",
        paste(head(x = missing_bias, n = 5L), collapse = ", "),
        if (length(x = missing_bias) > 5L) ", ..." else "",
        ". Use a bias file built for this genome, or drop these sequences ",
        "with 'exclude'.",
        call. = FALSE
      )
    }

    # Accumulate into plain vectors and attach once at the end. Assigning
    # through `gr[mask]$col <- value` goes through R's nested-replacement
    # machinery, copying the whole GRanges on every assignment.
    bg_seqnames <- as.character(x = seqnames(x = candidate_bg))
    mp_seqnames <- as.character(x = seqnames(x = motif.positions))
    bg_mean <- rep_len(x = NA_real_, length.out = length(x = candidate_bg))
    bg_max <- bg_mean
    mp_mean <- rep_len(x = NA_real_, length.out = length(x = motif.positions))
    mp_max <- mp_mean

    show_pb <- verbose && (length(x = chromosomes) > 1) && (n_motif_total > 1)
    if (show_pb) {
      pb <- txtProgressBar(
        min = 0,
        max = length(x = chromosomes),
        style = 3,
        file = stderr()
      )
      on.exit(expr = close(pb), add = TRUE)
    }
    invisible(x = gc())
    for (i in seq_along(along.with = chromosomes)) {
      ch <- chromosomes[i]
      b <- Rle(values = as.numeric(x = bias[[ch]][]))

      bg_idx <- which(bg_seqnames == ch)
      if (length(x = bg_idx) > 0) {
        v <- Views(subject = b, start = ranges(x = candidate_bg[bg_idx]))
        bg_mean[bg_idx] <- viewMeans(x = v)
        bg_max[bg_idx] <- viewMaxs(x = v)
      }

      mp_idx <- which(mp_seqnames == ch)
      if (length(x = mp_idx) > 0) {
        v <- Views(subject = b, start = ranges(x = motif.positions[mp_idx]))
        mp_mean[mp_idx] <- viewMeans(x = v)
        mp_max[mp_idx] <- viewMaxs(x = v)
      }

      if (show_pb) {
        setTxtProgressBar(pb, i)
      }
      invisible(x = gc(full = FALSE))
    }
    candidate_bg$meanBias <- bg_mean
    candidate_bg$maxBias <- bg_max
    motif.positions$meanBias <- mp_mean
    motif.positions$maxBias <- mp_max
  }
  
  # Break GC ties so the k-NN matcher returns distinct candidates per site
  if (("GC" %in% match_features) && gc_jitter_sd > 0) {
    candidate_bg$GC    <- candidate_bg$GC +
      rnorm(n = length(x = candidate_bg),    sd = gc_jitter_sd)
    motif.positions$GC <- motif.positions$GC +
      rnorm(n = length(x = motif.positions), sd = gc_jitter_sd)
  }

  if (!is.null(x = match_features)) {
    # convert to uncorrelated variables
    mat <- as.matrix(x = mcols(x = candidate_bg)[, match_features, drop = FALSE])
    # Validate both matrices up front. A single NA or non-finite value makes
    # cov() non-finite and chol() fail with an opaque LAPACK error, and a
    # non-finite motif value silently poisons the k-NN query far downstream.
    check_match_features(x = mat, what = "control", features = match_features)
    check_match_features(
      x = as.matrix(x = mcols(x = motif.positions)[, match_features, drop = FALSE]),
      what = "motif",
      features = match_features
    )
    chol_cov_mat <- tryCatch(
      expr = chol(x = cov(x = mat)),
      error = function(e) {
        stop(
          "Could not decompose the covariance matrix of the match features (",
          paste(match_features, collapse = ", "), "). They are likely ",
          "collinear; drop one of them from 'match_features'. Original error: ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
    trans_bg <- t(forwardsolve(t(chol_cov_mat), t(mat)))
    colnames(x = trans_bg) <- match_features
    invisible(x = gc())
  }

  if (verbose) {
    message("Matching control regions for each motif")
  }

  # iterate over motifs
  # For each feature, store the control rows selected for it
  n_bg <- length(x = candidate_bg)
  selected_list <- vector("list", length = n_motif_total)

  # Candidate controls are sliding windows that overlap their neighbours by
  # half a window, and padded OCRs can overlap each other, so two controls
  # picked for the same motif could share sequence and count the same
  # fragments twice. Sorting the candidates by position lets the controls
  # overlapping any one of them be found by walking outwards from it.
  bg_order <- order(
    as.integer(x = seqnames(x = candidate_bg)),
    start(x = candidate_bg)
  )
  bg_rank <- integer(length = n_bg)
  bg_rank[bg_order] <- seq_len(length.out = n_bg)
  bg_sorted_chrom <- as.integer(x = seqnames(x = candidate_bg))[bg_order]
  bg_sorted_start <- start(x = candidate_bg)[bg_order]
  bg_sorted_end <- end(x = candidate_bg)[bg_order]
  bg_max_width <- max(width(x = candidate_bg))

  if (n_motif_total > 1) {
    pb <- txtProgressBar(
      min = 1,
      max = n_motif_total,
      style = 3,
      file = stderr()
    )
  }

  for (i in seq_along(motif_set)) {
    motif_match <- motif.positions[motif.positions$motif == motif_set[i]]

    # exclude controls that overlap motif sites or flanks
    exclusion_regions <- c(
      granges(motif.sites[motif.sites$motif == motif_set[i]]),
      granges(motif_match)
    )
    olap <- findOverlaps(
      query = candidate_bg,
      subject = exclusion_regions,
      maxgap = mindist
    )
    excluded <- queryHits(x = olap)

    n_eligible_total <- n_bg - length(unique(excluded))
    n_requested <- length(x = motif_match) * n_background
    
    if (is.null(x = match_features)) {
      # random selection (no feature matching)
      n <- min(n_requested, n_eligible_total)
      if (n_eligible_total == 0) {
        warning("No eligible control regions for ", motif_set[i])
        selected <- integer(0)
      } else {
        used_controls <- rep(FALSE, n_bg)
        used_controls[excluded] <- TRUE
        selected <- integer(length = n)
        n_selected <- 0L
        # draw candidates in batches and keep each one that does not overlap a
        # control already kept for this motif, until enough have been found or
        # the pool of eligible candidates runs out
        pool <- which(x = !used_controls)
        while (n_selected < n && length(x = pool) > 0) {
          take <- sample.int(
            n = length(x = pool),
            size = min(length(x = pool), max(1000L, 2L * (n - n_selected)))
          )
          batch <- pool[take]
          pool <- pool[-take]
          for (ctrl in batch) {
            if (used_controls[ctrl]) next
            n_selected <- n_selected + 1L
            selected[n_selected] <- ctrl
            used_controls[overlapping_controls(
              index = bg_rank[ctrl],
              chrom = bg_sorted_chrom,
              starts = bg_sorted_start,
              ends = bg_sorted_end,
              max_width = bg_max_width,
              order = bg_order
            )] <- TRUE
            if (n_selected >= n) break
          }
        }
        selected <- if (n_selected == 0L) integer(0) else selected[seq_len(n_selected)]
      }
    } else {
      # per-instance k-nearest-neighbor matching using KD-tree
      trans_motif <- t(
        forwardsolve(t(chol_cov_mat), t(as.matrix(mcols(motif_match)[, match_features])))
      )
      eligible <- rep(TRUE, n_bg)
      if (length(excluded) > 0) eligible[excluded] <- FALSE
      eligible_idx <- which(eligible)
      bg_eligible <- trans_bg[eligible_idx, , drop = FALSE]

      n_cap <- as.integer(n_eligible_total * 0.2)
      n_instances <- nrow(trans_motif)
      k_search <- min(n_background, nrow(bg_eligible), max(1L, n_cap %/% n_instances))
      if (k_search < n_background && verbose) {
        message(sprintf("  %s: reducing k from %d to %d (cap %d controls for %d instances)",
                        motif_set[i], n_background, k_search, n_cap, n_instances))
      }
      if (k_search > 0) {
        k_query <- min(
          nrow(bg_eligible),
          max(
            as.integer(3L * k_search),
            20L
          )
        )
        nn_result <- RANN::nn2(
          data = bg_eligible,
          query = trans_motif,
          k = k_query
        )

        target_total <- min(length(eligible_idx), n_instances * k_search)
        selected <- integer(target_total)
        n_selected <- 0L
        used_controls <- rep(FALSE, n_bg)

        for (inst in seq_len(n_instances)) {
          inst_neighbors <- eligible_idx[nn_result$nn.idx[inst, ]]
          inst_count <- 0L

          for (ctrl in inst_neighbors) {
            if (n_selected >= target_total) break
            if (!used_controls[ctrl]) {
              n_selected <- n_selected + 1L
              selected[n_selected] <- ctrl
              # block this control and any overlapping it for this motif
              used_controls[overlapping_controls(
                index = bg_rank[ctrl],
                chrom = bg_sorted_chrom,
                starts = bg_sorted_start,
                ends = bg_sorted_end,
                max_width = bg_max_width,
                order = bg_order
              )] <- TRUE
              inst_count <- inst_count + 1L
              if (inst_count >= k_search) break
            }
          }

          if (n_selected >= target_total) {
            break
          }
        }
        if (n_selected == 0L) {
          selected <- integer(0)
        } else {
          selected <- selected[seq_len(n_selected)]
        }
      } else {
        selected <- integer(0)
      }
    }

    # record the control rows selected for this feature
    selected_list[[i]] <- selected
    if (verbose & (n_motif_total > 1)) {
      setTxtProgressBar(pb, i)
    }
    if (i %% 50 == 0) {
      invisible(x = gc(full = FALSE))
    }
  }

  # controls per feature, counted before the selections are inverted (each
  # `selected` holds distinct control rows, so this is just its length)
  n_controls <- setNames(lengths(x = selected_list), motif_set)

  # transpose to per-control feature lists, dropping controls no feature chose
  inverted <- invert_selections(
    selected_list = selected_list,
    n_bg = n_bg,
    id_offset = n_motif_total
  )
  rm(selected_list)
  motifs_list <- inverted$motif_id
  candidate_bg <- candidate_bg[inverted$used]
  rm(inverted)
  invisible(x = gc(full = FALSE))

  # get number of regions for each id
  n_motif_sites <- table(motif.positions$motif)[motif_set]

  if (verbose) {
    message("Writing results to index file")
  }

  # write results to index file
  write_index(
    motifs = motif.positions,
    controls = candidate_bg,
    motif_id = motifs_list,
    filename = output_path,
    n_motif = n_motif_sites,
    n_controls = n_controls,
    genome = genome,
    extra_chroms = union(x = exclude, y = seqlevels(x = open_regions)),
    motif_set_name = motif_set_name %||% "custom",
    build_info = build_info,
    feature_names = feature_names,
    matrices = feature_matrices
  )
}

# Matrix value types, in the order of their codes in the index (0-based).
MATRIX_TYPES <- c("none", "counts", "probabilities", "logodds", "other")

# Per-feature metadata for a set of PWMs: a unique ID, a display name and the
# matrix (A/C/G/T rows) with its value type, each named by ID.
#
# TFBSTools matrices supply their own ID and name (read from the slots, so
# TFBSTools need not be attached); for anything else both come from the list
# names. A missing display name falls back to the ID. IDs are made unique
# with `normalize_feature_names()`; display names may repeat.
#' @importFrom methods is slot
motif_metadata <- function(pwm) {
  elems <- as.list(x = pwm)
  n <- length(x = elems)
  list_names <- names(x = pwm) %||% rep(x = NA_character_, times = n)
  from_slot <- function(x, what) {
    if (is(object = x, class2 = "XMatrix")) {
      v <- slot(object = x, name = what)
      if (length(x = v) == 1L && !is.na(x = v) && nzchar(x = v)) {
        return(as.character(x = v))
      }
    }
    NA_character_
  }
  ids <- vapply(X = elems, FUN = from_slot, FUN.VALUE = "", what = "ID")
  ids[is.na(x = ids)] <- list_names[is.na(x = ids)]
  ids <- normalize_feature_names(x = unname(ids), n = n)
  disp <- vapply(X = elems, FUN = from_slot, FUN.VALUE = "", what = "name")
  disp[is.na(x = disp)] <- list_names[is.na(x = disp)]
  disp <- unname(disp)
  disp[is.na(x = disp) | !nzchar(x = disp)] <- ids[is.na(x = disp) | !nzchar(x = disp)]
  mats <- lapply(X = elems, FUN = motif_matrix)
  list(
    id = ids,
    name = setNames(object = disp, nm = ids),
    matrix = setNames(object = mats, nm = ids)
  )
}

# A motif's matrix as A/C/G/T rows by position columns, with its value type,
# or NULL when it cannot be read as one.
motif_matrix <- function(x) {
  type <- NULL
  if (is(object = x, class2 = "PWMatrix")) {
    type <- "logodds"
  } else if (is(object = x, class2 = "ICMatrix")) {
    type <- "other"
  }
  if (is(object = x, class2 = "XMatrix")) {
    x <- slot(object = x, name = "profileMatrix")
  }
  if (!is.matrix(x = x) || !is.numeric(x = x) || nrow(x = x) != 4L || ncol(x = x) == 0L ||
      ncol(x = x) > 65535L) {
    return(NULL)
  }
  bases <- c("A", "C", "G", "T")
  if (!is.null(x = rownames(x = x))) {
    if (!setequal(x = toupper(x = rownames(x = x)), y = bases)) {
      return(NULL)
    }
    x <- x[match(x = bases, table = toupper(x = rownames(x = x))), , drop = FALSE]
  }
  dimnames(x = x) <- list(bases, NULL)
  if (is.null(x = type)) {
    type <- if (anyNA(x = x)) {
      "other"
    } else if (any(x < 0)) {
      "logodds"
    } else if (all(abs(x = colSums(x = x) - 1) < 1e-3)) {
      "probabilities"
    } else if (all(x == round(x = x))) {
      "counts"
    } else {
      "other"
    }
  }
  list(type = type, values = x)
}

# Build provenance for the index header, as a named character vector of
# key/value strings in the order they are written. Arguments that do not
# apply to this build (NULL) are left out. `window` is not repeated here as
# the header already carries it; `seed = NULL` is recorded as
# "none" since an unseeded build is itself worth knowing about.
#' @importFrom utils packageVersion
index_build_info <- function(
    input,
    match_features,
    n_background,
    mindist,
    max_matches,
    p_cutoff,
    gc_jitter_sd,
    seed,
    exclude,
    bias,
    n_open_regions,
    n_mask_regions,
    genome = NULL
) {
  collapse <- function(x) {
    if (length(x = x) == 0L) "none" else paste(x, collapse = ",")
  }
  info <- list(
    build_date = format(x = Sys.time(), format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    bacha_version = as.character(x = packageVersion(pkg = "bacha")),
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    input = input,
    match_features = collapse(x = match_features),
    n_background = n_background,
    mindist = mindist,
    max_matches = max_matches,
    p_cutoff = p_cutoff,
    gc_jitter_sd = gc_jitter_sd,
    seed = seed %||% "none",
    exclude = collapse(x = exclude),
    bias = if (is.null(x = bias)) NULL else basename(path = bias),
    open_regions = n_open_regions,
    mask_regions = n_mask_regions,
    genome_version = genome_package_version(genome = genome)
  )
  info <- Filter(f = Negate(f = is.null), x = info)
  vapply(X = info, FUN = function(v) format(x = v, digits = 15L), FUN.VALUE = "")
}

# Version of the BSgenome package `genome` comes from, or NULL when it is not
# a BSgenome or its package is not installed.
genome_package_version <- function(genome) {
  if (!inherits(x = genome, what = "BSgenome")) {
    return(NULL)
  }
  tryCatch(
    expr = as.character(x = utils::packageVersion(pkg = genome@pkgname)),
    error = function(e) NULL
  )
}

# Append user-supplied `metadata` to the automatic build provenance. Checked
# up front so a bad entry fails before any index-building work runs.
add_build_metadata <- function(build_info, metadata) {
  if (is.null(x = metadata) || length(x = metadata) == 0L) {
    return(build_info)
  }
  keys <- names(x = metadata)
  if (is.null(x = keys) || anyNA(x = keys) || any(keys == "")) {
    stop("Every 'metadata' entry must be named", call. = FALSE)
  }
  if (!all(lengths(x = metadata) == 1L)) {
    stop("Each 'metadata' entry must be a single value", call. = FALSE)
  }
  clash <- intersect(x = keys, y = names(x = build_info))
  if (length(x = clash) > 0L) {
    stop(
      "'metadata' cannot reuse keys BachaIndex() records itself: ",
      paste(clash, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyDuplicated(x = keys) > 0L) {
    stop("'metadata' names must be unique", call. = FALSE)
  }
  values <- vapply(
    X = metadata,
    FUN = function(v) format(x = v, digits = 15L),
    FUN.VALUE = ""
  )
  if (any(grepl(pattern = "[\t\r\n]", x = c(keys, values)))) {
    stop("'metadata' names and values cannot contain tabs or newlines",
         call. = FALSE)
  }
  # `bacha view` starts each build line with its key, and its readers treat
  # a line starting with "#" as a header
  if (any(startsWith(x = keys, prefix = "#"))) {
    stop("'metadata' names cannot start with '#'", call. = FALSE)
  }
  out <- c(build_info, setNames(object = values, nm = keys))
  check_build_info(x = out)
  out
}

# Order the site windows so the two flanks of each site sit next to each other.
#
# Sorting by `site_id` alternates the two sides (up1, down1, up2, down2, ...)
# which removes the systematic upstream vs downstream bias.
interleave_flanks <- function(x) {
  x[order(x$site_id)]
}

# Find the candidate control regions overlapping the one at position-sorted
# index `index`, including itself.
#
# Candidates are sorted by (chromosome, start), so anything overlapping this
# region must start within `max_width` bp of its start; that bounds the walk in
# both directions, and the end coordinates decide which of those actually
# overlap. Indices are returned in the unsorted (original) numbering so they can
# index the per-motif vector of used controls directly.
overlapping_controls <- function(index, chrom, starts, ends, max_width, order) {
  n <- length(x = starts)
  lo <- index
  while (lo > 1L &&
         chrom[lo - 1L] == chrom[index] &&
         starts[index] - starts[lo - 1L] < max_width) {
    lo <- lo - 1L
  }
  hi <- index
  while (hi < n &&
         chrom[hi + 1L] == chrom[index] &&
         starts[hi + 1L] - starts[index] < max_width) {
    hi <- hi + 1L
  }
  window <- lo:hi
  hits <- window[
    ends[window] >= starts[index] & starts[window] <= ends[index]
  ]
  order[hits]
}

# Invert a per-feature list of selected control rows into a per-control list of
# feature IDs, dropping controls that no feature selected.
#
# Returned as a CompressedIntegerList
# `write_index()` only ever subsets, takes `lengths()` of, and `unlist()`s this
# object, all of which behave identically on an IntegerList.
#
#' @importFrom IRanges PartitioningByWidth relist
invert_selections <- function(selected_list, n_bg, id_offset) {
  n_feat <- length(x = selected_list)
  id_offset <- as.integer(x = id_offset)

  # Number of features each control row was selected for
  n_per_ctrl <- integer(n_bg)
  for (i in seq_len(length.out = n_feat)) {
    sel <- selected_list[[i]]
    if (length(x = sel) > 0L) {
      if (anyDuplicated(x = sel) > 0L) {
        stop(
          "Feature ", i, " selected the same control region more than once. ",
          "Control counts and the index layout would both be wrong.",
          call. = FALSE
        )
      }
      n_per_ctrl[sel] <- n_per_ctrl[sel] + 1L
    }
  }

  used <- which(x = n_per_ctrl > 0L)
  widths <- n_per_ctrl[used]
  total <- sum(as.double(x = widths))
  if (total > .Machine$integer.max) {
    stop(
      "Too many control associations to index (", format(x = total, big.mark = ","),
      "); reduce 'n_background' or 'max_matches'.",
      call. = FALSE
    )
  }
  total <- as.integer(x = total)

  # 0-based write cursor into `flat` for each retained control row, laid out in
  # ascending control order so that `flat` is already grouped by control
  cursor <- integer(n_bg)
  if (length(x = used) > 0L) {
    cursor[used] <- c(0L, cumsum(x = widths)[-length(x = widths)])
  }

  # Fill in ascending feature order, so each control's IDs come out ascending
  # (the index format requires this for the within-row delta encoding).
  flat <- integer(total)
  for (i in seq_len(length.out = n_feat)) {
    sel <- selected_list[[i]]
    if (length(x = sel) > 0L) {
      at <- cursor[sel] + 1L
      flat[at] <- i + id_offset
      cursor[sel] <- at
    }
  }

  list(
    motif_id = relist(
      flesh = flat, skeleton = PartitioningByWidth(x = widths)
    ),
    used = used
  )
}

# Open a writer connection that pipes binary bytes through the htslib
# `bgzip` command and writes BGZF (Blocked GZIP) output to `filename`.
# bgzip is shipped with htslib (and bundled with samtools/bcftools/tabix
# installs); on macOS it is available via `brew install htslib`, on
# Debian/Ubuntu via `apt install tabix`. A clear error is raised at this
# point if the command is not on the user's PATH, so the failure happens
# before any heavy index-building work runs.
open_bgzf_writer <- function(filename) {
  bgzip_path <- Sys.which(names = "bgzip")
  if (bgzip_path == "") {
    stop(
      "The `bgzip` command from htslib is required to write the Bacha ",
      "index file. Install htslib (e.g. `brew install htslib` on macOS, ",
      "`apt install tabix` on Debian/Ubuntu, or build from source at ",
      "https://www.htslib.org) and ensure `bgzip` is on your PATH.",
      call. = FALSE
    )
  }
  cmd <- paste(shQuote(string = bgzip_path), "-c >", shQuote(string = filename))
  pipe(description = cmd, open = "wb")
}

# Validate build provenance before anything is written. `bacha view` prints
# each entry as a `key<TAB>value` line, so neither part may contain a tab or
# newline, and each must fit a u16 length prefix.
check_build_info <- function(x) {
  if (length(x = x) == 0L) {
    return(invisible(x = NULL))
  }
  keys <- names(x = x)
  if (!is.character(x = x) || is.null(x = keys) || anyNA(x = x) ||
      any(!nzchar(x = keys)) || anyDuplicated(x = keys) > 0L) {
    stop("'build_info' must be a character vector with unique, non-empty names",
         call. = FALSE)
  }
  if (length(x = x) > 65535L) {
    stop("'build_info' holds at most 65535 entries", call. = FALSE)
  }
  fields <- c(keys, unname(x))
  if (any(grepl(pattern = "[\t\r\n]", x = fields))) {
    stop("'build_info' keys and values cannot contain tabs or newlines",
         call. = FALSE)
  }
  if (any(startsWith(x = keys, prefix = "#"))) {
    stop("'build_info' keys cannot start with '#'", call. = FALSE)
  }
  if (any(nchar(x = fields, type = "bytes") > 65535L)) {
    stop("'build_info' keys and values must be under 64 KB each", call. = FALSE)
  }
  invisible(x = NULL)
}

# Write index in custom binary format
# Format: BACH magic (4 bytes) + version (u8) + header + columnar data,
# BGZF (Blocked GZIP) compressed via the htslib `bgzip` command. Position
# arrays are delta-encoded within each chromosome and the flattened control
# TF indices are delta-encoded within each row to reduce post-compression
# size. bacha's reader auto-detects gzip / xz / BGZF / uncompressed indexes
# from magic bytes, so files written by older versions of this function
# (xz) continue to load correctly.
#
#' @importFrom BiocGenerics sort start end
#' @importFrom GenomeInfoDb seqnames seqlengths
#' @importFrom IRanges width
write_index <- function(
    motifs,
    controls,
    motif_id,
    filename,
    n_motif,
    n_controls,
    genome,
    motif_set_name,
    extra_chroms = character(0),
    build_info = character(0),
    feature_names = NULL,
    matrices = NULL
) {
  check_build_info(x = build_info)

  # extract genome name from BSgenome object
  genome_name <- if (inherits(x = genome, what = "BSgenome")) {
    genome@pkgname
  } else {
    as.character(x = genome)
  }

  # determine window size from data
  win <- unique(x = width(x = motifs))
  if (length(x = win) != 1) {
    stop("All motif regions must have the same width")
  }

  # Check window sizes are the same for control regions
  ctrl_win <- unique(x = width(x = controls))
  if (length(x = ctrl_win) != 1 || !identical(x = ctrl_win, y = win)) {
    stop(
      "All control regions must have the same width as motif regions (", win,
      "); got ", paste(head(x = sort(ctrl_win), n = 5L), collapse = ", "),
      call. = FALSE
    )
  }

  # TF names and 1-based indices
  tf_names <- names(x = n_motif)
  n_tfs <- length(x = tf_names)
  if (n_tfs == 0L) {
    stop("No features to write to the index", call. = FALSE)
  }
  motif_index <- seq_len(length.out = n_tfs)
  names(x = motif_index) <- tf_names

  # Display names, by feature ID; a feature without one is shown by its ID.
  display_names <- unname(feature_names[tf_names])
  if (length(x = display_names) != n_tfs) {
    display_names <- rep(x = NA_character_, times = n_tfs)
  }
  display_names[is.na(x = display_names)] <- tf_names[is.na(x = display_names)]
  # `bacha view` prints IDs and names as tab-separated text
  if (any(grepl(pattern = "[\t\r\n]", x = c(tf_names, display_names)))) {
    stop("Feature IDs and names cannot contain tabs or newlines", call. = FALSE)
  }
  # `bacha view` starts feature, site and control lines with the feature ID,
  # and its readers treat a line starting with "#" as a header
  if (any(startsWith(x = tf_names, prefix = "#"))) {
    stop("Feature IDs cannot start with '#'", call. = FALSE)
  }

  # get chromosome sizes from the genome
  all_lengths <- seqlengths(x = genome)
  all_lengths <- all_lengths[!is.na(x = all_lengths)]
  chroms_used <- unique(x = c(
    as.character(x = seqnames(x = motifs)),
    as.character(x = seqnames(x = controls))
  ))

  # Every chromosome carrying a region must have a length in the genome,
  # otherwise chrom_to_id[] below yields NA
  missing_len <- setdiff(x = chroms_used, y = names(x = all_lengths))
  if (length(x = missing_len) > 0L) {
    stop(
      "No sequence length available for chromosome(s): ",
      paste(head(x = missing_len, n = 5L), collapse = ", "),
      if (length(x = missing_len) > 5L) ", ..." else "",
      ". Ensure seqlengths(genome) covers every region, or drop these ",
      "sequences with 'exclude'.",
      call. = FALSE
    )
  }

  # The table also names chromosomes that carry no site or control: those
  # dropped by 'exclude' or those whose regions were all filtered out
  # this is for genome mismatch detection in scoring
  extra_chroms <- setdiff(
    x = intersect(x = names(x = all_lengths), y = extra_chroms),
    y = chroms_used
  )

  # Keep the genome's own seqlevel order across the whole table.
  chrom_names <- names(x = all_lengths)[
    names(x = all_lengths) %in% c(chroms_used, extra_chroms)
  ]
  chrom_lengths <- all_lengths[chrom_names]

  n_chroms <- length(x = chrom_names)
  chrom_to_id <- setNames(
    object = seq_along(along.with = chrom_names) - 1L,
    nm = chrom_names
  )

  # helper: write length-prefixed UTF-8 string (u16 length + bytes)
  write_string <- function(s, con) {
    raw_bytes <- charToRaw(x = s)
    writeBin(as.integer(x = length(x = raw_bytes)), con, size = 2L, endian = "little")
    if (length(x = raw_bytes) > 0) {
      writeBin(raw_bytes, con)
    }
  }

  con <- open_bgzf_writer(filename = filename)
  on.exit(expr = tryCatch(close(con), error = function(e) NULL), add = TRUE)

  # magic + version
  writeBin(charToRaw(x = "BACH"), con)
  writeBin(1L, con, size = 1L)

  # header
  write_string(genome_name, con)
  write_string(motif_set_name, con)
  writeBin(as.integer(x = win), con, size = 4L, endian = "little")
  writeBin(as.integer(x = n_tfs), con, size = 4L, endian = "little")
  writeBin(as.integer(x = n_chroms), con, size = 4L, endian = "little")

  # build provenance: entry count, then key/value string pairs
  writeBin(length(x = build_info), con, size = 2L, endian = "little")
  for (key in names(x = build_info)) {
    write_string(key, con)
    write_string(build_info[[key]], con)
  }

  # feature table: ID then display name for each feature
  for (i in seq_len(length.out = n_tfs)) {
    write_string(tf_names[i], con)
    write_string(display_names[i], con)
  }

  # instance counts: site_counts [u32; n_tfs] then control_counts [u32; n_tfs]
  writeBin(as.integer(x = n_motif), con, size = 4L, endian = "little")
  writeBin(as.integer(x = n_controls), con, size = 4L, endian = "little")

  # chromosome table
  for (i in seq_along(along.with = chrom_names)) {
    write_string(chrom_names[i], con)
    writeBin(as.integer(x = chrom_lengths[i]), con, size = 4L, endian = "little")
  }

  # matrix section: presence flag, then per feature a type code, a width and
  # the values position by position (A, C, G, T within each position)
  if (is.null(x = matrices)) {
    writeBin(0L, con, size = 1L)
  } else {
    writeBin(1L, con, size = 1L)
    for (tf in tf_names) {
      m <- matrices[[tf]]
      if (is.null(x = m)) {
        writeBin(0L, con, size = 1L)
        writeBin(0L, con, size = 2L, endian = "little")
        next
      }
      writeBin(match(x = m$type, table = MATRIX_TYPES) - 1L, con, size = 1L)
      writeBin(ncol(x = m$values), con, size = 2L, endian = "little")
      writeBin(as.numeric(x = m$values), con, size = 4L, endian = "little")
    }
  }

  # sort motifs by chromosome then start
  motifs <- sort(x = motifs)

  # site positions (columnar): n_rows, chrom_ids, starts (delta), tf_indices
  n_sites <- length(x = motifs)
  writeBin(as.integer(x = n_sites), con, size = 4L, endian = "little")

  site_chroms <- as.integer(x = chrom_to_id[as.character(x = seqnames(x = motifs))])
  site_starts <- as.integer(x = start(x = motifs))
  site_tf_idx <- as.integer(x = motif_index[motifs$motif])
  site_starts_delta <- delta_per_chrom(chroms = site_chroms, starts = site_starts)

  write_varint_column(con = con, x = site_chroms)
  write_varint_column(con = con, x = site_starts_delta)
  write_varint_column(con = con, x = site_tf_idx)

  # Sort controls by chromosome then start
  controls$.motif_id_idx <- seq_along(controls)
  controls <- sort(x = controls)
  ctrl_order <- controls$.motif_id_idx
  controls$.motif_id_idx <- NULL

  # control positions (columnar): n_rows, chrom_ids, starts (delta), n_per_row,
  # total_indices, tf_indices (within-row delta)
  n_ctrl <- length(x = controls)
  writeBin(as.integer(x = n_ctrl), con, size = 4L, endian = "little")

  ctrl_chroms <- as.integer(x = chrom_to_id[as.character(x = seqnames(x = controls))])
  ctrl_starts <- as.integer(x = start(x = controls))
  n_per_row <- as.integer(x = lengths(x = motif_id))[ctrl_order]
  total_indices <- sum(as.double(x = n_per_row))
  if (total_indices > .Machine$integer.max) {
    stop(
      "Too many control associations to index (",
      format(x = total_indices, big.mark = ","),
      "); reduce 'n_background' or 'max_matches'.",
      call. = FALSE
    )
  }
  ctrl_starts_delta <- delta_per_chrom(chroms = ctrl_chroms, starts = ctrl_starts)

  write_varint_column(con = con, x = ctrl_chroms)
  write_varint_column(con = con, x = ctrl_starts_delta)
  write_varint_column(con = con, x = n_per_row)
  writeBin(as.integer(x = total_indices), con, size = 4L, endian = "little")
  write_ctrl_indices(
    con = con,
    motif_id = motif_id,
    ctrl_order = ctrl_order,
    n_per_row = n_per_row
  )

  close(con)
  invisible(x = filename)
}

# Write the flattened, within-row delta encoded control index stream as one
# varint column, in batches of whole rows.
#
# `motif_id` is in unsorted control order; `ctrl_order` maps sorted row -> the
# original row, and `n_per_row` is already in sorted order. The column's byte
# length precedes it, so the batches are walked twice: once to size the
# column, once to write it.
write_ctrl_indices <- function(con, motif_id, ctrl_order, n_per_row,
                               max_batch = 2e7) {
  n_rows <- length(x = n_per_row)
  # cumulative entry count so batch edges can be found without scanning
  cum <- cumsum(as.double(x = n_per_row))
  for_each_batch <- function(f) {
    start_row <- 1L
    while (start_row <= n_rows) {
      base <- if (start_row == 1L) 0 else cum[start_row - 1L]
      end_row <- findInterval(x = base + max_batch, vec = cum)
      # a single row wider than the batch budget still has to go out on its own
      if (end_row < start_row) {
        end_row <- start_row
      }
      rows <- seq.int(from = start_row, to = end_row)
      vals <- as.integer(
        x = unlist(x = motif_id[ctrl_order[rows]], use.names = FALSE)
      )
      f(delta_within_row(x = vals, n_per_row = n_per_row[rows]))
      start_row <- end_row + 1L
    }
  }
  n_bytes <- 0
  for_each_batch(f = function(d) n_bytes <<- n_bytes + varint_nbytes(x = d))
  write_u64(con = con, x = n_bytes)
  for_each_batch(f = function(d) writeBin(encode_varint(x = d), con))
  invisible(NULL)
}

# Total number of bytes `x` takes as unsigned LEB128 varints. `x` must be
# non-negative and below 2^31 (R integers); varints of up to 5 bytes would
# hold any u32.
varint_nbytes <- function(x) {
  check_varint_values(x = x)
  # double arithmetic: a long column's byte count can exceed integer range
  as.double(x = length(x = x)) + sum(x >= 128L) + sum(x >= 16384L) + sum(x >= 2097152L) +
    sum(x >= 268435456L)
}

check_varint_values <- function(x) {
  if (anyNA(x = x) || any(x < 0L)) {
    stop(
      "Index columns must hold non-negative values; positions or control ",
      "feature indices are out of order",
      call. = FALSE
    )
  }
}

# Encode `x` as unsigned LEB128 varints: 7 bits per byte, least significant
# group first, high bit set on every byte except a value's last.
encode_varint <- function(x) {
  x <- as.integer(x = x)
  check_varint_values(x = x)
  if (length(x = x) == 0L) {
    return(raw(length = 0L))
  }
  nb <- 1L + (x >= 128L) + (x >= 16384L) + (x >= 2097152L) + (x >= 268435456L)
  ends <- cumsum(x = as.double(x = nb))
  out <- raw(length = ends[length(x = ends)])
  first <- ends - nb
  for (j in 0:4) {
    m <- which(x = nb > j)
    if (length(x = m) == 0L) {
      break
    }
    byte <- bitwAnd(a = bitwShiftR(a = x[m], n = 7L * j), b = 127L)
    more <- nb[m] > j + 1L
    byte[more] <- byte[more] + 128L
    out[first[m] + j + 1] <- as.raw(x = byte)
  }
  out
}

# Write a non-negative whole number below 2^53 as a u64 LE.
write_u64 <- function(con, x) {
  halves <- c(x %% 2^32, x %/% 2^32)
  # writeBin takes signed 32-bit integers; reinterpret each half
  halves[halves >= 2^31] <- halves[halves >= 2^31] - 2^32
  writeBin(as.integer(x = halves), con, size = 4L, endian = "little")
}

# Write `x` as a varint column: its byte length (u64 LE), then the values.
# Encoded in chunks to bound memory on long columns.
write_varint_column <- function(con, x, chunk = 1e7) {
  x <- as.integer(x = x)
  write_u64(con = con, x = varint_nbytes(x = x))
  n <- length(x = x)
  from <- 1
  while (from <= n) {
    to <- min(from + chunk - 1, n)
    writeBin(encode_varint(x = x[from:to]), con)
    from <- to + 1
  }
  invisible(NULL)
}

# Per-chromosome delta encoding for sorted positions. Within each chromosome
# run the first start stays absolute and subsequent starts become differences
# from the previous start. The reader recovers absolute starts via a
# cumulative sum that resets at each chromosome boundary.
delta_per_chrom <- function(chroms, starts) {
  n <- length(x = starts)
  if (n <= 1L) return(starts)
  same_chrom <- chroms[-1L] == chroms[-n]
  out <- starts
  out[-1L] <- starts[-1L] - as.integer(x = same_chrom) * starts[-n]
  out
}

# Within-row delta encoding for the flattened control-index stream. Indices
# within a row are sorted ascending; the first stays absolute and subsequent
# indices become differences from the previous index. The reader recovers
# absolute indices via a cumulative sum that resets at each row boundary.
delta_within_row <- function(x, n_per_row) {
  n <- length(x = x)
  if (n <= 1L) return(x)
  row_id <- rep.int(x = seq_along(along.with = n_per_row), times = n_per_row)
  same_row <- row_id[-1L] == row_id[-n]
  out <- x
  out[-1L] <- x[-1L] - as.integer(x = same_row) * x[-n]
  out
}

# Check matching-feature matrix
#' @importFrom stats var
check_match_features <- function(x, what, features) {
  if (nrow(x = x) == 0L) {
    stop("No ", what, " regions remain to match on", call. = FALSE)
  }
  bad <- features[apply(X = x, MARGIN = 2L, FUN = function(v) any(!is.finite(x = v)))]
  if (length(x = bad) > 0L) {
    stop(
      "Non-finite values in ", what, " match feature(s): ",
      paste(bad, collapse = ", "),
      ". This usually means a region fell outside the data used to compute ",
      "the feature (for example a bias file shorter than the chromosome).",
      call. = FALSE
    )
  }
  flat <- features[apply(X = x, MARGIN = 2L, FUN = var) == 0]
  if (length(x = flat) > 0L) {
    stop(
      "Match feature(s) with zero variance across ", what, " regions: ",
      paste(flat, collapse = ", "),
      ". Remove them from 'match_features'.",
      call. = FALSE
    )
  }
  invisible(x = NULL)
}

normalize_feature_names <- function(x, n, prefix = "motif") {
  if (is.null(x) || length(x = x) != n) {
    x <- rep(NA_character_, n)
  }
  x <- as.character(x)
  missing_name <- is.na(x) | !nzchar(x)
  x[missing_name] <- paste0(prefix, "_", which(missing_name))
  # Disambiguate duplicates so that PWM variants sharing a gene symbol
  # (e.g. HOCOMOCO cluster variants loaded with `Signac::ReadPWM`'s
  # `short_names = TRUE`) don't get merged into one entry. make.unique
  # appends ".1", ".2", ... to all but the first occurrence.
  dup_count <- sum(duplicated(x = x))
  if (dup_count > 0L) {
    warning(
      "Found ", dup_count,
      " duplicate feature name(s); appending '.1', '.2', ... ",
      "to keep them distinct. Affected names: ",
      paste(utils::head(unique(x[duplicated(x = x)]), 5L), collapse = ", "),
      if (dup_count > 5L) ", ..." else "",
      call. = FALSE
    )
    x <- make.unique(names = x)
  }
  x
}

# Compute a per-region sequence statistic without materialising every sequence
# at once. `getSeq()` over a genome-wide motif flank set (tens of millions of
# regions) builds a DNAStringSet of one byte per base -- several GB before any
# arithmetic -- and the per-base matrices in add_entropy() multiply that again.
# Each region's statistic depends only on its own sequence, so processing in
# blocks gives bit-identical results with a bounded working set.
#
# `fun(seqs, len)` receives one block's DNAStringSet and widths, and returns one
# numeric value per sequence.
#
#' @importFrom Biostrings getSeq
#' @importFrom IRanges width
seq_stat_chunked <- function(regions, genome, fun, max_bases = 2e8) {
  n <- length(x = regions)
  if (n == 0L) {
    return(numeric(0))
  }
  len <- width(x = regions)
  # Size blocks by total bases so a large `window` shrinks the block rather
  # than inflating the working set.
  per_block <- max(1L, min(n, as.integer(x = max_bases %/% max(max(len), 1L))))
  out <- numeric(length = n)
  block_start <- seq.int(from = 1L, to = n, by = per_block)
  for (s in block_start) {
    e <- min(s + per_block - 1L, n)
    idx <- seq.int(from = s, to = e)
    out[idx] <- fun(getSeq(x = genome, regions[idx]), len[idx])
  }
  out
}

#' @importFrom Biostrings letterFrequency
#' @importFrom S4Vectors mcols elementNROWS
add_gc <- function(regions, genome) {
  gc_of <- function(seqs, len) {
    as.vector(x = letterFrequency(x = seqs, "CG")) / len * 100
  }
  if (inherits(x = regions, what = "GRangesList")) {
    all_gr <- unlist(x = regions, use.names = FALSE)
    mcols(all_gr)$GC <- seq_stat_chunked(
      regions = all_gr, genome = genome, fun = gc_of
    )
    regions <- split(
      x = all_gr,
      f = rep(x = names(x = regions), elementNROWS(x = regions))
    )
    return(regions)
  } else if (inherits(x = regions, what = "GRanges")) {
    mcols(x = regions)$GC <- seq_stat_chunked(
      regions = regions, genome = genome, fun = gc_of
    )
    return(regions)
  } else {
    stop("Input must be a GRanges or GRangesList")
  }
}

#' @importFrom Biostrings vcountPattern
#' @importFrom S4Vectors mcols
add_cpg <- function(regions, genome) {
  mcols(x = regions)$CpG <- seq_stat_chunked(
    regions = regions,
    genome = genome,
    fun = function(seqs, len) {
      vcountPattern(pattern = "CG", subject = seqs) / (len - 1L) * 100
    }
  )
  regions
}

#' @importFrom Biostrings letterFrequency
#' @importFrom S4Vectors mcols
add_entropy <- function(regions, genome) {
  if (!inherits(x = regions, what = "GRanges")) {
    stop("Input must be a GRanges")
  }
  mcols(x = regions)$entropy <- seq_stat_chunked(
    regions = regions,
    genome = genome,
    fun = function(seqs, len) {
      counts <- letterFrequency(x = seqs, letters = c("A", "C", "G", "T"))
      totals <- rowSums(x = counts)
      totals[totals == 0L] <- 1L  # avoid div-by-zero; entropy stays 0
      p <- counts / totals
      log_p <- log2(x = p)
      log_p[!is.finite(x = log_p)] <- 0  # 0 * log(0) := 0
      -rowSums(x = p * log_p)
    }
  )
  regions
}

#' @importFrom GenomicRanges distanceToNearest resize
#' @importFrom S4Vectors mcols queryHits
add_tss_distance <- function(regions, tss_positions) {
  if (!inherits(x = regions, what = "GRanges")) {
    stop("'regions' must be a GRanges object")
  }
  if (!inherits(x = tss_positions, what = "GRanges")) {
    stop("'tss_positions' must be a GRanges object")
  }

  region_centers <- resize(x = regions, width = 1, fix = "center")
  tss_centers <- resize(x = tss_positions, width = 1, fix = "center")

  distances <- rep(Inf, length(region_centers))
  overlap_seqlevels <- intersect(seqlevelsInUse(x = region_centers), seqlevelsInUse(x = tss_centers))
  region_idx <- which(as.character(seqnames(region_centers)) %in% overlap_seqlevels)

  if (length(overlap_seqlevels) == 0) {
    stop("No overlapping seqlevels between 'regions' and 'tss_positions'")
  }

  if (length(region_idx) > 0) {
    region_subset <- keepSeqlevels(region_centers[region_idx], overlap_seqlevels, pruning.mode = "coarse")
    tss_subset <- keepSeqlevels(tss_centers, overlap_seqlevels, pruning.mode = "coarse")
    if (length(tss_subset) == 0) {
      stop("No TSS positions remain on seqlevels shared with 'regions'")
    }
    hits <- distanceToNearest(x = region_subset, subject = tss_subset, ignore.strand = TRUE)
    distances[region_idx[queryHits(hits)]] <- mcols(hits)$distance
  }

  finite_mask <- is.finite(distances)
  if (any(!finite_mask)) {
    if (!any(finite_mask)) {
      stop("Unable to compute finite TSS distances for any regions")
    }
    distances[!finite_mask] <- max(distances[finite_mask]) * 2
  }

  mcols(regions)$tssDist <- log1p(distances)
  regions
}
