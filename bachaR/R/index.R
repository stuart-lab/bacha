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
#' - 2 bytes (u16 LE): number of TFs (`n_tfs`)
#' - 1 byte (u8): number of chromosomes (`n_chroms`)
#'
#' **TF names**
#'
#' - `n_tfs` length-prefixed strings: TF names
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
#' **Site positions** (sorted by chromosome then start)
#'
#' - 4 bytes (u32 LE): number of site rows (`n_sites`)
#' - `[u8; n_sites]`: chromosome IDs (0-based index into chromosome table)
#' - `[u32 LE; n_sites]`: start positions (1-based), delta-encoded within each
#'   chromosome run. The first start of each chromosome is absolute;
#'   subsequent starts within the same chromosome are differences from the
#'   previous start. Reader recovers absolute starts via cumulative sum
#'   reset at each chromosome boundary, then converts to 0-based half-open
#'   intervals for fragment endpoint counting.
#' - `[u16 LE; n_sites]`: TF indices (1-based index into TF names)
#'
#' **Control positions** (sorted by chromosome then start)
#'
#' - 4 bytes (u32 LE): number of control rows (`n_ctrl`)
#' - `[u8; n_ctrl]`: chromosome IDs
#' - `[u32 LE; n_ctrl]`: 1-based start positions, delta-encoded within each
#'   chromosome run (same scheme as site starts).
#' - `[u16 LE; n_ctrl]`: number of TF associations per control (`n_per_row`)
#' - 4 bytes (u32 LE): total number of TF index entries (`total_indices`)
#' - `[u16 LE; total_indices]`: flattened TF indices for each control row,
#'   delta-encoded within each row. Within a row the indices are sorted
#'   ascending; the first index of each row is absolute and subsequent
#'   indices are differences from the previous index. Each absolute index
#'   is encoded as `tf_index + n_tfs` to distinguish controls from sites
#'   in the shared ID space.
#'
#' @param pwm A list of position weight matrices or position frequency matrices.
#'   If provided, the `positions` argument should not be set.
#' @param genome A BSgenome object.
#' @param bias Optional path to a H5 file containing Tn5 insertion bias scores
#'   for each base in the genome. Required only when `match_features` includes
#'   `"meanBias"` or `"maxBias"`. See <https://zenodo.org/records/15224770>;
#'   <https://doi.org/10.1038/s41586-024-08443-4>.
#' @param open_regions A GRanges object containing open chromatin regions.
#' @param mask_regions A GRanges object containing regions to mask from the
#'   motif scan (e.g. blacklist regions or repeat regions).
#' @param output_path Path to write the output index file.
#' @param positions A [GenomicRanges::GRangesList] object containing
#'   pre-computed genomic positions of interest for a set of features. Each
#'   element should be a different feature (e.g. a TF motif or GWAS trait) and
#'   the genomic regions correspond to regions relevant to that feature to
#'   aggregate information across. If `positions` are supplied, `pwm` should
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
#' @importFrom BiocGenerics start<- end
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
    names(x = pwm) <- normalize_feature_names(x = names(x = pwm), n = length(x = pwm))
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
    motif.names <- rep(
      x = names(x = positions), elementNROWS(x = positions)
    )
    positions <- unlist(x = positions, use.names = FALSE)
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
    motif_set_name = motif_set_name %||% "custom"
  )
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
    extra_chroms = character(0)
) {
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
  # Feature indices are written as u16, and controls are stored as
  # `tf_index + n_tfs`, so the largest value written is 2 * n_tfs.
  if (n_tfs == 0L) {
    stop("No features to write to the index", call. = FALSE)
  }
  if (2L * n_tfs > 65535L) {
    stop("Index format supports at most ", 65535L %/% 2L, " features, got ",
         n_tfs, call. = FALSE)
  }
  motif_index <- seq_len(length.out = n_tfs)
  names(x = motif_index) <- tf_names

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
  n_room <- 255L - length(x = chroms_used)
  if (length(x = extra_chroms) > n_room) {
    warning(
      "The chromosome table holds at most 255 entries; omitting ",
      length(x = extra_chroms) - max(n_room, 0L),
      " region-free chromosome(s). Fragments on those chromosomes will be ",
      "reported as absent from the index.",
      call. = FALSE
    )
    extra_chroms <- if (n_room > 0L) extra_chroms[seq_len(n_room)] else character(0)
  }

  # Keep the genome's own seqlevel order across the whole table.
  chrom_names <- names(x = all_lengths)[
    names(x = all_lengths) %in% c(chroms_used, extra_chroms)
  ]
  chrom_lengths <- all_lengths[chrom_names]

  n_chroms <- length(x = chrom_names)
  if (n_chroms > 255L) {
    stop("Index format supports at most 255 chromosomes, got ", n_chroms,
         ". Consider filtering scaffolds with the 'exclude' parameter.")
  }
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
  writeBin(as.integer(x = n_tfs), con, size = 2L, endian = "little")
  writeBin(as.integer(x = n_chroms), con, size = 1L)

  # TF names
  for (tf in tf_names) write_string(tf, con)

  # instance counts: site_counts [u32; n_tfs] then control_counts [u32; n_tfs]
  writeBin(as.integer(x = n_motif), con, size = 4L, endian = "little")
  writeBin(as.integer(x = n_controls), con, size = 4L, endian = "little")

  # chromosome table
  for (i in seq_along(along.with = chrom_names)) {
    write_string(chrom_names[i], con)
    writeBin(as.integer(x = chrom_lengths[i]), con, size = 4L, endian = "little")
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

  writeBin(site_chroms, con, size = 1L)
  writeBin(site_starts_delta, con, size = 4L, endian = "little")
  writeBin(site_tf_idx, con, size = 2L, endian = "little")

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

  writeBin(ctrl_chroms, con, size = 1L)
  writeBin(ctrl_starts_delta, con, size = 4L, endian = "little")
  writeBin(n_per_row, con, size = 2L, endian = "little")
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

# Write the flattened, within-row delta encoded control index stream in batches
# of whole rows.
#
# `motif_id` is in unsorted control order; `ctrl_order` maps sorted row -> the
# original row, and `n_per_row` is already in sorted order.
write_ctrl_indices <- function(con, motif_id, ctrl_order, n_per_row,
                               max_batch = 2e7) {
  n_rows <- length(x = n_per_row)
  if (n_rows == 0L) {
    return(invisible(NULL))
  }
  # cumulative entry count so batch edges can be found without scanning
  cum <- cumsum(as.double(x = n_per_row))
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
    writeBin(
      delta_within_row(x = vals, n_per_row = n_per_row[rows]),
      con,
      size = 2L,
      endian = "little"
    )
    start_row <- end_row + 1L
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
