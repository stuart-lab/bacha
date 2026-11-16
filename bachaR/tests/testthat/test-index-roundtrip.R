# Writes a small index and reads it back through the bacha CLI. Set
# BACHA_PATH to test a binary that is not first on PATH.
bacha_bin <- function() {
  path <- Sys.getenv("BACHA_PATH", unset = unname(Sys.which("bacha")))
  if (!nzchar(path)) skip("bacha executable not found")
  if (!nzchar(Sys.which("bgzip"))) skip("bgzip not found")
  path
}

write_toy_index <- function(path) {
  si <- GenomeInfoDb::Seqinfo(
    seqnames = c("chr1", "chr2", "chr3"), seqlengths = c(1e5, 5e4, 2e4)
  )
  gr <- function(chr, starts) {
    GenomicRanges::GRanges(chr, IRanges::IRanges(starts, width = 10), seqinfo = si)
  }
  motifs <- gr(c("chr1", "chr1", "chr2"), c(100, 300, 500))
  motifs$motif <- c("MA1", "MA2", "MA1")
  pfm <- matrix(c(5, 0, 0, 5, 0, 5, 5, 0), nrow = 4,
                dimnames = list(c("A", "C", "G", "T"), NULL))
  write_index(
    motifs = motifs,
    controls = gr("chr1", c(1000, 2000)),
    motif_id = list(c(1L, 2L), 1L),
    filename = path,
    n_motif = c(MA1 = 2L, MA2 = 1L),
    n_controls = c(MA1 = 2L, MA2 = 1L),
    genome = gr("chr1", 1),
    motif_set_name = "toy",
    extra_chroms = "chr3",
    build_info = c(seed = "1"),
    feature_names = c(MA1 = "Alpha", MA2 = "Beta"),
    matrices = list(MA1 = list(type = "counts", values = pfm))
  )
  pfm
}

test_that("an index written by R reads back through the CLI", {
  bin <- bacha_bin()
  path <- tempfile(fileext = ".bi")
  on.exit(unlink(path))
  pfm <- write_toy_index(path)

  view <- system2(bin, c("view", path), stdout = TRUE)
  section <- function(name) {
    i <- match(name, view) + 1L
    if (startsWith(view[i], "# ")) i <- i + 1L
    j <- i
    while (j <= length(view) && !startsWith(view[j], "#")) j <- j + 1L
    view[seq_len(j - i) + i - 1L]
  }
  expect_true("## motif set: toy" %in% view)
  expect_equal(section("##[build]"), "seed\t1")
  expect_equal(section("##[features]"), c("MA1\tAlpha", "MA2\tBeta"))
  # chr3 has no sites but is kept for assembly-mismatch detection
  expect_equal(
    section("##[chromosome sizes]"),
    c("chr1\t100000", "chr2\t50000", "chr3\t20000")
  )
  expect_equal(
    section("##[site positions]"),
    c("chr1\t100\t109\t1", "chr1\t300\t309\t2", "chr2\t500\t509\t1")
  )
  expect_equal(
    section("##[control positions]"),
    c("chr1\t1000\t1009\t1,2", "chr1\t2000\t2009\t1")
  )

  capture.output(
    tab <- IndexSummary(path, bacha.path = bin, return.table = TRUE)
  )
  expect_equal(tab$name[match(c("MA1", "MA2"), tab$tf)], c("Alpha", "Beta"))
  expect_equal(tab$n_sites[match(c("MA1", "MA2"), tab$tf)], c(2L, 1L))

  # only MA1 was given a matrix
  motifs <- IndexMotifs(path, bacha.path = bin)
  expect_named(motifs, "MA1")
  expect_equal(attr(motifs$MA1, "type"), "counts")
  expect_equal(as.vector(motifs$MA1), as.vector(pfm))
})
