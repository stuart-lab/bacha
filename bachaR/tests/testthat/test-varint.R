test_that("encode_varint writes unsigned LEB128", {
  expect_equal(encode_varint(0L), as.raw(0x00))
  expect_equal(encode_varint(127L), as.raw(0x7f))
  expect_equal(encode_varint(128L), as.raw(c(0x80, 0x01)))
  expect_equal(encode_varint(300L), as.raw(c(0xac, 0x02)))
  expect_equal(
    encode_varint(.Machine$integer.max),
    as.raw(c(0xff, 0xff, 0xff, 0xff, 0x07))
  )
  expect_length(encode_varint(integer(0)), 0L)
})

test_that("varint_nbytes matches the encoded length", {
  x <- c(0L, 127L, 128L, 16383L, 16384L, 2097152L, 268435456L)
  expect_equal(varint_nbytes(x), length(encode_varint(x)))
})

test_that("negative values are rejected", {
  expect_error(encode_varint(c(1L, -1L)), "non-negative")
})

test_that("starts are delta-encoded within each chromosome", {
  chroms <- c(1L, 1L, 1L, 2L, 2L)
  starts <- c(10L, 15L, 40L, 7L, 9L)
  expect_equal(delta_per_chrom(chroms, starts), c(10L, 5L, 25L, 7L, 2L))
})

test_that("control indices are delta-encoded within each row", {
  expect_equal(
    delta_within_row(c(2L, 5L, 9L, 1L, 4L), n_per_row = c(3L, 2L)),
    c(2L, 3L, 4L, 1L, 3L)
  )
})
