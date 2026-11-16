test_that("metadata is added to the build block", {
  out <- add_build_metadata(c(seed = "1"), list(source = "JASPAR", n = 3))
  expect_equal(out, c(seed = "1", source = "JASPAR", n = "3"))
  expect_equal(add_build_metadata(c(seed = "1"), NULL), c(seed = "1"))
})

test_that("invalid metadata is rejected", {
  base <- c(seed = "1")
  expect_error(add_build_metadata(base, c("x")), "must be named")
  expect_error(add_build_metadata(base, list(a = 1:2)), "single value")
  expect_error(add_build_metadata(base, c(seed = "2")), "cannot reuse")
  expect_error(add_build_metadata(base, c(a = "1\t2")), "tabs or newlines")
  expect_error(add_build_metadata(base, c(`#a` = "1")), "start with '#'")
})
