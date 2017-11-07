context("Testing normalization methods output")

test_that("libsize_factors are not null", {
  expect_true(is.numeric(libsize_factors))
})

test_that("cpm output", {
  expect_equal(dim(cpm), dim(counts))
})

