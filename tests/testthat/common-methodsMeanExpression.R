context("Testing fit output")

test_that("pvalues are not all NAs", {
  expect_true(sum(is.na(as.character(output$pvalue))) != length(output$pvalue))
})

