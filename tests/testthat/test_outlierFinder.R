library(doppelgangR)

test_that("outlierFinder handles Edge Cases and Error Handling", {
  set.seed(1)
  
  # Create a dummy similarity matrix
  mat <- matrix(runif(100, min = 0, max = 1), nrow = 10, ncol = 10)
  diag(mat) <- 1
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  rownames(mat) <- colnames(mat) <- paste0("Sample", 1:10)
  
  # Edge Case: transFun = NULL
  res_null_trans <- outlierFinder(mat, transFun = NULL)
  expect_true(is.list(res_null_trans) && "outlierFinder.res" %in% names(res_null_trans))
  
  # Edge Case: tail = "lower"
  res_lower <- outlierFinder(mat, tail = "lower")
  expect_true(is.list(res_lower) && "outlierFinder.res" %in% names(res_lower))
  
  # Edge Case: tail = "both"
  res_both <- outlierFinder(mat, tail = "both")
  expect_true(is.list(res_both) && "outlierFinder.res" %in% names(res_both))
  
  # Error Handling: tail = "invalid"
  expect_error(outlierFinder(mat, tail = "invalid"), "tail argument should be upper, lower, or both")
  
  # Error Handling: bonf.prob = NULL and normal.upper.thresh = NULL
  res_both_null <- outlierFinder(mat, bonf.prob = NULL, normal.upper.thresh = NULL)
  expect_null(res_both_null)
})
