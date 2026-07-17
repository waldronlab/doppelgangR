library(doppelgangR)
library(Biobase)

test_that("corFinder logic handles inputs, errors, and missing/infinite values", {
  set.seed(123)
  
  # Create a dummy expression matrix
  mat1 <- matrix(rnorm(30), ncol=3)
  rownames(mat1) <- paste0("Gene", 1:10)
  colnames(mat1) <- paste0("Sample", 1:3)
  
  mat2 <- matrix(rnorm(30), ncol=3)
  rownames(mat2) <- paste0("Gene", 1:10)
  colnames(mat2) <- paste0("Sample", 4:6)
  
  eset1 <- ExpressionSet(assayData = mat1)
  eset2 <- ExpressionSet(assayData = mat2)
  
  # Error Handling: Not a list or not length 2
  expect_error(corFinder(eset1), "eset.pair should be a list of two ExpressionSets")
  
  # Edge Case: Introduce NA/Inf to test non-finite filtering
  mat1_inf <- mat1
  mat1_inf[1, 1] <- Inf
  mat1_inf[2, 2] <- NA
  eset1_inf <- ExpressionSet(assayData = mat1_inf)
  
  # When one row has Inf and one has NA, only 8 rows are left, which should be fine
  res_inf <- corFinder(list(A=eset1_inf, B=eset2), use.ComBat = TRUE)
  expect_true(is.matrix(res_inf))
  
  # Error Handling: Not enough finite rows left
  mat1_all_inf <- mat1
  mat1_all_inf[1:9, 1] <- Inf # 9 out of 10 rows have Inf
  eset1_all_inf <- ExpressionSet(assayData = mat1_all_inf)
  expect_error(corFinder(list(A=eset1_all_inf, B=eset2), use.ComBat = TRUE), "Fewer than two genes without all finite values")
})
