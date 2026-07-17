library(Biobase)
library(doppelgangR)

test_that("smokingGunFinder handles input correctly (Error Handling)", {
  mat <- matrix(runif(20), 4, 5)
  colnames(mat) <- paste0("Sample", 1:5)
  pdat <- data.frame(id = 1:5, var2 = letters[1:5])
  rownames(pdat) <- colnames(mat)
  eset <- ExpressionSet(assayData = mat, phenoData = AnnotatedDataFrame(pdat))
  
  expect_error(smokingGunFinder(list(eset), "id"), "eset.pair should be a list of two ExpressionSets")
  expect_error(smokingGunFinder(eset, "id"), "eset.pair should be a list of two ExpressionSets")
})

test_that("smokingGunFinder works with identical ExpressionSets (Normal Use)", {
  mat <- matrix(runif(20), 4, 5)
  colnames(mat) <- paste0("Sample", 1:5)
  pdat <- data.frame(id = c(1, 2, 3, 2, 4), var2 = letters[1:5])
  rownames(pdat) <- colnames(mat)
  eset <- ExpressionSet(assayData = mat, phenoData = AnnotatedDataFrame(pdat))
  
  res <- smokingGunFinder(list(eset, eset), "id")
  expect_true(is.matrix(res))
  expect_equal(dim(res), c(5, 5))
  # Should be upper triangular with NAs in lower
  expect_true(all(is.na(res[!upper.tri(res)])))
  # Check if it found the duplicate id (2) at indices 2 and 4
  expect_equal(res[2, 4], 1)
})

test_that("smokingGunFinder works with different ExpressionSets (Normal Use)", {
  mat1 <- matrix(runif(20), 4, 5)
  colnames(mat1) <- paste0("Sample", 1:5)
  pdat1 <- data.frame(id = 1:5, var2 = letters[1:5])
  rownames(pdat1) <- colnames(mat1)
  eset1 <- ExpressionSet(assayData = mat1, phenoData = AnnotatedDataFrame(pdat1))
  
  mat2 <- matrix(runif(12), 4, 3)
  colnames(mat2) <- paste0("OtherSample", 1:3)
  pdat2 <- data.frame(id = c(1, 6, 3), var2 = c("A", "C", "D"))
  rownames(pdat2) <- colnames(mat2)
  eset2 <- ExpressionSet(assayData = mat2, phenoData = AnnotatedDataFrame(pdat2))
  
  res <- smokingGunFinder(list(eset1, eset2), "id")
  expect_true(is.matrix(res))
  expect_equal(dim(res), c(5, 3))
  # Check if it found matches: 1 and 3 are present in both
  expect_equal(res[1, 1], 1) # Sample1 vs OtherSample1 (id 1)
  expect_equal(res[3, 3], 1) # Sample3 vs OtherSample3 (id 3)
  expect_equal(sum(res), 2)
})

test_that("smokingGunFinder handles transformation function (Edge Cases)", {
  mat1 <- matrix(runif(20), 4, 5)
  colnames(mat1) <- paste0("Sample", 1:5)
  pdat1 <- data.frame(id = c("A", "B", "C", "D", "E"))
  rownames(pdat1) <- colnames(mat1)
  eset1 <- ExpressionSet(assayData = mat1, phenoData = AnnotatedDataFrame(pdat1))
  
  mat2 <- matrix(runif(12), 4, 3)
  colnames(mat2) <- paste0("OtherSample", 1:3)
  pdat2 <- data.frame(id = c("a", "f", "c"))
  rownames(pdat2) <- colnames(mat2)
  eset2 <- ExpressionSet(assayData = mat2, phenoData = AnnotatedDataFrame(pdat2))
  
  res_no_trans <- smokingGunFinder(list(eset1, eset2), "id")
  expect_equal(sum(res_no_trans), 0) # No matches because case differs
  
  res_trans <- smokingGunFinder(list(eset1, eset2), "id", transFun = toupper)
  expect_equal(sum(res_trans), 2) # Found A and C after toupper
})

test_that("smokingGunFinder handles missing columns and NAs", {
  mat1 <- matrix(runif(20), 4, 5)
  colnames(mat1) <- paste0("Sample", 1:5)
  pdat1 <- data.frame(id = c(1, 2, NA, 4, 5)) # Contains NA
  rownames(pdat1) <- colnames(mat1)
  eset1 <- ExpressionSet(assayData = mat1, phenoData = AnnotatedDataFrame(pdat1))
  
  mat2 <- matrix(runif(12), 4, 3)
  colnames(mat2) <- paste0("OtherSample", 1:3)
  pdat2 <- data.frame(id = c(1, NA, 3), other_col = c("A", "B", "C")) # Contains NA and missing a col
  rownames(pdat2) <- colnames(mat2)
  eset2 <- ExpressionSet(assayData = mat2, phenoData = AnnotatedDataFrame(pdat2))
  
  # Test missing column (skips and returns 0 matrix if only that column is provided)
  res_missing <- smokingGunFinder(list(eset1, eset2), "other_col")
  expect_equal(sum(res_missing), 0)
  
  # Test NA handling: NAs shouldn't match with anything
  res_na <- smokingGunFinder(list(eset1, eset2), "id")
  expect_equal(sum(res_na), 1) # Only id 1 matches. NAs do not match.
})
