library(doppelgangR)
library(Biobase)

test_that("phenoFinder logic handles inputs, errors, and missing column names", {
  set.seed(123)
  
  # Create a dummy dataframe with some clinical attributes
  pdata1 <- data.frame(
    age = c(50, 60, 55),
    stage = c("I", "II", "III"),
    stringsAsFactors = FALSE
  )
  rownames(pdata1) <- c("Sample1", "Sample2", "Sample3")
  
  pdata2 <- data.frame(
    age = c(51, 62, 53),
    stage = c("I", "II", "IV"),
    stringsAsFactors = FALSE
  )
  rownames(pdata2) <- c("Sample4", "Sample5", "Sample6")
  
  # Create ExpressionSets
  mat1 <- matrix(rnorm(30), ncol=3)
  colnames(mat1) <- rownames(pdata1)
  eset1 <- ExpressionSet(assayData = mat1, phenoData = AnnotatedDataFrame(pdata1))
  
  mat2 <- matrix(rnorm(30), ncol=3)
  colnames(mat2) <- rownames(pdata2)
  eset2 <- ExpressionSet(assayData = mat2, phenoData = AnnotatedDataFrame(pdata2))
  
  # Error Handling: Not a list or not length 2
  expect_error(phenoFinder(eset1), "eset.pair should be a list of length 2")
  expect_error(phenoFinder(list(eset1, eset2, eset1)), "eset.pair should be a list of length 2")
  
  # Error Handling: Different column names
  pdata_diff <- pdata2
  colnames(pdata_diff)[1] <- "patient_age"
  mat_diff <- matrix(rnorm(30), ncol=3)
  colnames(mat_diff) <- rownames(pdata_diff)
  eset_diff <- ExpressionSet(assayData = mat_diff, phenoData = AnnotatedDataFrame(pdata_diff))
  expect_error(phenoFinder(list(eset1, eset_diff)), "pData slots of esets must have identical column names")
  
  # Edge Case: Missing rownames get filled with make.names(1:nrow)
  rownames(pData(eset1)) <- NULL
  rownames(pData(eset2)) <- NULL
  res_missing_rows <- phenoFinder(list(eset1, eset2))
  expect_true(is.matrix(res_missing_rows))
  expect_equal(rownames(res_missing_rows), c("X1", "X2", "X3"))
})
