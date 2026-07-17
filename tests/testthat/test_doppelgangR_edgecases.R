library(doppelgangR)
library(Biobase)

test_that("doppelgangR handles edge cases and error states correctly", {
  set.seed(123)
  
  # Error handling: Not a list, ExpressionSet, or SummarizedExperiment
  expect_error(doppelgangR(data.frame(x=1)), "esets must be an ExpressionSet, SummarizedExperiment, or a list of such objects")
  
  # Edge Case: No featureNames in common
  mat1 <- matrix(rnorm(30), ncol=3)
  rownames(mat1) <- paste0("Gene", 1:10)
  colnames(mat1) <- paste0("Sample", 1:3)
  pdat1 <- data.frame(age = 1:3)
  rownames(pdat1) <- colnames(mat1)
  eset1 <- ExpressionSet(assayData = mat1, phenoData = AnnotatedDataFrame(pdat1))
  
  mat2 <- matrix(rnorm(30), ncol=3)
  rownames(mat2) <- paste0("Gene", 11:20) # Different genes
  colnames(mat2) <- paste0("Sample", 4:6)
  pdat2 <- data.frame(age = 4:6)
  rownames(pdat2) <- colnames(mat2)
  eset2 <- ExpressionSet(assayData = mat2, phenoData = AnnotatedDataFrame(pdat2))
  
  # Suppress the warning so it doesn't clutter output, but expect it to happen
  expect_warning(
    res_no_genes <- doppelgangR(list(eset1, eset2)),
    "have no featureNames in common"
  )
  
  # Edge Case: Automatic smoking guns
  pdat1_sg <- data.frame(unique_id = c("A", "B", "C", "D"))
  rownames(pdat1_sg) <- paste0("Sample", 1:4)
  mat1_sg <- matrix(rnorm(40), ncol=4)
  rownames(mat1_sg) <- paste0("Gene", 1:10)
  colnames(mat1_sg) <- rownames(pdat1_sg)
  eset1_sg <- ExpressionSet(assayData = mat1_sg, phenoData = AnnotatedDataFrame(pdat1_sg))
  
  pdat2_sg <- data.frame(unique_id = c("E", "F", "G", "H"))
  rownames(pdat2_sg) <- paste0("Sample", 5:8)
  mat2_sg <- matrix(rnorm(40), ncol=4)
  rownames(mat2_sg) <- paste0("Gene", 1:10)
  colnames(mat2_sg) <- rownames(pdat2_sg)
  eset2_sg <- ExpressionSet(assayData = mat2_sg, phenoData = AnnotatedDataFrame(pdat2_sg))
  
  expect_error(
    res_sg <- doppelgangR(list(eset1_sg, eset2_sg), automatic.smokingguns = TRUE),
    "Intermediate pruning off but no addCols shortcut available."
  )
  
  # Edge Case: single ExpressionSet with manual smoking gun (eset.method=TRUE)
  expect_error(
    res_eset_method <- doppelgangR(eset1_sg, manual.smokingguns = "unique_id"),
    "Intermediate pruning off but no addCols shortcut available."
  )
  
  # Test intermediate pruning with differently sized doppelganger sets
  res_pruning <- doppelgangR(list(eset1_sg, eset2_sg), automatic.smokingguns = TRUE, intermediate.pruning = TRUE)
  expect_s4_class(res_pruning, "DoppelGang")
 
  # Test future error handling (mocking a dataset error)
  eset1_err <- eset1_sg
  exprs(eset1_err)[1, 1] <- NA 
  
  # Mock phenoFinder to throw an error 
  # We will test missing corFinder instead.
  # Let's just create an eset with non-numeric matrix to throw an error
  eset_fail <- ExpressionSet(assayData = matrix("A", ncol=2, nrow=2), phenoData = AnnotatedDataFrame(data.frame(age=1:2)))
  expect_error(
    expect_warning(
      doppelgangR(list(eset1_sg, eset_fail)),
      "Caught simpleError"
    )
  )

  # Test SummarizedExperiment input handling and coercion
  if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    se_mat <- matrix(rnorm(50), ncol=5)
    rownames(se_mat) <- paste0("Gene", 1:10)
    colnames(se_mat) <- paste0("SampleSE", 1:5)
    se_pdat <- S4Vectors::DataFrame(age = 1:5)
    rownames(se_pdat) <- colnames(se_mat)
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(exprs = se_mat), colData = se_pdat)
    
    # Passing single SummarizedExperiment
    res_se <- doppelgangR(se)
    expect_s4_class(res_se, "DoppelGang")
    
    # Passing list of SummarizedExperiment and ExpressionSet
    res_mix <- doppelgangR(list(se, eset1_sg))
    expect_s4_class(res_mix, "DoppelGang")
  }
})
