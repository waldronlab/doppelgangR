library(SummarizedExperiment)
library(doppelgangR)

test_that("colFinder handles input correctly (Error Handling)", {
  se <- SummarizedExperiment(matrix(runif(20), 4, 5))
  colData(se) <- DataFrame(var1 = 1:5, var2 = letters[1:5])
  
  expect_error(colFinder(list(se)), "list should be a list of length 2")
  expect_error(colFinder(se), "list should be a list of length 2")
  
  se2 <- se
  colData(se2) <- DataFrame(var3 = 1:5, var4 = letters[1:5])
  expect_error(colFinder(list(se, se2)), "Slots of list must have identical column names")
})

test_that("colFinder works with identical SummarizedExperiments (Normal Use)", {
  se1 <- SummarizedExperiment(matrix(runif(20), 4, 5))
  colData(se1) <- DataFrame(var1 = 1:5, var2 = c("A", "A", "B", "B", "C"))
  colnames(se1) <- paste0("Sample", 1:5)
  
  res <- colFinder(list(se1, se1))
  expect_true(is.matrix(res))
  expect_equal(dim(res), c(5, 5))
  # When identical, it should return upper triangular matrix with NAs in lower
  expect_true(all(is.na(res[!upper.tri(res)])))
})

test_that("colFinder works with different SummarizedExperiments (Normal Use)", {
  se1 <- SummarizedExperiment(matrix(runif(20), 4, 5))
  colData(se1) <- DataFrame(var1 = 1:5, var2 = c("A", "A", "B", "B", "C"))
  colnames(se1) <- paste0("Sample", 1:5)
  
  se2 <- SummarizedExperiment(matrix(runif(12), 4, 3))
  colData(se2) <- DataFrame(var1 = c(1, 2, 6), var2 = c("A", "A", "C"))
  colnames(se2) <- paste0("OtherSample", 1:3)
  
  res <- colFinder(list(se1, se2))
  expect_true(is.matrix(res))
  expect_equal(dim(res), c(5, 3))
})

test_that("colFinder handles missing rownames in colData", {
  se1 <- SummarizedExperiment(matrix(runif(20), 4, 5))
  colData(se1) <- DataFrame(var1 = 1:5, var2 = c("A", "A", "B", "B", "C"), row.names = NULL)
  colnames(se1) <- NULL
  
  se2 <- SummarizedExperiment(matrix(runif(12), 4, 3))
  colData(se2) <- DataFrame(var1 = c(1, 2, 6), var2 = c("A", "A", "C"), row.names = NULL)
  colnames(se2) <- NULL
  
  res <- colFinder(list(se1, se2))
  expect_true(is.matrix(res))
  expect_equal(dim(res), c(5, 3))
  expect_equal(rownames(res), make.names(1:5))
  expect_equal(colnames(res), make.names(1:3))
})
