library(Biobase)
library(doppelgangR)

# Setup dummy data exactly as in test_main.R to ensure doppelgangR succeeds
set.seed(1)
m1 <- matrix(rnorm(1100), ncol = 11)
colnames(m1) <- paste("m", 1:11, sep = "")
rownames(m1) <- make.names(1:nrow(m1))
n1 <- matrix(rnorm(1000), ncol = 10)
colnames(n1) <- paste("n", 1:10, sep = "")
rownames(n1) <- make.names(1:nrow(n1))
m.pdata <- matrix(letters[sample(1:26, size = 110, replace = TRUE)], ncol = 10)
rownames(m.pdata) <- colnames(m1)
n.pdata <- matrix(letters[sample(1:26, size = 100, replace = TRUE)], ncol = 10)
rownames(n.pdata) <- colnames(n1)
m.eset <- ExpressionSet(assayData = m1, phenoData = AnnotatedDataFrame(data.frame(m.pdata)))
m.eset$id <- toupper(colnames(m1))
n.eset <- ExpressionSet(assayData = n1, phenoData = AnnotatedDataFrame(data.frame(n.pdata)))
n.eset$id <- toupper(colnames(n1))
esets <- list(m = m.eset, n = n.eset)

test_that("DoppelGang class methods work (Normal Use)", {
  suppressWarnings(suppressMessages(
    res <- doppelgangR(esets, BPPARAM = BiocParallel::SerialParam())
  ))
  
  # Correctness / Normal Use: Test class
  expect_s4_class(res, "DoppelGang")
  
  # Normal Use: Test show method
  output_show <- capture.output(show(res))
  expect_true(any(grepl("S4 object of class: DoppelGang", output_show)))
  expect_true(any(grepl("Number of potential doppelgangers:", output_show)))
  
  # Normal Use: Test print method
  output_print <- capture.output(print(res))
  # Should print a data.frame with the summary results
  expect_true(any(grepl("sample1", output_print)))
  expect_true(any(grepl("sample2", output_print)))
  
  # Normal Use: Test summary method
  summ <- summary(res)
  expect_s3_class(summ, "data.frame")
  expect_equal(nrow(summ), nrow(res@summaryresults))
  
  # Normal Use: Test plot method
  # Save plot to pdf to avoid opening window
  tmp <- tempfile()
  pdf(tmp)
  on.exit({ dev.off(); unlink(tmp) }, add = TRUE)
  expect_no_error(plot(res))
})

test_that("DoppelGang plot method skip.no.doppels and plot.pair arguments work (Edge Cases)", {
  suppressWarnings(suppressMessages(
    res <- doppelgangR(esets, BPPARAM = BiocParallel::SerialParam())
  ))
  
  tmp <- tempfile()
  pdf(tmp)
  on.exit({ dev.off(); unlink(tmp) }, add = TRUE)
  # Test skip.no.doppels
  expect_no_error(plot(res, skip.no.doppels = TRUE))
  
  # Test plot.pair
  expect_no_error(plot(res, plot.pair = c("m", "n")))
  
  # Error handling: plot.pair with wrong names
  expect_error(plot(res, plot.pair = c("m", "q")), "One or both of plot.pair do not match names")
  
  # Error handling: plot.pair with wrong length
  expect_error(plot(res, plot.pair = c("m")), "plot.pair must be a character vector of length two")
})
