library(doppelgangR)

test_that(".outer2df inner logic handles edge cases and parameters correctly", {
  x <- c("A", "B", "C")
  y <- c("D", "E", "F")
  
  # Error handling: Not both vectors, and not a matrix
  expect_error(doppelgangR:::.outer2df(data.frame(x=1), y="B"), "Require either x to be a matrix")
  
  # Edge Case: bidirectional = FALSE, diag = FALSE
  res_ff <- doppelgangR:::.outer2df(x, y, bidirectional = FALSE, diag = FALSE)
  expect_equal(nrow(res_ff), 3) # upper triangle of 3x3 matrix has 3 elements
  
  # Edge Case: bidirectional = FALSE, diag = TRUE
  res_ft <- doppelgangR:::.outer2df(x, y, bidirectional = FALSE, diag = TRUE)
  expect_equal(nrow(res_ft), 6) # upper triangle + diag of 3x3 matrix has 6 elements
  
  # Edge Case: bidirectional = TRUE, diag = FALSE
  res_tf <- doppelgangR:::.outer2df(x, y, bidirectional = TRUE, diag = FALSE)
  expect_equal(nrow(res_tf), 6) # upper + lower triangle has 6 elements
})
