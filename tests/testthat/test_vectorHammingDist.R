library(doppelgangR)

test_that("vectorHammingDist works with no NAs (Normal Use)", {
  mat1 <- matrix(c("A", "B", "C", "D", "E"), nrow = 1)
  mat2 <- matrix(c("A", "B", "X", "Y", "E"), nrow = 1)
  
  res <- vectorHammingDist(mat1, mat2, 1, 1)
  # 2 differences out of 5 elements -> 0.4
  expect_equal(res, 0.4)
})

test_that("vectorHammingDist works with identical vectors (Correctness)", {
  mat1 <- matrix(c("A", "B", "C", "D", "E"), nrow = 1)
  
  res <- vectorHammingDist(mat1, mat1, 1, 1)
  # 0 differences
  expect_equal(res, 0)
})

test_that("vectorHammingDist handles NAs appropriately (Edge Cases)", {
  mat1 <- matrix(c("A", "B", NA, "D", "E"), nrow = 1)
  mat2 <- matrix(c("A", "B", "X", "Y", NA), nrow = 1)
  
  res <- vectorHammingDist(mat1, mat2, 1, 1)
  # Positions:
  # 1: A vs A -> Match
  # 2: B vs B -> Match
  # 3: NA vs X -> Ignored
  # 4: D vs Y -> Diff
  # 5: E vs NA -> Ignored
  # So we have 3 valid comparisons, 1 difference -> 1/3
  expect_equal(res, 1/3)
})

test_that("vectorHammingDist handles matrices with multiple rows (Normal Use)", {
  mat1 <- matrix(c(
    "A", "B", "C",
    "D", "E", "F"
  ), nrow = 2, byrow = TRUE)
  
  mat2 <- matrix(c(
    "A", "X", "C",
    "Y", "E", "Z"
  ), nrow = 2, byrow = TRUE)
  
  # Row 1 vs Row 1: A,B,C vs A,X,C -> 1 diff out of 3 = 1/3
  expect_equal(vectorHammingDist(mat1, mat2, 1, 1), 1/3)
  
  # Row 2 vs Row 2: D,E,F vs Y,E,Z -> 2 diff out of 3 = 2/3
  expect_equal(vectorHammingDist(mat1, mat2, 2, 2), 2/3)
  
  # Row 1 vs Row 2: A,B,C vs Y,E,Z -> 3 diff out of 3 = 1
  expect_equal(vectorHammingDist(mat1, mat2, 1, 2), 1)
})
