library(doppelgangR)

test_that("phenoDist handles vector inputs (Normal Use)", {
  x <- c(1, 2, 3)
  y <- c(1, 4, 3)
  dist_res <- phenoDist(x, y)
  expect_true(is.numeric(dist_res))
  expect_equal(length(dist_res), 1)
})

test_that("phenoDist .discretizeDataFrame inner logic handles bins correctly", {
  # We test the internal discretization indirectly by passing numeric dataframes
  # that trigger the .discretizeRow logic when distinct values > bins
  
  # DataFrame with more levels than bins (bins=2)
  set.seed(123)
  df1 <- data.frame(
    var1 = 1:10,       # 10 levels, > 2 bins -> cut() logic
    var2 = rep(1:2, 5) # 2 levels, <= 2 bins -> as.factor() logic
  )
  rownames(df1) <- paste0("Sample", 1:10)
  
  # Compare it to itself to trigger single-argument logic
  dist_res1 <- phenoDist(df1, bins = 2)
  expect_true(is.matrix(dist_res1))
  expect_equal(dim(dist_res1), c(10, 10))
  
  # DataFrame 2
  df2 <- data.frame(
    var1 = 11:20,
    var2 = rep(1:2, 5)
  )
  rownames(df2) <- paste0("OtherSample", 1:10)
  
  # Compare df1 to df2 to trigger two-argument logic with continuous variables
  dist_res2 <- phenoDist(df1, df2, bins = 2)
  expect_true(is.matrix(dist_res2))
  expect_equal(dim(dist_res2), c(10, 10))
})
