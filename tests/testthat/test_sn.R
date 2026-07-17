library(doppelgangR)

test_that("Skew-t distribution functions work correctly", {
  # Test dst
  x <- c(-1, 0, 1)
  expect_type(dst(x, location=0, scale=1, shape=0, df=Inf), "double")
  expect_type(dst(x, location=0, scale=1, shape=0, df=5), "double")
  expect_type(dst(x, dp=c(0, 1, 0, Inf)), "double")
  expect_type(dst(x, location=0, scale=1, shape=0, df=5, log=TRUE), "double")
  expect_error(dst(x, location=0, dp=c(0, 1, 0, Inf), shape=0), "You cannot set both component parameters and dp")

  # Test pst
  expect_type(pst(x, location=0, scale=1, shape=0, df=Inf), "double")
  expect_type(pst(x, location=0, scale=1, shape=0, df=5), "double")
  expect_type(pst(x, dp=c(0, 1, 0, Inf)), "double")
  expect_error(pst(x, location=0, dp=c(0, 1, 0, Inf), shape=0), "You cannot set both component parameters and dp")

  # Test qst
  p <- c(0.1, 0.5, 0.9)
  expect_type(qst(p, location=0, scale=1, shape=0, df=Inf), "double")
  expect_type(qst(p, location=0, scale=1, shape=0, df=5), "double")
  expect_type(qst(p, dp=c(0, 1, 0, Inf)), "double")
  expect_error(qst(p, location=0, dp=c(0, 1, 0, Inf), shape=0), "You cannot set both component parameters and dp")

  # Test rst
  expect_type(rst(5, location=0, scale=1, shape=0, df=Inf), "double")
  expect_type(rst(5, location=0, scale=1, shape=0, df=5), "double")
  expect_type(rst(5, dp=c(0, 1, 0, Inf)), "double")
  expect_error(rst(5, location=0, dp=c(0, 1, 0, Inf), shape=0), "You cannot set both component parameters and dp")
})

test_that("mst.mle works", {
  # Generate some mock data for mst.mle
  set.seed(1)
  X <- matrix(rnorm(200), ncol=2)
  # Test mst.mle execution (this is complex, but let's try a simple fit)
  # doppelgangR's code uses st.mle internally for corFinder.
  res <- st.mle(y = X[, 1])
  expect_type(res, "list")
  
  # mst.mle is exported, so we test it via doppelgangR::mst.mle()
  res_mst <- doppelgangR::mst.mle(y = X)
  expect_type(res_mst, "list")
  
  res_mst_fixed <- doppelgangR::mst.mle(y = X, fixed.df = 5)
  expect_type(res_mst_fixed, "list")

  # Test with freq
  res_mst_w <- doppelgangR::mst.mle(y = X, freq = rep(1, 100))
  expect_type(res_mst_w, "list")
})

test_that("st.mle works with various arguments", {
  set.seed(1)
  y <- rnorm(100)
  X <- matrix(1, nrow=100, ncol=1)
  
  # Cover the trace = TRUE argument block in st.mle
  expect_error(doppelgangR:::st.mle(X=X, y=y, trace=TRUE), NA)
  
  # Test fixed.df
  expect_error(doppelgangR:::st.mle(X=X, y=y, fixed.df=5), NA)
})
