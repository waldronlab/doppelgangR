library(doppelgangR)

test_that("Internal nonexports math and plot functions run without error", {
  
  # Plotting functions
  x <- seq(-3, 3, length=10)
  y <- seq(-3, 3, length=10)
  
  expect_error(doppelgangR:::dsn2.plot(x, y, dp=list(xi=c(0,0), Omega=diag(2), alpha=c(0,0))), NA)
  expect_error(doppelgangR:::dst2.plot(x, y, dp=list(xi=c(0,0), Omega=diag(2), alpha=c(0,0), df=5)), NA)

  # Probability functions
  expect_type(doppelgangR:::pmsn(c(0,0), dp=list(xi=c(0,0), Omega=diag(2), alpha=c(0,0))), "double")
  expect_type(doppelgangR:::pmst(c(0,0), dp=list(xi=c(0,0), Omega=diag(2), alpha=c(0,0), df=5)), "double")

  # msn.quantities
  expect_type(doppelgangR:::msn.quantities(dp=list(xi=c(0,0), Omega=diag(2), alpha=c(0,0))), "list")
  
  # st.cumulants
  expect_type(doppelgangR:::st.cumulants(dp=c(0, 1, 0, 5), n=4), "double")
  
  # st.SFscore
  expect_type(doppelgangR:::st.SFscore(shape=0, df=5, z=0), "double")

  # T.Owen
  expect_type(doppelgangR:::T.Owen(1, 1), "double")
  
  # msn.mle and msn.fit
  set.seed(1)
  X <- matrix(1, nrow=100, ncol=1)
  y <- matrix(rnorm(200), ncol=2)
  expect_type(doppelgangR:::msn.mle(y=y), "list")
  expect_type(doppelgangR:::msn.fit(X=X, y=y), "list")

  # st.mmle
  expect_type(doppelgangR:::st.mmle(X=X, y=rnorm(100), df=5), "list")
  # sn.mle
  expect_type(doppelgangR:::sn.mle(y=rnorm(100)), "list")
})
