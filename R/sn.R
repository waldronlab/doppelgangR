#' Skew-t Distribution
#' 
#' Density function, distribution function and random number generation for the
#' skew-\eqn{t} (ST) distribution.  Functions copied from \code{sn} CRAN
#' library v0.4.18 for argument name compatibility with \code{st.mle} function
#' from the same version.
#' 
#' 
#' @aliases dst pst qst rst
#' @param x vector of quantiles. Missing values (\code{NA}s) are allowed.
#' @param p vector of probabililities
#' @param location vector of location parameters.
#' @param scale vector of (positive) scale parameters.
#' @param shape vector of shape parameters. With \code{pst} and \code{qst}, it
#' must be of length 1.
#' @param df degrees of freedom (scalar); default is \code{df=Inf} which
#' corresponds to the skew-normal distribution.
#' @param dp a vector of length 4, whose elements represent location, scale
#' (positive), shape and df, respectively.  If \code{dp} is specified, the
#' individual parameters cannot be set.
#' @param n sample size.
#' @param log logical; if TRUE, densities are given as log-densities.
#' @param tol a scalar value which regulates the accuracy of the result of
#' \code{qsn}.
#' @param ... additional parameters passed to \code{integrate}.
#' @return Density (\code{dst}), probability (\code{pst}), quantiles
#' (\code{qst}) and random sample (\code{rst}) from the skew-\eqn{t}
#' distribution with given \code{location}, \code{scale}, \code{shape} and
#' \code{df} parameters.
#' @section Details: Typical usages are \preformatted{% dst(x, location=0,
#' scale=1, shape=0, df=Inf, log=FALSE) dst(x, dp=, log=FALSE) pst(x,
#' location=0, scale=1, shape=0, df=Inf, ...) pst(x, dp=, log=FALSE) qst(p,
#' location=0, scale=1, shape=0, df=Inf, tol=1e-8, ...) qst(x, dp=, log=FALSE)
#' rst(n=1, location=0, scale=1, shape=0, df=Inf) rst(x, dp=, log=FALSE) }
#' @seealso \code{\link{st.mle}}
#' @references Azzalini, A. and Capitanio, A. (2003). Distributions generated
#' by perturbation of symmetry with emphasis on a multivariate skew-\emph{t}
#' distribution. \emph{J.Roy. Statist. Soc. B} \bold{65}, 367--389.
#' @keywords distribution
#' @examples
#' 
#' pdf <- dst(seq(-4,4,by=0.1), shape=3, df=5)
#' rnd <- rst(100, 5, 2, -5, 8)
#' q <- qst(c(0.25,0.5,0.75), shape=3, df=5)
#' stopifnot(identical(all.equal(pst(q, shape=3, df=5), c(0.25,0.5,0.75)), TRUE))
#' 
#' @export dst
dst <-
  function (x,
            location = 0,
            scale = 1,
            shape = 0,
            df = Inf,
            dp = NULL,
            log = FALSE)
  {
    if (!is.null(dp)) {
      if (!missing(shape))
        stop("You cannot set both component parameters and dp")
      location <- dp[1]
      scale <- dp[2]
      shape <- dp[3]
      df <- dp[4]
    }
    if (df > 1e6)
      return(dsn(x, location, scale, shape, log = log))
    z   <- (x - location) / scale
    pdf <- dt(z, df = df, log = log)
    cdf <- pt(shape * z * sqrt((df + 1) / (z ^ 2 + df)), df = df + 1, log.p =
                log)
    if (log)
      logb(2) + pdf + cdf - logb(scale)
    else
      2 * pdf * cdf / scale
  }


rst <-
  function (n = 1,
            location = 0,
            scale = 1,
            shape = 0,
            df = Inf,
            dp = NULL)
  {
    if (!is.null(dp)) {
      if (!missing(shape))
        stop("You cannot set both component parameters and dp")
      location <- dp[1]
      scale <- dp[2]
      shape <- dp[3]
      df <- dp[4]
    }
    z <- rsn(n, location = 0, scale, shape)
    if (df == Inf)
      return(z + location)
    v <- rchisq(n, df) / df
    y <- z / sqrt(v) + location
    attr(y, "parameters") <- c(location, scale, shape, df)
    return(y)
  }

pst <-
  function (x,
            location = 0,
            scale = 1,
            shape = 0,
            df = Inf,
            dp = NULL,
            ...)
  {
    if (!is.null(dp)) {
      if (!missing(shape))
        stop("You cannot set both component parameters and dp")
      location <- dp[1]
      scale <- dp[2]
      shape <- dp[3]
      df <- dp[4]
    }
    fp <- function(v, shape, df, t.value)
      psn(sqrt(v) * t.value, 0, 1, shape) * dchisq(v * df, df = df) * df
    if (df > 1e6)
      # (== Inf)
      p <- psn(x, location, scale, shape)
    else
    {
      if (df <= 0)
        stop("df must be non-negative")
      z <- (x - location) / scale
      p <- numeric(length(z))
      for (i in 1:length(z)) {
        p[i]  <-
          if (round(df) == df)
            pmst(z[i], 0, matrix(1, 1, 1), shape, df, ...)
        else{
          if (abs(z[i]) == Inf)
            (1 + sign(z[i])) / 2
          else{
            if (z[i] < 0)
              integrate(dst,-Inf, z[i], shape = shape, df = df, ...)$value
            else
              integrate(
                fp,
                0,
                Inf,
                shape = shape,
                df = df,
                t.value = z[i],
                ...
              )$value
          }
        }
      }
      pmax(0, pmin(1, p))
    }
  }

qst <- function (p,
                 location = 0,
                 scale = 1,
                 shape = 0,
                 df = Inf,
                 tol = 1e-06,
                 dp = NULL,
                 ...)
{
  if (!is.null(dp)) {
    if (!missing(shape))
      stop("You cannot set both component parameters and dp")
    location <- dp[1]
    scale <- dp[2]
    shape <- dp[3]
    df <- dp[4]
  }
  if (df > 1e4)
    # (== Inf)
    return(qsn(p, location, scale, shape))
  max.q <- sqrt(qf(p, 1, df))
  min.q <- -sqrt(qf(1 - p, 1, df))
  if (shape == Inf)
    return(location + scale * max.q)
  if (shape == -Inf)
    return(location + scale * min.q)
  na <- is.na(p) | (p < 0) | (p > 1)
  zero <- (p == 0)
  one <- (p == 1)
  p <- replace(p, (na | zero | one), 0.5)
  cum <- st.cumulants(0, 1, shape, max(df, 5), n = 4)
  g1 <- cum[3] / cum[2] ^ (3 / 2)
  g2 <- cum[4] / cum[2] ^ 2
  x <- qnorm(p)
  x <- (x + (x ^ 2 - 1) * g1 / 6 + x * (x ^ 2 - 3) * g2 / 24 -
          x * (2 *  x ^ 2 - 5) * g1 ^ 2 / 36)
  x <- cum[1] + sqrt(cum[2]) * x
  max.err <- 1
  while (max.err > tol) {
    x1 <- x - (pst(x, 0, 1, shape, df, ...) - p) / dst(x, 0, 1, shape, df)
    x1 <- pmin(x1, max.q)
    x1 <- pmax(x1, min.q)
    max.err <- max(abs(x1 - x) / (1 + abs(x)))
    x <- x1
  }
  x <- replace(x, na, NA)
  x <- replace(x, zero,-Inf)
  x <- replace(x, one, Inf)
  return(as.numeric(location + scale * x))
}

st.mle <- function(X,
                   y,
                   freq,
                   start,
                   fixed.df = NA,
                   trace = FALSE,
                   algorithm = c("nlminb", "Nelder-Mead", "BFGS", "CG", "SANN"),
                   control = list())
{
  y.name  <- deparse(substitute(y))
  y <- data.matrix(y)
  if (missing(X))
    X <- matrix(1, nrow = length(y), ncol = 1)
  dimnames(y)[[2]] <- list(y.name)
  if (missing(start)) {
    cp0 <- sn.mle(
      X = X,
      y = y,
      plot.it = FALSE,
      trace = trace
    )$cp
    m <- length(cp0) - 2
    cp0[m + 2] <- cp0[m + 2] * 0.9
    mle0 <- cp.to.dp(cp0)
    start <- list(
      beta = mle0[1:m],
      Omega = matrix(mle0[m + 1] ^ 2, 1, 1),
      alpha = mle0[m + 2],
      df = 10
    )
  }
  else {
    m <- length(start) - 3
    if (m < 1)
      stop("bad start vector")
    start <-
      list(
        beta = start[1:m],
        Omega = matrix(start[m + 1] ^ 2, 1, 1),
        alpha = start[m + 2],
        df = start[m + 3]
      )
  }
  fit <-
    mst.mle(
      X,
      y,
      freq,
      start = start,
      fixed.df = fixed.df,
      trace = trace,
      algorithm = algorithm,
      control = control
    )
  mle <- list()
  mle$call <- match.call()
  dp <- fit$dp
  se <- fit$se
  p  <- length(dp$beta)
  dp.names <- c(if (p == 1)
    "location"
    else
      dimnames(dp$beta)[[1]],
    "scale", "shape", "df")
  mle$dp  <-
    c(dp$beta, sqrt(as.vector(dp$Omega)), dp$alpha, dp$df)
  names(mle$dp) <- dp.names
  mle$se <- if (all(is.na(se)))
    NA
  else
    c(se$beta, mle$dp[p + 1] * se$internal[p + 1],
      se$alpha, dp$df * se$internal[p + 3])
  mle$logL <- fit$logL
  mle$algorithm <- fit$algorithm
  mle
}




#' Maximum likelihood estimation for a (multivariate) skew-t distribution
#' 
#' Fits a skew-t (ST) or multivariate skew-t (MST) distribution to data, or
#' fits a linear regression model with (multivariate) skew-t errors, using
#' maximum likelihood estimation.  Functions copied from \code{sn} CRAN library
#' v0.4.18 because they were later deprecated in that library.
#' 
#' If \code{y} is a vector and it is supplied to \code{mst.mle}, then it is
#' converted to a one-column matrix, and a scalar skew-t distribution is
#' fitted. This is also the mechanism used by \code{st.mle} which is simply an
#' interface to \code{mst.mle}.
#' 
#' The parameter \code{freq} is intended for use with grouped data, setting the
#' values of \code{y} equal to the central values of the cells; in this case
#' the resulting estimate is an approximation to the exact maximum likelihood
#' estimate. If \code{freq} is not set, exact maximum likelihood estimation is
#' performed.
#' 
#' % To fit a scalar skew-t distribution to grouped data by exact % maximum
#' likelihood estimation, use \code{st.mle.grouped}.
#' 
#' Numerical search of the maximum likelihood estimates is performed in a
#' suitable re-parameterization of the original parameters with aid of the
#' selected optimizer (\code{nlminb} or \code{optim}) which is supplied with
#' the derivatives of the log-likelihood function. Notice that, in case the
#' optimizer is \code{optim}), the gradient may or may not be used, depending
#' on which specific method has been selected.  On exit from the optimizer, an
#' inverse transformation of the parameters is performed. For a specific
#' description on the re-parametrization adopted, see Section 5.1 and Appendix
#' B of Azzalini \& Capitanio (2003).
#' 
#' @aliases mst.mle st.mle
#' @param y a matrix (for \code{mst.mle}) or a vector (for \code{st.mle}). If
#' \code{y} is a matrix, rows refer to observations, and columns to components
#' of the multivariate distribution.
#' @param X a matrix of covariate values.  If missing, a one-column matrix of
#' 1's is created; otherwise, it must have the same number of rows of \code{y}.
#' If \code{X} is supplied, then it must include a column of 1's.
#' @param freq a vector of weights. If missing, a vector of 1's is created;
#' otherwise it must have length equal to the number of rows of \code{y}.
#' @param start for \code{mst.mle}, a list contaning the components
#' \code{beta},\code{Omega}, \code{alpha}, \code{df} of the type described
#' below; for \code{st.mle}, a vector whose components contain analogous
#' ingredients as before, with the exception that the scale parameter is the
#' square root of \code{Omega}.  In both cases, the \code{dp} component of the
#' returned list from a previous call has the required format and it can be
#' used as a new \code{start}. If the \code{start} parameter is missing,
#' initial values are selected by the function.
#' @param fixed.df a scalar value containing the degrees of freedom (df), if
#' these must be taked as fixed, or \code{NA} (default value) if \code{df} is a
#' parameter to be estimated.
#' @param trace logical value which controls printing of the algorithm
#' convergence. If \code{trace=TRUE}, details are printed. Default value is
#' \code{FALSE}.
#' @param algorithm a character string which selects the numerical optimization
#' procedure used to maximize the loglikelihood function. If this string is set
#' equal to \code{"nlminb"}, then this function is called; in all other cases,
#' \code{optim} is called, with \code{method} set equal to the given string.
#' Default value is \code{"nlminb"}.
#' @param control this parameter is passed to the chose optimizer, either
#' \code{nlminb} or \code{optim}; see the documentation of this function for
#' its usage.
#' @return A list containing the following components:
#' 
#' \item{call}{ a string containing the calling statement. } \item{dp}{ for
#' \code{mst.mle}, this is a list containing the direct parameters \code{beta},
#' \code{Omega}, \code{alpha}. Here, \code{beta} is a matrix of regression
#' coefficients with \code{dim(beta)=c(ncol(X),ncol(y))}, \code{Omega} is a
#' covariance matrix of order \code{ncol(y)}, \code{alpha} is a vector of shape
#' parameters of length \code{ncol(y)}.  For \code{st.mle}, \code{dp} is a
#' vector of length \code{ncol(X)+3}, containing \code{c(beta, omega, alpha,
#' df)}, where \code{omega} is the square root of \code{Omega}. } \item{se}{ a
#' list containing the components \code{beta}, \code{alpha}, \code{info}. Here,
#' \code{beta} and \code{alpha} are the standard errors for the corresponding
#' point estimates; \code{info} is the observed information matrix for the
#' working parameter, as explained below. } \item{algorithm}{ the list returned
#' by the chose optimizer, either \code{nlminb} or \code{optim}, plus an item
#' with the \code{name} of the selected algorithm; see the documentation of
#' either \code{nlminb} or \code{optim} for explanation of the other
#' components. }
#' @section Background: The family of multivariate skew-t distributions is an
#' extension of the multivariate Student's t family, via the introduction of a
#' \code{shape} parameter which regulates skewness; when \code{shape=0}, the
#' skew-t distribution reduces to the usual t distribution. When \code{df=Inf}
#' the distribution reduces to the multivariate skew-normal one; see
#' \code{dmsn}. See the reference below for additional information.
#' @seealso \code{\link{dst}}
#' @references Azzalini, A. and Capitanio, A. (2003).  Distributions generated
#' by perturbation of symmetry with emphasis on a multivariate skew \emph{t}
#' distribution.  The full version of the paper published in abriged form in
#' \emph{J.Roy. Statist. Soc. B} \bold{65}, 367--389, is available at
#' \url{http://azzalini.stat.unipd.it/SN/se-ext.ps}
#' @keywords distribution regression
#' @examples
#' 
#' dat <- rt(100, df=5, ncp=100)
#' fit <- st.mle(y=dat)
#' fit
#' 
#' @export mst.mle
mst.mle <-
  function (X,
            y,
            freq,
            start,
            fixed.df = NA,
            trace = FALSE,
            algorithm = c("nlminb", "Nelder-Mead", "BFGS", "CG", "SANN"),
            control = list())
  {
    algorithm <- match.arg(algorithm)
    y.name <- deparse(substitute(y))
    y.names <- dimnames(y)[[2]]
    y <- data.matrix(y)
    X <- if (missing(X))
      matrix(rep(1, nrow(y)), ncol = 1)
    else
      data.matrix(X)
    if (missing(freq))
      freq <- rep(1, nrow(y))
    x.names <- dimnames(X)[[2]]
    d <- ncol(y)
    n <- sum(freq)
    m <- ncol(X)
    if (missing(start)) {
      qrX <- qr(X)
      beta <- as.matrix(qr.coef(qrX, y))
      Omega <- matrix(var(qr.resid(qrX, y)), d, d)
      omega <- sqrt(diag(Omega))
      alpha <- rep(0, d)
      df <- ifelse(is.na(fixed.df), 10, fixed.df)
      if (trace) {
        cat("mst.mle: dp=", "\n")
        print(c(beta, Omega, alpha))
        cat("df:", df, "\n")
      }
    }
    else {
      if (!is.na(fixed.df))
        start$df <- fixed.df
      if (all(names(start) == c("beta", "Omega", "alpha", "df"))) {
        beta <- start$beta
        Omega <- start$Omega
        alpha <- start$alpha
        df <- start$df
      }
      else
        stop("start parameter is not in the form that I expected")
    }
    eta <- alpha / sqrt(diag(Omega))
    Oinv <- solvePD(Omega)
    upper <- chol(Oinv)
    D <- diag(upper)
    A <- upper / D
    D <- D ^ 2
    if (d > 1)
      param <-
      c(beta,-log(D) / 2, A[!lower.tri(A, diag = TRUE)], eta)
    else
      param <- c(beta,-log(D) / 2, eta)
    if (is.na(fixed.df))
      param <- c(param, log(df))
    if (algorithm == "nlminb") {
      opt <- nlminb(
        param,
        objective = mst.dev,
        gradient = mst.dev.grad,
        control = control,
        X = X,
        y = y,
        freq = freq,
        trace = trace,
        fixed.df = fixed.df
      )
      info <- num.deriv2(
        opt$par,
        FUN = "mst.dev.grad",
        X = X,
        y = y,
        freq = freq,
        fixed.df = fixed.df
      ) / 2
      opt$value <-  opt$objective
    }
    else{
      opt <- optim(
        param,
        fn = mst.dev,
        gr = mst.dev.grad,
        method = algorithm,
        control = control,
        hessian = TRUE,
        X = X,
        y = y,
        freq = freq,
        trace = trace,
        fixed.df = fixed.df
      )
      info <- opt$hessian / 2
    }
    dev   <- opt$value
    param <- opt$par
    opt$name <- algorithm
    if (trace) {
      cat("Message from optimization routine:", opt$message, "\n")
      cat("deviance:", dev, "\n")
    }
    beta <- matrix(param[1:(m * d)], m, d)
    D <- exp(-2 * param[(m * d + 1):(m * d + d)])
    A <- diag(d)
    i0 <- m * d + d * (d + 1) / 2
    if (d > 1)
      A[!lower.tri(A, diag = TRUE)] <- param[(m * d + d + 1):i0]
    eta <- param[(i0 + 1):(i0 + d)]
    if (is.na(fixed.df))
      df <- exp(param[i0 + d + 1])
    else
      df <- fixed.df
    Oinv <- t(A) %*% diag(D, d, d) %*% A
    Omega <- solvePD(Oinv)
    omega <- sqrt(diag(Omega))
    alpha <- eta * omega
    dimnames(beta) <- list(x.names, y.names)
    dimnames(Omega) <- list(y.names, y.names)
    if (length(y.names) > 0)
      names(alpha) <- y.names
    if (all(is.finite(info))) {
      qr.info <- qr(info)
      info.ok <- (qr.info$rank == length(param))
    }
    else
      info.ok <- FALSE
    if (info.ok) {
      se2 <- diag(solve(qr.info))
      if (min(se2) < 0)
        se <- NA
      else {
        se <- sqrt(se2)
        se.beta <- matrix(se[1:(m * d)], m, d)
        se.alpha <- se[(i0 + 1):(i0 + d)] * omega
        dimnames(se.beta)[2] <- list(y.names)
        dimnames(se.beta)[1] <- list(x.names)
        names(se.alpha) <- y.names
        se.df <- df * se[i0 + d + 1]
        se <- list(
          beta = se.beta,
          alpha = se.alpha,
          df = se.df,
          internal = se,
          info = info
        )
      }
    }
    else
      se <- NA
    dp <- list(
      beta = beta,
      Omega = Omega,
      alpha = alpha,
      df = df
    )
    list(
      call = match.call(),
      logL = -dev / 2,
      deviance = dev,
      dp = dp,
      se = se,
      algorithm = opt
    )
  }
