\name{vectorWeightedDist}
\alias{vectorWeightedDist}
\title{Calculate a weighted distance between two vectors, using pairwise complete observations.}
\description{Simple function to count the fraction of different elements (in
the same position) between two vectors of the same length, after
removing elements from both vectors corresponding to positions
that are NA in either vector. Distance is the probability for observing
the matches and mismatches in two random patients.}
\usage{vectorWeightedDist(x, y, k, l)}
\arguments{
  \item{x}{a matrix }
  \item{y}{a matrix with the same number of columns as x}
  \item{k}{row in x to test for differences}
  \item{l}{row in y to test for differences}
}

\value{Returns a numeric value, the log of the probability of observing the
matches in x and y}

\author{Levi Waldron, Markus Riester, Marcel Ramos}


\examples{
mymat1 <- matrix(rnorm(20), ncol = 5)
mymat1[1, 4] <- NA
mymat2 <- matrix(rnorm(20), ncol = 5)
vectorWeightedDist(mymat1, mymat2, 1, 2)
}


