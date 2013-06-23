\name{vectorHammingDist}
\alias{vectorHammingDist}
\title{Calculate Hamming Distance between two vectors, using pairwise complete observations.}
\description{Simple function to count the fraction of different elements (in
the same position) between two vectors of the same length, after
removing elements from both vectors corresponding to positions
that are NA in either vector.}
\usage{vectorHammingDist(x, y)}
\arguments{
  \item{x}{a vector (normally character)}
  \item{y}{a vector the same length as x}
}

\value{Returns a numeric value, the Hamming Distance (the number of
non-equal values between x and y).}

\author{Levi Waldron and Markus Riester}





