\name{phenoDist}
\alias{phenoDist}
\title{##Calculate distance between two vectors, rows of one matrix/dataframe, or rows of two matrices/dataframes.}

\usage{phenoDist(x, y = NULL, vectorDistFun = vectorHammingDist, ...)}
\arguments{
  \item{x}{A vector, matrix or dataframe}
  \item{y}{NULL, a vector, matrix, or dataframe.  If x is a vector, y must also be specified.}
  \item{vectorDistFun}{A function of two vectors that returns the distance between those vectors.}
  \item{\dots}{Extra arguments passed on to vectorDistFun}
}

\value{a matrix of distances between pairs of rows of x (if y is
unspecified), or between all pairs of rows between x and y (if
both are provided).}

\author{Levi Waldron and Markus Riester}





