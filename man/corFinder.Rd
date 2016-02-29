\name{corFinder}
\alias{corFinder}
\title{Calculate pair-wise correlations between samples using the expr() slots of a list of two ExpressionSets.}
\description{This function acts as a wrapper around ComBat (sva package) and
cor(), to calculate pairwise correlations within one or between
two ExpressionSets.}
\usage{corFinder(eset.pair, separator = ":", use.ComBat = TRUE, ...)}
\arguments{
  \item{eset.pair}{a list of ExpressionSets, with two elements.  If the two
elements are identical, return the correlation matrix for pairs
of samples in the first element.  If not identical, return pairs
between the two elements.}
  \item{separator}{Separator between dataset name and sample name.  Dataset names are
added to sample names to keep track of dataset of origin.}
  \item{use.ComBat}{Use the sva::ComBat function for batch correction of the expr()
data between the two datasets.}
  \item{\dots}{Extra arguments passed to the cor() function.}
}

\value{Returns a matrix of sample-wise Pearson Correlations.}

\author{Levi Waldron, Markus Riester, Marcel Ramos}

\examples{
m1 <- matrix(rnorm(1100), ncol=11,
             dimnames = list(make.names(1:100),
                             paste0("m", 1:11)))

n1 <- matrix(rnorm(1000), ncol=10,
             dimnames = list(make.names(1:100),
                             paste0("n", 1:10)))

m1[, 1] <- n1[, 1] + rnorm(100, sd=0.25)
m1[, 2] <- m1[, 3] + rnorm(100, sd=0.25)
n1[, 2] <- n1[, 3] + rnorm(100, sd=0.25)
n1[, 4] <- m1[, 6] + rnorm(100, sd=0.25)

m.eset <- ExpressionSet(assayData=m1)
n.eset <- ExpressionSet(assayData=n1)

esets <- list(m=m.eset, n=n.eset)
corFinder(esets)
}



