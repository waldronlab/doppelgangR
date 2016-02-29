\name{phenoFinder}
\alias{phenoFinder}
\title{Calculate pairwise similarities of phenoData between samples for a list containing two ExpressionSets}
\description{This function acts as a wrapper to phenoDist to handle cases of
one ExpressionSet, a list of two identical ExpressionSets, or a
list of two different ExpressionSets.}
\usage{phenoFinder(eset.pair, separator = ":", ...)}
\arguments{
  \item{eset.pair}{input: a list of ExpressionSets with two elements, or an
ExpressionSet.  If the two elements are identical, return the
correlation matrix for pairs of samples in the first element.  If
not identical, return pairs between the two elements.}
  \item{separator}{a separator between dataset name (taken from the list names) and
sample name (taken from sampleNames(eset), to keep track of which
samples come from which dataset.}
  \item{\dots}{Extra arguments passed on to phenoDist}
}

\value{A matrix of similarities between the phenotypes of pairs of samples.}

\author{Levi Waldron, Markus Riester, Marcel Ramos}

\examples{
example("corFinder")

## Create pData datasets
m.pdata <- matrix(letters[sample(1:26, size=110, replace=TRUE)], ncol=10)
n.pdata <- matrix(letters[sample(1:26, size=100, replace=TRUE)], ncol=10)

## Induce phenotype doppelgangers
n.pdata[10, ] <- m.pdata[10, ]

n.eset$id <- toupper(colnames(n1))
m.eset$id <- toupper(colnames(m1))
pData(m.eset) <- data.frame(c(pData(m.eset), data.frame(m.pdata)), stringsAsFactors = FALSE)

n.eset$id[4] <- "gotcha"
m.eset$id[5] <- "gotcha"
pData(n.eset) <- data.frame(c(pData(n.eset), data.frame(n.pdata)))

phenoFinder(esets)
}

