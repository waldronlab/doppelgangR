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
library(curatedOvarianData)
data(GSE32063_eset)
data(GSE17260_eset)
esets2 <- list(JapaneseB=GSE32063_eset,
                Yoshihara2010=GSE17260_eset)

esets2 <- lapply(esets2, function(X){
    ## standardize the sample ids to improve matching based on clinical annotation
    X$alt_sample_name <- paste(X$sample_type, gsub("[^0-9]", "", X$alt_sample_name), sep="_")
    pData(X) <- pData(X)[, !grepl("uncurated_author_metadata", colnames(pData(X)))]
    X <- X[, 1:20]  ##speed computations
    return(X) })

phenoFinder(esets2)
}

