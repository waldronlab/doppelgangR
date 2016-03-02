\name{phenoDist}
\alias{phenoDist}
\title{Calculate distance between two vectors, rows of one matrix/dataframe, or rows of two matrices/dataframes.}
\description{This function does some simple looping to allow x and y to be
various combinations of vectors and matrices/dataframes.}
\usage{phenoDist(x, y = NULL, bins = 10, vectorDistFun = vectorWeightedDist, 
    ...)}
\arguments{
  \item{x}{A vector, matrix or dataframe}
  \item{y}{NULL, a vector, matrix, or dataframe.  If x is a vector, y must also be specified.}
  \item{bins}{discretize continuous fields in the specified number of bins}
  \item{vectorDistFun}{A function of two vectors that returns the distance between those vectors.}
  \item{\dots}{Extra arguments passed on to vectorDistFun}
}

\value{a matrix of distances between pairs of rows of x (if y is
unspecified), or between all pairs of rows between x and y (if
both are provided).}

\author{Levi Waldron, Markus Riester, Marcel Ramos}




\examples{
library(curatedOvarianData)
data(GSE32063_eset)
data(GSE17260_eset)
pdat1 <- pData(GSE32063_eset)
pdat2 <- pData(GSE17260_eset)
## Curation of the alternative sample identifiers makes duplicates stand out more:
pdat1$alt_sample_name <- paste(pdat1$sample_type, gsub("[^0-9]", "", pdat1$alt_sample_name), sep="_")
pdat2$alt_sample_name <- paste(pdat2$sample_type, gsub("[^0-9]", "", pdat2$alt_sample_name), sep="_")
## Removal of columns that cannot possibly match also helps duplicated patients to stand out
pdat1 <- pdat1[, !grepl("uncurated_author_metadata", colnames(pdat1))]
pdat2 <- pdat2[, !grepl("uncurated_author_metadata", colnames(pdat2))]
## Use phenoDist() to calculate a weighted distance matrix
distmat <- phenoDist(as.matrix(pdat1), as.matrix(pdat2))
## Note outliers with identical clinical data, these are probably the same patients:
graphics::boxplot(distmat)
}
