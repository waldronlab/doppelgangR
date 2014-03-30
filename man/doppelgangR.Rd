\name{doppelgangR}
\alias{doppelgangR}
\title{doppelgangR}
\description{Identify samples with suspiciously high correlations and phenotype similarities}
\usage{doppelgangR(esets, separator = ":", corFinder.args = list(separator = separator, 
    use.ComBat = TRUE, method = "pearson"), phenoFinder.args = list(separator = separator), 
    phenoDist.args = list(vectorDistFun = vectorWeightedDist), 
    outlierFinder.expr.args = list(bonf.pvalue = 0.005, transFun = atanh, 
        tail = "upper"), outlierFinder.pheno.args = list(normal.upper.thresh = 0.99, 
        bonf.pvalue = NULL, tail = "upper"), smokingGunFinder.args = list(transFun = I), 
    impute.knn.args = list(k = 10, rowmax = 0.5, colmax = 0.8, 
        maxp = 1500, rng.seed = 362436069), manual.smokingguns = NULL, 
    automatic.smokingguns = TRUE, within.datasets.only = FALSE, 
    cache.dir = "cache", verbose = TRUE)}
\arguments{
  \item{esets}{a list of ExpressionSets, containing the numeric and phenotypic data to be analyzed.}
  \item{separator}{a delimitor to use between dataset names and sample names}
  \item{corFinder.args}{a list of arguments to be passed to the corFinder function.}
  \item{phenoFinder.args}{a list of arguments to be passed to the phenoFinder function.  If
NULL, samples with similar phenotypes will not be searched for.}
  \item{phenoDist.args}{a list of arguments to be passed to the phenoDist function}
  \item{outlierFinder.expr.args}{a list of arguments to be passed to outlierFinder when called for expression data}
  \item{outlierFinder.pheno.args}{a list of arguments to be passed to outlierFinder when called for phenotype data}
  \item{smokingGunFinder.args}{a list of arguments to be passed to smokingGunFinder}
  \item{impute.knn.args}{a list of arguments to be passed to impute::impute.knn.  Set to
NULL to do no knn imputation.}
  \item{manual.smokingguns}{a character vector of phenoData columns that, if identical, will
be considered evidence of duplication}
  \item{automatic.smokingguns}{automatically look for "smoking guns."  If TRUE, look for
phenotype variables that are unique to each patient in dataset 1,
also unique to each patient in dataset 2, but contain exact
matches between datasets 1 and 2.}
  \item{within.datasets.only}{If TRUE, only search within each dataset for doppelgangers.}
  \item{cache.dir}{The name of a directory in which to cache or look up results to save
re-calculating correlations.  Set to NULL for no caching.}
  \item{verbose}{Print progress information}
}

\value{Returns an object of S4-class "DoppelGang".  See ?DoppelGang-class.}

\author{Levi Waldron and Markus Riester}




\examples{
library(doppelgangR)
library(curatedOvarianData)

data(GSE32062.GPL6480_eset)
data(GSE32063_eset)
data(GSE12470_eset)
data(GSE17260_eset)

testesets <- list(JapaneseA=GSE32062.GPL6480_eset,
                  JapaneseB=GSE32063_eset,
                  Yoshihara2009=GSE12470_eset,
                  Yoshihara2010=GSE17260_eset)

testesets <- lapply(testesets, function(X){
    # standardize the sample ids to improve matching based on clinical annotation
    sampleNames(X) <- make.names(paste(X$sample_type,
        gsub("\\\\D","",X$alt_sample_name), sep="_"))
    X$alt_sample_name <- sampleNames(X)
    pData(X) <- pData(X)[, !grepl("uncurated_author_metadata", colnames(pData(X)))]
    X })

results1 <- doppelgangR(testesets, corFinder.args=list(use.ComBat=TRUE), cache.dir=NULL)
summary(results1)
plot(results1)
## Set phenoFinder.args to ignore similar phenotypes:
results2 <- doppelgangR(testesets, corFinder.args=list(use.ComBat=FALSE), phenoFinder.args=NULL, cache.dir=NULL)
summary(results2)
}
