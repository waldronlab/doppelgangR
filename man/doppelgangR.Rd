\name{doppelgangR}
\alias{doppelgangR}
\title{doppelgangR}
\description{Identify samples with suspiciously high correlations and phenotype similarities}
\usage{doppelgangR(esets, separator = ":", corFinder.args = list(separator = separator, 
    use.ComBat = TRUE, method = "pearson"), phenoFinder.args = list(separator = separator), 
    phenoDist.args = list(vectorDistFun = vectorHammingDist), 
    outlierFinder.expr.args = list(bonf.pvalue = 0.05, transFun = atanh), 
    outlierFinder.pheno.args = list(normal.upper.thresh = 0.99, 
        bonf.pvalue = NULL), smokingGunFinder.args = list(transFun = I), 
    manual.smokingguns = NULL, automatic.smokingguns = TRUE)}
\arguments{
  \item{esets}{a list of ExpressionSets, containing the numeric and phenotypic data to be analyzed.}
  \item{separator}{a delimitor to use between dataset names and sample names}
  \item{corFinder.args}{a list of arguments to be passed to the corFinder function}
  \item{phenoFinder.args}{a list of arguments to be passed to the phenoFinder function}
  \item{phenoDist.args}{a list of arguments to be passed to the phenoDist function}
  \item{outlierFinder.expr.args}{a list of arguments to be passed to outlierFinder when called for expression data}
  \item{outlierFinder.pheno.args}{a list of arguments to be passed to outlierFinder when called for phenotype data}
  \item{smokingGunFinder.args}{a list of arguments to be passed to smokingGunFinder}
  \item{manual.smokingguns}{a character vector of phenoData columns that, if identical, will
be considered evidence of duplication}
  \item{automatic.smokingguns}{automatically look for "smoking guns."  If TRUE, look for
phenotype variables that are unique to each patient in dataset 1,
also unique to each patient in dataset 2, but contain exact
matches between datasets 1 and 2.}
}

\value{Returns a dataframe with columns: sample1, sample2, correlation,
phenotype.similarity, numeric.doppel, pheno.doppel, smokinggun,
sample1.keep, sample2.keep.}

\author{Levi Waldron and Markus Riester}




\examples{
library(curatedOvarianData)
data(GSE32062.GPL6480_eset)
data(GSE32063_eset)
data(GSE12470_eset)
testesets <- list(JapaneseA=GSE32062.GPL6480_eset, JapaneseB=GSE32063_eset, Yoshihara2010=GSE12470_eset)
doppelgangR(testesets, corFinder.args=list(use.ComBat=FALSE, method="pearson"))
}
