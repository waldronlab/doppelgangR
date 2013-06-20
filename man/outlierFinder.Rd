\name{outlierFinder}
\alias{outlierFinder}
\title{##Identifies outliers in a similarity matrix. }
\description{By default uses the
Fisher z-transform for Pearson correlation (atanh), and
identifies outliers as those above the quantile of a normal
distribution with mean and standard deviation estimated from the
z-transformed matrix.  The quantile is calculated from the
Bonferroni-corrected cumulative probability of the upper tail.}
\usage{outlierFinder(similarity.mat, bonf.pvalue = 0.05, transFun = atanh, 
    normal.upper.thresh = NULL, prune.output = TRUE)}
\arguments{
  \item{similarity.mat}{A matrix of similarities - larger values mean more similar.}
  \item{bonf.pvalue}{Bonferroni-corrected p-value.  A raw.pvalue is calculated by
dividing this by the number of non-missing values in
similarity.mat, and the rejection threshold is qnorm(1-raw.pvalue,
mean, sd) where mean and sd are estimated from the
transFun-transformed similarity.mat.}
  \item{transFun}{A function applied to the numeric values of similarity.mat, that
should result in normally-distributed values.}
  \item{normal.upper.thresh}{Instead of specifying bonf.pvalue and transFun, an upper
similarity threshold can be set, and values above this will be
considered likely duplicates.}
  \item{prune.output}{If prune.output=TRUE, only return likely doppelgangers.}
}

\value{Returns either NULL or a dataframe with three columns: sample1, sample2, and similarity.}

\author{Levi Waldron and Markus Riester}





