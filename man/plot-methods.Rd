\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,DoppelGang-method}
\alias{plot.DoppelGang}
\alias{plot.doppelgangR}
\title{ Histograms of all pairwise sample correlations, showing
  identified doppelgangers.}
\description{Identified doppelgangers are shown with a red vertical line
  overlaid on a histogram of pairwise sample correlations.  One plot is
  made per pair of datasets.}
\arguments{
  \item{x}{An object of class \code{\link{DoppelGang}}}
  \item{skip.no.doppels}{(default FALSE) If TRUE, do not plot histograms where no
    doppelgangers were identified.}
  \item{plot.pair}{An optional character vector of length two, providing
    the names of two datasets.  If provided, only the comparison of
    these two datasets will be plotted.}
  \item{\dots}{ Additional arguments passed on to \code{\link{hist}}. }
}
\value{None}
\examples{
library(curatedOvarianData)
data(TCGA_eset)
data(GSE26712_eset)
testesets <- list(Bonome08=GSE26712_eset, TCGA=TCGA_eset)
library(doppelgangR)
results1 <- doppelgangR(testesets, corFinder.args=list(use.ComBat=FALSE), phenoFinder.args=NULL)
plot(results1)
}
\author{Levi Waldron}
\section{Methods}{
\describe{
\item{\code{signature(x = "DoppelGang")}}{
  Histograms of all pairwise sample correlations, showing
  identified doppelgangers.
}
}}
\keyword{methods}

