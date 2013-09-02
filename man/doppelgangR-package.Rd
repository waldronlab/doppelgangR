\name{doppelgangR-package}
\alias{doppelgangR-package}
\alias{doppelgangR}
\docType{package}
\title{Identify likely duplicate samples from genomic or meta-data.}
\description{The main function is doppelgangR(), which takes as minimal
input a list of ExpressionSet object, and searches all list pairs for
duplicated samples.  The search is based on the genomic data (exprs(eset)),
phenotype/clinical data (pData(eset)), and "smoking guns" - supposedly unique
identifiers found in pData(eset).}
\details{
\tabular{ll}{Package: \tab doppelgangR\cr
Maintainer: \tab Levi Waldron <levi@jimmy.harvard.edu>\cr
Depends: \tab R (>= 2.10), affy, methods, Biobase\cr
Suggests: \tab curatedOvarianData\cr
Imports: \tab sva, impute\cr
Author: \tab Levi Waldron and Markus Riester\cr
Version: \tab 0.5.0\cr
Date: \tab 2013-09-02\cr
License: \tab GPL (>=2.0)\cr
Title: \tab Identify likely duplicate samples from genomic or meta-data.\cr
URL: \tab https://github.com/lwaldron/doppelgangR\cr
BugReports: \tab https://github.com/lwaldron/doppelgangR/issues\cr}
}
\author{Levi Waldron and Markus Riester}

\keyword{ package }


