\name{DoppelGang-class}
\docType{class}
\alias{DoppelGang-class}
\alias{DoppelGang}
\alias{show,DoppelGang-method}
\alias{print,DoppelGang-method}
\alias{summary,DoppelGang-method}
\title{"DoppelGang-class and related methods"}
\description{An object returned from function \code{\link{doppelgangR}}
  containing likely doppelgangers, pairwise sample correlations between
  each pair of datasets, and "smoking gun" variables for each pair of
  datasets.}
\section{Slots}{
	 \describe{
    \item{\code{results}:}{ A dataframe of doppelganger search results.
    It may contain columns: sample1 (label), sample2 (label),
    expr.similarity (numeric), expr.doppel (TRUE for expression-evidence
  duplicates), pheno.similarity (numeric), pheno.doppel (TRUE for
  phenotype-evidence duplicates), smokinggun.similarity (numeric),
  smokinggun.doppel (TRUE for duplicates based on "smoking gun" evidence.)}
    \item{\code{correlations}:}{Pairwise expression similarity matrix (by
      default correlations) for each pair of datasets and within each dataset.}
    \item{\code{smokingguns}:}{Column names of pData that were specified
      or inferred to be "smoking gun" phenotypes: where an identical
      value implies duplicated samples.}
  }
}
\section{Methods}{
 \itemize{
   \item{print}{Use \code{print(DoppelGang-object)} for brief information.}
   \item{summary}{Use \code{summary(DoppelGang-object)} for brief information.}
   \item{show}{Use \code{DoppelGang-object} for brief information.}
 }
}

\section{Details}{
  Note that "smoking guns" refers to phenotype labels, such as sample
  IDs, that should be unique.  If these labels are duplicated between
  two samples, they are flagged under the smokinggun.doppel column of
  the results slot (e.g. output@results$smokinggun.doppel).
}
\author{ Levi Waldron }

\seealso{\code{\link{doppelgangR}}}
