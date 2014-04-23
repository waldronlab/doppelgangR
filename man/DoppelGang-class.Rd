\name{DoppelGang-class}
\Rdversion{1.1}
\docType{class}
\alias{DoppelGang-class}
\alias{DoppelGang}
%% \alias{plot}
%% \alias{print}
%% \alias{show}
%% \alias{summary}

\title{DoppelGang S4 class}
\description{S4 class containing results of doppelgangR() function.}
\section{Objects from the Class}{Objects can be created by calls of the form \code{new(DoppelGang ...)}}
\section{Slots}{
  \describe{
    \item{\code{fullresults}:}{Object of class \code{"list"}}
    \item{\code{summaryresults}:}{Object of class \code{"data.frame"}}
    \item{\code{inputargs}:}{\code{"list"} object containing input arguments to doppelgangR() function.}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "DoppelGang")}: See ?plot,DoppelGang-method}
    \item{print}{\code{signature(x = "DoppelGang")}: ... }
    \item{show}{\code{signature(object = "DoppelGang")}: ... }
    \item{summary}{\code{signature(object = "DoppelGang")}: ... }
	 }
}

\seealso{\code{\link{plot,DoppelGang-method}}}

\author{Levi Waldron and Markus Riester}





