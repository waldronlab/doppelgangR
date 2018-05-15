#' DoppelGang S4 class
#' 
#' S4 class containing results of doppelgangR() function.
#' 
#' 
#' @name DoppelGang-class
#' @aliases DoppelGang-class DoppelGang
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new(DoppelGang ...)}
#' @param x,object A DoppelGang class object
#' @author Levi Waldron and Markus Riester
#' @seealso \code{\link{plot,DoppelGang-method}}
#' @export
setClass(
  Class = "DoppelGang",
  representation(
    fullresults = "list",
    summaryresults = "data.frame",
    inputargs = "list"
  )
)

setGeneric("print")
setGeneric("summary")
setGeneric("plot")

#' @name DoppelGang-class
#' @aliases summary,DoppelGang-method
#' @export
setMethod("summary", signature(object = "DoppelGang"),
          function(object)
            object@summaryresults)
 
#' @name DoppelGang-class
#' @aliases show,DoppelGang-method
#' @export
setMethod("show", signature(object = "DoppelGang"),
          function(object) {
            cat ("S4 object of class:" , class (object) , "\n")
            cat (
              "Number of potential doppelgangers:" ,
              nrow(object@summaryresults),
              ": ",
              sum(object@summaryresults$expr.doppel),
              "expression, ",
              sum(object@summaryresults$pheno.doppel),
              "phenotype, ",
              sum(object@summaryresults$smokinggun.doppel),
              "smoking gun. \n"
            )
            cat (" \n ")
            cat ("Use summary(object) to obtain a data.frame of potential doppelgangrs.",
                 "\n")
            cat (" \n ")
          })

#' @name DoppelGang-class
#' @aliases print,DoppelGang-method
#' @export
setMethod("print", signature(x = "DoppelGang"),
          function(x)
            print(x@summaryresults))

#' Histograms of all pairwise sample correlations, showing identified
#' doppelgangers.
#' 
#' Identified doppelgangers are shown with a red vertical line overlaid on a
#' histogram of pairwise sample correlations.  One plot is made per pair of
#' datasets.
#' 
#' 
#' @name plot-methods
#' @aliases plot-methods plot,DoppelGang plot,DoppelGang-method
#' plot,DoppelGang,ANY-method plot.DoppelGang plot.doppelgangR
#' @docType methods
#' @param x An object of class \code{\link{DoppelGang}}
#' @param skip.no.doppels (default FALSE) If TRUE, do not plot histograms where
#' no doppelgangers were identified.
#' @param plot.pair An optional character vector of length two, providing the
#' names of two datasets.  If provided, only the comparison of these two
#' datasets will be plotted.
#' @param \dots Additional arguments passed on to \code{\link{hist}}.
#' @return None
#' @section Methods: \describe{ \item{list("signature(x = \"DoppelGang\")")}{
#' Histograms of all pairwise sample correlations, showing identified
#' doppelgangers. } }
#' @author Levi Waldron
#' @keywords methods
#' @examples
#' 
#' library(curatedOvarianData)
#' data(TCGA_eset)
#' data(GSE26712_eset)
#' ## Remove some TCGA samples to speed computation:
#' keep.tcga <-
#' c("TCGA.13.2060", "TCGA.24.2290", "TCGA.25.2392", "TCGA.25.2404",
#' "TCGA.59.2349", "TCGA.09.2044", "TCGA.24.2262", "TCGA.24.2293",
#' "TCGA.25.2393", "TCGA.25.2408", "TCGA.59.2350", "TCGA.09.2045",
#' "TCGA.24.2267", "TCGA.59.2351", "TCGA.09.2048", "TCGA.24.2271",
#' "TCGA.24.2298", "TCGA.25.2398", "TCGA.59.2354", "TCGA.09.2050",
#' "TCGA.24.2281", "TCGA.09.2051", "TCGA.29.2428", "TCGA.09.2055",
#' "TCGA.24.2289", "TCGA.29.2414", "TCGA.59.2352", "TCGA.36.2532",
#' "TCGA.36.2529", "TCGA.36.2551", "TCGA.42.2590", "TCGA.13.2071",
#' "TCGA.29.2432", "TCGA.36.2537", "TCGA.36.2547", "TCGA.04.1369",
#' "TCGA.42.2591", "TCGA.23.2641", "TCGA.29.2434", "TCGA.36.2538",
#' "TCGA.36.2548", "TCGA.04.1516", "TCGA.42.2593", "TCGA.36.2549",
#' "TCGA.04.1644", "TCGA.13.2057", "TCGA.23.2647", "TCGA.36.2530",
#' "TCGA.36.2552", "TCGA.42.2587", "TCGA.13.2061", "TCGA.42.2588",
#' "TCGA.36.2544", "TCGA.42.2589", "TCGA.13.2066", "TCGA.61.2613",
#' "TCGA.61.2614", "TCGA.24.1852", "TCGA.29.1704", "TCGA.13.1819"
#' )
#' keep.tcga <- unique(c(keep.tcga, sampleNames(TCGA_eset)[1:200]))
#' testesets <- list(Bonome08=GSE26712_eset, TCGA=TCGA_eset[, keep.tcga])
#' results1 <- doppelgangR(testesets,
#'                         corFinder.args=list(use.ComBat=FALSE), phenoFinder.args=NULL, cache.dir=NULL)
#' plot(results1)
#' @export
setMethod("plot", signature(x = "DoppelGang"),
          function(x,
                   skip.no.doppels = FALSE,
                   plot.pair = NULL,
                   ...) {
            if (!is.null(plot.pair)) {
              if (is.character(plot.pair) & length(plot.pair) == 2) {
                study.names <-
                  unique(unlist(strsplit(
                    names(x@fullresults), split = x@inputargs$separator
                  )))
                if (!all(plot.pair %in% study.names))
                  stop("One or both of plot.pair do not match names(esets)")
                iplot <-
                  which(names(x@fullresults) %in% c(
                    paste(plot.pair, collapse = ":"),
                    paste(rev(plot.pair), collapse = ":")
                  ))
              } else{
                stop(
                  "plot.pair must be a character vector of length two, containing the names of two datasets seen in names(esets)"
                )
              }
            } else{
              iplot <- 1:length(x@fullresults)
            }
            op <- par()
            for (i in iplot) {
              cors <- x@fullresults[[i]]$correlations
              cors <- na.omit(as.numeric(cors))
              expr.doppels <-
                x@fullresults[[i]]$expr.doppels$outlierFinder.res
              expr.doppels <-
                expr.doppels[expr.doppels$doppel,]
              if (skip.no.doppels &
                  (nrow(expr.doppels) == 0 | is.null(expr.doppels)))
                next
              stfit <- x@fullresults[[i]]$expr.doppels$stfit
              hist.stats <-
                hist(
                  cors,
                  main = paste(names(x@fullresults)[i], stfit$algorithm$message, sep = "\n"),
                  freq = TRUE,
                  xlab = "Pairwise Correlations",
                  breaks = "FD",
                  ...
                )
              abline(v = expr.doppels$similarity,
                     col = "red",
                     lw = 0.5)
              xvals <- hist.stats$breaks
              transFun <-
                x@inputargs$outlierFinder.expr.args$transFun
              suppressWarnings(xvals <-
                                 xvals[!is.na(transFun(xvals))])
              dens.vals <-
                diff(
                  pst(
                    transFun(xvals),
                    location = stfit$dp["location"],
                    scale = stfit$dp["scale"],
                    shape = stfit$dp["shape"],
                    df = stfit$dp["df"]
                  )
                )
              xvals <- (xvals[-length(xvals)] + xvals[-1]) / 2
              freq.vals <- dens.vals * sum(hist.stats$counts)
              lines(xvals, freq.vals)
              par(ask = prod(par("mfcol")) < length(x@fullresults) &&
                    dev.interactive())
            }
            suppressWarnings(par(op))
          })
