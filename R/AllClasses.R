setClass(Class = "DoppelGang", representation(fullresults = "list", summaryresults="data.frame", inputargs="list"))

setGeneric("print")
setGeneric("summary")
setGeneric("plot")

setMethod("print", signature(x="DoppelGang"),
          function(x) print(x@summaryresults))

setMethod("summary", signature(object="DoppelGang"),
          function(object) summary(object@summaryresults[, 3:min(ncol(object@summaryresults), 8)]))

setMethod("show", signature="DoppelGang",
          function(object){
              cat ("S4 object of class:" , class ( object ) , "\n")
              cat ("Number of potential doppelgangers:" , nrow(object@summaryresults), ": ",
                   sum(object@summaryresults$expr.doppel), "expression, ",
                   sum(object@summaryresults$pheno.doppel), "phenotype, ",
                   sum(object@summaryresults$smokinggun.doppel), "smoking gun. \n")
              cat (" \n ")
              cat ("See object@summaryresults for a data.frame of doppelgangrs.", "\n")
              cat (" \n ")
          })

setMethod("plot", signature(x="DoppelGang"),
          function(x, skip.no.doppels=FALSE, plot.pair=NULL, ...){
              results <- x@fullresults
              if(!is.null(plot.pair)){
                  if(is.character(plot.pair) & length(plot.pair) == 2){
                      study.names <- unique(unlist(strsplit(names(results1@fullresults), split=x@inputargs$separator)))
                      if(!all(plot.pair %in% study.names))
                          stop("One or both of plot.pair do not match names(esets)")
                      iplot <- na.omit(unique(match(c(paste(plot.pair, collapse=":"), paste(rev(plot.pair), collapse=x@inputargs$separator)), names(results))))
                  }else{
                      stop("plot.pair must be a character vector of length two, containing the names of two datasets seen in names(esets)")
                  }
              }else{
                  iplot <- 1:length(results)
              }
              for (i in iplot){
                  cors <- results[[i]]$correlations
                  cors <- na.omit(as.numeric(cors))
                  expr.doppels <- results[[i]]$expr.doppels
                  expr.doppels <- expr.doppels[expr.doppels$doppel, ]
                  if(skip.no.doppels & (nrow(expr.doppels) == 0 | is.null(expr.doppels)))
                      next
                  hist(cors, main = names(results)[i],
                       xlab = "Pairwise Correlations", breaks="FD", ...)
                  abline(v=expr.doppels$similarity, col="red", lw=0.5)
              }})
