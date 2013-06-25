setClass(Class = "DoppelGang", representation(fullresults = "list", summaryresults="data.frame"))

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
              if(!is.null(plot.pair) & is.character(plot.pair) &
                 length(plot.pair) == 2){
                  iplot <- match(plot.pair[1], names(results))
                  jplot <- match(plot.pair[2], names(results[[iplot]]))
                 }else{
                     iplot <- 1:length(results)
                     jplot <- 1:length(results[[1]])
                 }
                 for (i in 1:length(results)){
                     for (j in 1:length(results[[i]])){
                         if(!(i %in% iplot ) | !(j %in% jplot))
                             next
                         cors <- results[[i]][[j]]$correlations
                         cors <- na.omit(as.numeric(cors))
                         expr.doppels <- results[[i]][[j]]$expr.doppels
                         expr.doppels <- expr.doppels[expr.doppels$doppel, ]
                         if(skip.no.doppels & (nrow(expr.doppels) == 0 | is.null(expr.doppels)))
                             next
                         hist(cors, main = paste(names(results)[i],
                                    names(results[[i]])[j], sep = " / "),
                              xlab = "Pairwise Correlations", breaks="FD", ...)
                         abline(v=expr.doppels$similarity, col="red", lw=0.5)
                     }}})
