setClass(Class = "DoppelGang", representation(fullresults = "list", summaryresults="data.frame"))

## setGeneric("DoppelGang", function(object) standardGeneric("DoppelGang"))

setGeneric("print")
setGeneric("summary")
setGeneric("plot")

setMethod("print", signature(x="DoppelGang"),
          function(x) print(x@summaryresults))

setMethod("summary", signature(object="DoppelGang"),
          function(object) summary(object@summaryresults))

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
          function(x){
              for (i in 1:length(x@fullresults)){
                  for (j in 1:length(x@fullresults[[i]])){
                      cors <- x@fullresults[[i]][[j]]$correlations
                      cors <- cors[upper.tri(cors)]
                      hist(cors, main = paste(names(x@fullresults)[i],
                                              names(x@fullresults[[i]])[j], sep = " / "),
                           xlab = "Pairwise Correlations", breaks="FD")
                      expr.doppels <- x@fullresults[[i]][[j]]$expr.doppels
                      abline(v=expr.doppels$similarity, col="red")
                      par(ask=TRUE)
                  }}})
