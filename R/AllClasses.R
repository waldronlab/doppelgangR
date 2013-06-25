setClass(Class = "DoppelGang", representation(fullresults = "list", summaryresults="data.frame"))

## setGeneric("DoppelGang", function(object) standardGeneric("DoppelGang"))

setGeneric("print")
setGeneric("summary")

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

