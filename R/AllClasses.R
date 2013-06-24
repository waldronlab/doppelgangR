setClass(Class = "DoppelGang", representation(results = "data.frame", 
                correlations = "list", smokingguns = "list"))

## setGeneric("DoppelGang", function(object) standardGeneric("DoppelGang"))

setGeneric("print")
setGeneric("summary")

setMethod("print", signature(x="DoppelGang"),
          function(x) print(x@results))

setMethod("summary", signature(object="DoppelGang"),
          function(object) summary(object@results))

setMethod("show", signature="DoppelGang",
          function(object){
              cat ("S4 object of class:" , class ( object ) , "\n")
              cat ("Number of potential doppelgangers:" , nrow(object@results), ": ",
                   sum(object@results$expr.doppel), "expression, ",
                   sum(object@results$pheno.doppel), "phenotype, ",
                   sum(object@results$smokinggun.doppel), "smoking gun. \n")
              cat (" \n ")
              cat ("See object@results for a data.frame of results.", "\n")
              cat (" \n ")
          })

