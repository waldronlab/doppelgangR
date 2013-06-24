setClass(Class = "DoppelGang", representation(results = "data.frame", 
                correlations = "list", smokingguns = "list"))

## setGeneric("DoppelGang", function(object) standardGeneric("DoppelGang"))

setGeneric("print")
setGeneric("summary")

setMethod("print", signature(x="DoppelGang"),
          function(x) print(x@results))

setMethod("summary", signature(object="DoppelGang"),
          function(object) summary(object@results))

