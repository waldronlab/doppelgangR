smokingGunFinder <- function  ##Find doppelgangers based on "smoking gun" phenotypes - those that should be unique to each patient.
(eset.pair,
### a list of ExpressionSets, with two elements.  If the two elements
### are identical, the function will check for duplicate IDs within
### one element. If not identical, it will check for duplicate IDs
### between elements.
smokingguns,
### phenoData column names found in multiple elements of eset.pair
### that may contain "smoking guns" such as identifiers that should be
### unique to each sample.
transFun=I,
### a function to apply to IDs before comparing.  By default apply no transformation.
separator=":"
### Separator between dataset name and sample name.  Dataset names are
### added to sample names to keep track of dataset of origin.
){
    if(class(eset.pair) == "ExpressionSet"){
        eset.pair <- list(eset.pair, eset.pair)
        separator=""
    }
    if(class(eset.pair) != "list" | length(eset.pair) > 2)
        stop("eset.pair should be a list of two esets")
    smokinggun.doppels <- sapply(smokingguns, function(x){
        if(!(x %in% colnames(pData(eset.pair[[1]]))) | !(x %in% colnames(pData(eset.pair[[2]])))){
            warning(paste(x, "not found in one of the pData(eset)"))
            return(NULL)
        }
        pdat.vec1 <- transFun(pData(eset.pair[[1]])[, x])
        pdat.vec2 <- transFun(pData(eset.pair[[2]])[, x])
        if(identical(pdat.vec1, pdat.vec2)){
            nonunique.elements <- names(table(pdat.vec1))[table(pdat.vec1) > 1]
            nonunique.samplenames <- sampleNames(eset.pair[[1]])[pdat.vec1 %in% nonunique.elements]
            if(length(nonunique.samplenames) > 0){
                nonunique.samplenames <- paste(names(eset.pair)[1], nonunique.samplenames, sep=separator)
                return(as.matrix(.outer2df(nonunique.samplenames, nonunique.samplenames, FALSE, FALSE)))
            }else{
                return(NULL)
            }
        }else{
            if(any(pdat.vec1 %in% pdat.vec2)){
                t(sapply(pdat.vec1[pdat.vec1 %in% pdat.vec2], function(y){
                    samplenames1 <- paste(names(eset.pair)[1], sampleNames(eset.pair[[1]])[pdat.vec1 %in% y], sep=separator)
                    samplenames2 <- paste(names(eset.pair)[2], sampleNames(eset.pair[[2]])[pdat.vec2 %in% y], sep=separator)
                    return(as.matrix(.outer2df(samplenames1, samplenames2, bidirectional=FALSE, diag=TRUE)))
                }))
            }else{
                return(NULL)
            }
        }})
    smokinggun.doppels <- data.frame(do.call(rbind, smokinggun.doppels), stringsAsFactors=FALSE)
    if(identical(nrow(smokinggun.doppels) > 0, TRUE)){
        colnames(smokinggun.doppels)[1:2] <- c("sample1", "sample2")
        smokinggun.doppels$identifier <- rownames(smokinggun.doppels)
        smokinggun.doppels$doppel <- TRUE
        return(smokinggun.doppels)
    }else{
        return(NULL)
    }
### Returns a four-column dataframe with columns: sample1, sample2,
### identifier, doppel.  "identifier" is the value of the supposedly
### unique ID that these samples had in common, doppel is TRUE.
}
