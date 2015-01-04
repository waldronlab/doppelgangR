smokingGunFinder <- function  #Find doppelgangers based on "smoking gun" phenotypes - those that should be unique to each patient.
### Checks all pairwise combinations of samples for values of the
### "smoking" gun phenotypes that are identical.
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
        separator <- ""
    }
    if(class(eset.pair) != "list" | length(eset.pair) > 2)
        stop("eset.pair should be a list of two esets")
    smokingmat <- matrix(0, nrow=ncol(eset.pair[[1]]), ncol=ncol(eset.pair[[2]]))
    rownames(smokingmat) <- paste(names(eset.pair)[1], sampleNames(eset.pair[[1]]), sep=separator)
    colnames(smokingmat) <- paste(names(eset.pair)[2], sampleNames(eset.pair[[2]]), sep=separator)
    for (x in smokingguns){
        if(!(x %in% colnames(pData(eset.pair[[1]]))) | !(x %in% colnames(pData(eset.pair[[2]]))))
            next
        pdat.vec1 <- transFun(pData(eset.pair[[1]])[, x])
        pdat.vec2 <- transFun(pData(eset.pair[[2]])[, x])
        for (i in 1:length(pdat.vec1))
            smokingmat[i, pdat.vec2 %in% pdat.vec1[i]] <- smokingmat[i, pdat.vec2 %in% pdat.vec1[i]] + 1
    }
    if(identical(pdat.vec1, pdat.vec2))
        smokingmat[!upper.tri(smokingmat) & !lower.tri(smokingmat)] <- NA
    return(smokingmat)
### Returns an adjacency matrix for samples where matches have value
### 1, non-matches have value zero.  Value for a sample against itself
### is NA.
}
