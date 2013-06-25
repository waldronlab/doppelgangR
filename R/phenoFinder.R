phenoFinder <- function # Calculate pairwise similarities of phenoData between samples for a list containing two ExpressionSets
### This function acts as a wrapper to phenoDist to handle cases of
### one ExpressionSet, a list of two identical ExpressionSets, or a
### list of two different ExpressionSets.
(eset.pair,
### input: a list of ExpressionSets with two elements, or an
### ExpressionSet.  If the two elements are identical, return the
### correlation matrix for pairs of samples in the first element.  If
### not identical, return pairs between the two elements.
separator=":",
### a separator between dataset name (taken from the list names) and
### sample name (taken from sampleNames(eset), to keep track of which
### samples come from which dataset.
...
### Extra arguments passed on to phenoDist
){
    if(class(eset.pair) == "ExpressionSet"){
        eset.pair <- list(eset.pair, eset.pair)
        names(eset.pair) <- 1:2
    }
    if((!identical(class(eset.pair), "list") | length(eset.pair) != 2))
        stop("eset.pair should be a list of length 2")
    if(!identical(colnames(pData(eset.pair[[1]])), colnames(pData(eset.pair[[2]]))))
        stop("Both ExpressionSets should have the same colnames in their phenoData slots.")
    matrix.one <- as.matrix(pData(eset.pair[[1]]))
    ##This part removes columns that are all NA.  The complication
    ##with keep.col is necessary because Surv objects get turned
    ##into two elements in keep.col.
    keep.col <- apply(matrix.one, 2, function(x) !all(is.na(x)))
    keep.col <- keep.col[names(keep.col) %in% colnames(matrix.one)]
    matrix.one <- matrix.one[, match(names(keep.col), colnames(matrix.one))]
    matrix.one <- matrix.one[, keep.col]
    rownames(matrix.one) <- paste(names(eset.pair)[1], rownames(matrix.one), sep=separator)
    if( identical(eset.pair[[1]], eset.pair[[2]]) ){
        ##Calculate similarity matrix for a single ExpressionSet:
        similarity.mat <- 1 - phenoDist(matrix.one, ...)
        similarity.mat[!upper.tri(similarity.mat)] <- NA  ##NA for all but upper triangle.
    }else{
        ##Calculate similarity matrix for two distinct ExpressionSets:
        matrix.two <- as.matrix(pData(eset.pair[[2]]))
        matrix.two <- matrix.two[, match(colnames(matrix.one), colnames(matrix.two))]
        rownames(matrix.two) <- paste(names(eset.pair)[2], rownames(matrix.two), sep=separator)
        similarity.mat <- 1 - phenoDist(matrix.one, matrix.two, ...)
    }
    return(similarity.mat)
### A matrix of similarities between the phenotypes of pairs of samples.
}
