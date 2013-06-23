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
    if((class(eset.pair) != "ExpressionSet")
       & (class(eset.pair) != "list" || length(eset.pair) > 2))
        stop("eset.pair should be a list of two esets")
    if( identical(class(eset.pair), "list") & !identical(eset.pair[[1]], eset.pair[[2]]) ){
        if(!identical(colnames(pData(eset.pair[[1]])), colnames(pData(eset.pair[[2]]))))
        if(!identical(colnames(pData(eset.pair[[1]])),
        colnames(pData(eset.pair[[2]]))))
            stop("Both ExpressionSets should have the same columnames in their phenoData slots.")
        matrix.one <- as.matrix(pData(eset.pair[[1]]))
        matrix.two <- as.matrix(pData(eset.pair[[2]]))
        matrix.one <- matrix.one[, apply(matrix.one, 2, function(x) !all(is.na(x)))]
        matrix.two <- matrix.two[, match(colnames(matrix.one), colnames(matrix.two))]
        rownames(matrix.one) <- paste(names(eset.pair)[1], rownames(matrix.one), sep=separator)
        rownames(matrix.two) <- paste(names(eset.pair)[2], rownames(matrix.two), sep=separator)
        similarity.mat <- 1 - phenoDist(matrix.one, matrix.two, ...)
    }else{
        ##Calculate similarity matrix for a single ExpressionSet:
        if(identical(class(eset.pair), "list")){
            matrix.one <- pData(eset.pair[[1]])
            rownames(matrix.one) <- paste(names(eset.pair)[1], rownames(matrix.one), sep=separator)
        }else{
            matrix.one <- pData(eset.pair)
        }
        ##get rid of columns that are all NA:
        matrix.one <- matrix.one[, apply(matrix.one, 2, function(x) !all(is.na(x)))]
        similarity.mat <- 1 - phenoDist(matrix.one, ...)
        similarity.mat[!upper.tri(similarity.mat)] <- NA  ##NA for all but upper triangle.
    }
    return(similarity.mat)
### A matrix of similarities between the phenotypes of pairs of samples.
}
