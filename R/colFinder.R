.checkSameSets <- function(sets){
  matrix.one <- assay(sets[[1]])
  dimnames(matrix.one) <- NULL
  matrix.two <- assay(sets[[2]])
  dimnames(matrix.two) <- NULL
  coldata.one <- colData(sets[[1]])
  rownames(coldata.one) <- NULL
  coldata.two <- colData(sets[[2]])
  rownames(coldata.two) <- NULL
  output <- 
    identical(matrix.one, matrix.two) &&
    identical(coldata.one, coldata.two)
  return(output)
}

#' Calculate pairwise similarities of colData between samples for a list
#' containing two DataFrame
#' 
#' This function acts as a wrapper to colData to handle cases of one
#' DataFrame, a list of two identical DataFrame, or a list of two
#' different DataFrame
#' 
#' 
#' @param summex.list input: a list of DataFrame with two elements, or a
#' DataFrame. If the two elements are identical, return the correlation
#' matrix for pairs of samples in the first element.  If not identical, return
#' pairs between the two elements.
#' @param \dots Extra arguments passed on to colFinder
#' @return A matrix of similarities between the colData of pairs of samples.
#' @author Fabio Da Col, Marcel Ramos
#' 
#' @export colFinder
colFinder <- function(summex.list, ...) {
    if (!is(summex.list, "list") | length(summex.list) != 2)
        stop("list should be a list of length 2")

  if (!identical(colnames(colData(summex.list[[1]])), colnames(colData(summex.list[[2]])))) {
    stop(
      "Slots of list must have identical column names, e.g.
      identical(colnames(colData(summex.list[[1]])), colnames(colData(summex.list[[2]])))"
    )
  }
  matrix.one <- as.matrix(colData(summex.list[[1]]))
  if(is.null(rownames(matrix.one)))
    rownames(matrix.one) <- make.names(1:nrow(matrix.one))
  keep.col <- apply(matrix.one, 2, function(x){
    !all(is.na(x))
  })
  keep.col <- keep.col[names(keep.col) %in% colnames(matrix.one)]
  matrix.one <- 
    subset(matrix.one, select = match(names(keep.col), colnames(matrix.one)))
  matrix.one <- subset(matrix.one, select = keep.col)
  if(.checkSameSets(summex.list)){
    similarity.mat <- 1 - phenoDist(matrix.one, ...)
    similarity.mat[!upper.tri(similarity.mat)] <- NA
  }else{
    matrix.two <- as.matrix(summex.list[[2]])
    matrix.two <-
      matrix.two[, match(colnames(matrix.one), colnames(matrix.two)), drop = FALSE]
    if(is.null(rownames(matrix.two)))
      rownames(matrix.two) <- make.names(1:nrow(matrix.two))
    similarity.mat <- 1 - phenoDist(matrix.one, matrix.two, ...)
  }
  return(similarity.mat)
}