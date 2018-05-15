#' Calculate pairwise similarities of phenoData between samples for a list
#' containing two ExpressionSets
#' 
#' This function acts as a wrapper to phenoDist to handle cases of one
#' ExpressionSet, a list of two identical ExpressionSets, or a list of two
#' different ExpressionSets.
#' 
#' 
#' @param eset.pair input: a list of ExpressionSets with two elements, or an
#' ExpressionSet.  If the two elements are identical, return the correlation
#' matrix for pairs of samples in the first element.  If not identical, return
#' pairs between the two elements.
#' @param separator a separator between dataset name (taken from the list
#' names) and sample name (taken from sampleNames(eset), to keep track of which
#' samples come from which dataset.
#' @param \dots Extra arguments passed on to phenoDist
#' @return A matrix of similarities between the phenotypes of pairs of samples.
#' @author Levi Waldron, Markus Riester, Marcel Ramos
#' @examples
#' 
#' library(curatedOvarianData)
#' data(GSE32063_eset)
#' data(GSE17260_eset)
#' esets2 <- list(JapaneseB=GSE32063_eset,
#'                 Yoshihara2010=GSE17260_eset)
#' 
#' ## standardize the sample ids to improve matching based on clinical annotation
#' esets2 <- lapply(esets2, function(X){
#'     X$alt_sample_name <- paste(X$sample_type, gsub("[^0-9]", "", X$alt_sample_name), sep="_")
#' 
#' ## Removal of columns that cannot possibly match also helps duplicated patients to stand out
#'     pData(X) <- pData(X)[, !grepl("uncurated_author_metadata", colnames(pData(X)))]
#'     X <- X[, 1:20]  ##speed computations
#'     return(X) })
#' 
#' ## See first six samples in both rows and columns
#' phenoFinder(esets2)[1:6, 1:6]
#' 
#' @export phenoFinder
phenoFinder <-
  function # Calculate pairwise similarities of phenoData between samples for a list containing two ExpressionSets
### This function acts as a wrapper to phenoDist to handle cases of
### one ExpressionSet, a list of two identical ExpressionSets, or a
### list of two different ExpressionSets.
(eset.pair,
 ### input: a list of ExpressionSets with two elements, or an
 ### ExpressionSet.  If the two elements are identical, return the
 ### correlation matrix for pairs of samples in the first element.  If
 ### not identical, return pairs between the two elements.
 separator = ":",
 ### a separator between dataset name (taken from the list names) and
 ### sample name (taken from sampleNames(eset), to keep track of which
 ### samples come from which dataset.
 ...
 ### Extra arguments passed on to phenoDist
) {
  if ((!is(eset.pair, "list") | length(eset.pair) != 2))
    stop("eset.pair should be a list of length 2")
  if (!identical(colnames(pData(eset.pair[[1]])), colnames(pData(eset.pair[[2]]))))
    stop(
      "pData slots of esets must have identical column names, e.g.
     identical(colnames(pData(eset[[1]])), colnames(pData(eset[[2]])))"
    )
  matrix.one <- as.matrix(pData(eset.pair[[1]]))
  if (is.null(rownames(matrix.one)))
    rownames(matrix.one) <- make.names(1:nrow(matrix.one))
  ##This part removes columns that are all NA.  The complication
  ##with keep.col is necessary because Surv objects get turned
  ##into two elements in keep.col.
  keep.col <- apply(matrix.one, 2, function(x)
    ! all(is.na(x)))
  keep.col <- keep.col[names(keep.col) %in% colnames(matrix.one)]
  matrix.one <-
    subset(matrix.one, select = match(names(keep.col), colnames(matrix.one)))
  matrix.one <- subset(matrix.one, select = keep.col)
  if (.checkSameEsets(eset.pair)) {
    ##Calculate similarity matrix for a single ExpressionSet:
    similarity.mat <- 1 - phenoDist(matrix.one, ...)
    similarity.mat[!upper.tri(similarity.mat)] <-
      NA  ##NA for all but upper triangle.
  } else{
    ##Calculate similarity matrix for two distinct ExpressionSets:
    matrix.two <- as.matrix(pData(eset.pair[[2]]))
    matrix.two <-
      matrix.two[, match(colnames(matrix.one), colnames(matrix.two)), drop =
                   FALSE]
    if (is.null(rownames(matrix.two)))
      rownames(matrix.two) <- make.names(1:nrow(matrix.two))
    similarity.mat <- 1 - phenoDist(matrix.one, matrix.two, ...)
  }
  return(similarity.mat)
  ### A matrix of similarities between the phenotypes of pairs of samples.
}
