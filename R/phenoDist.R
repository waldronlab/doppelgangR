

phenoDist <-
  structure(
    function #Calculate distance between two vectors, rows of one matrix/dataframe, or rows of two matrices/dataframes.
    ### This function does some simple looping to allow x and y to be
    ### various combinations of vectors and matrices/dataframes.
    (
      x,
      ### A vector, matrix or dataframe
      y = NULL,
      ### NULL, a vector, matrix, or dataframe.  If x is a vector, y must also be specified.
      bins = 10,
      ### discretize continuous fields in the specified number of bins
      vectorDistFun = vectorWeightedDist,
      ### A function of two vectors that returns the distance between those vectors.
      ...
      ### Extra arguments passed on to vectorDistFun
    ) {
      if (is.vector(x) && is.vector(y)) {
        z <- vectorDistFun(matrix(x, nrow = 1), matrix(y, nrow = 1), 1, 1, ...)
      }
      else {
        x <- .discretizeDataFrame(x, bins)
        if (is.null(y)) {
          z <- matrix(0, nrow = nrow(x), ncol = nrow(x))
          for (k in 1:(nrow(x) - 1)) {
            for (l in (k + 1):nrow(x)) {
              z[k, l] <- vectorDistFun(x, x, k, l, ...)
              z[l, k] <- z[k, l]
            }
          }
          dimnames(z) <- list(rownames(x), rownames(x))
        } else{
          y <- .discretizeDataFrame(y, bins)
          z <- matrix(0, nrow = nrow(x), ncol = nrow(y))
          for (k in 1:(nrow(x))) {
            for (l in 1:nrow(y)) {
              z[k, l] <- vectorDistFun(x, y, k, l, ...)
            }
          }
          dimnames(z) <- list(rownames(x), rownames(y))
        }
      }
      z
      ### a matrix of distances between pairs of rows of x (if y is
      ### unspecified), or between all pairs of rows between x and y (if
      ### both are provided).
    },
    ex = function() {
      library(curatedOvarianData)
      data(GSE32063_eset)
      data(GSE17260_eset)
      pdat1 <- pData(GSE32063_eset)
      pdat2 <- pData(GSE17260_eset)
      ## Curation of the alternative sample identifiers makes duplicates stand out more:
      pdat1$alt_sample_name <-
        paste(pdat1$sample_type,
              gsub("[^0-9]", "", pdat1$alt_sample_name),
              sep = "_")
      pdat2$alt_sample_name <-
        paste(pdat2$sample_type,
              gsub("[^0-9]", "", pdat2$alt_sample_name),
              sep = "_")
      ## Removal of columns that cannot possibly match also helps duplicated patients to stand out
      pdat1 <-
        pdat1[,!grepl("uncurated_author_metadata", colnames(pdat1))]
      pdat2 <-
        pdat2[,!grepl("uncurated_author_metadata", colnames(pdat2))]
      ## Use phenoDist() to calculate a weighted distance matrix
      distmat <- phenoDist(as.matrix(pdat1), as.matrix(pdat2))
      ## Note outliers with identical clinical data, these are probably the same patients:
      graphics::boxplot(distmat)
    }
  )

.discretizeDataFrame <- function(X, bins) {
  .discretizeRow <- function(x) {
    if (length(levels(as.factor(x))) > bins)
      return(cut(x, breaks = bins))
    as.factor(x)
  }
  idx <- apply(X, 2L, is.numeric)
  if (sum(idx) == 0)
    return(X)
  X[, idx] <-
    as.data.frame(apply(X[, idx, drop = FALSE], 2, .discretizeRow))
  X
}
