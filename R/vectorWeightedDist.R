#' Calculate a weighted distance between two vectors, using pairwise complete
#' observations.
#' 
#' Simple function to count the fraction of different elements (in the same
#' position) between two vectors of the same length, after removing elements
#' from both vectors corresponding to positions that are NA in either vector.
#' Distance is the probability for observing the matches and mismatches in two
#' random patients.
#' 
#' 
#' @param x a matrix
#' @param y a matrix with the same number of columns as x
#' @param k row in x to test for differences
#' @param l row in y to test for differences
#' @return Returns a numeric value, the log of the probability of observing the
#' matches in x and y
#' @author Levi Waldron, Markus Riester, Marcel Ramos
#' @examples
#' 
#' mymat1 <- matrix(rnorm(20), ncol = 5)
#' mymat1[1, 4] <- NA
#' mymat2 <- matrix(rnorm(20), ncol = 5)
#' vectorWeightedDist(mymat1, mymat2, 1, 2)
#' 
#' @export vectorWeightedDist
vectorWeightedDist <-
  function #Calculate a weighted distance between two vectors, using pairwise complete observations.
### Simple function to count the fraction of different elements (in
### the same position) between two vectors of the same length, after
### removing elements from both vectors corresponding to positions
### that are NA in either vector. Distance is the probability for observing
### the matches and mismatches in two random patients.
(x,
 ### a matrix
 y,
 ### a matrix with the same number of columns as x
 k,
 ### row in x to test for differences
 l,
 ### row in y to test for differences
 col_freqs_x = NULL,
 ### optional precomputed column value proportions for x
 col_freqs_y = NULL
 ### optional precomputed column value proportions for y
) {
  idx <- !(is.na(x[k, ]) | is.na(y[l, ]))
  
  if (sum(idx) < 2)
    return(1)
  
  orig_indices <- which(idx)
  
  if (!is.null(col_freqs_x) && !is.null(col_freqs_y)) {
    p.x <- vapply(orig_indices, function(i) {
      as.numeric(col_freqs_x[[i]][x[k, i]])
    }, FUN.VALUE = numeric(1))
    p.y <- vapply(orig_indices, function(i) {
      as.numeric(col_freqs_y[[i]][y[l, i]])
    }, FUN.VALUE = numeric(1))
    
    idx_diff <- x[k, orig_indices] != y[l, orig_indices]
  } else {
    x_sub <- x[, idx, drop = FALSE]
    y_sub <- y[, idx, drop = FALSE]
    
    p.x <- vapply(seq_len(ncol(x_sub)), function(i)
      sum(x_sub[k, i] == x_sub[, i],
          na.rm = TRUE) / sum(!is.na(x_sub[, i])),
      FUN.VALUE = numeric(1))
    p.y <- vapply(seq_len(ncol(y_sub)), function(i)
      sum(y_sub[l, i] == y_sub[, i],
          na.rm = TRUE) / sum(!is.na(y_sub[, i])),
      FUN.VALUE = numeric(1))
      
    idx_diff <- x_sub[k, ] != y_sub[l, ]
  }
  
  w <- 1 - p.x * p.y
  if (sum(w) < 2)
    return(1)
  1 - (sum(ifelse(idx_diff, -1, 1) * w) / sum(w) + 1) / 2
  ### Returns a numeric value, the log of the probability of observing the
  ### matches in x and y
}
