
phenoDist <- function ###Calculate distance between two vectors, rows of one matrix/dataframe, or rows of two matrices/dataframes.
(x,
### A vector, matrix or dataframe
y=NULL,
### NULL, a vector, matrix, or dataframe.  If x is a vector, y must also be specified.
vectorDistFun=vectorHammingDist,
### A function of two vectors that returns the distance between those vectors.
...
### Extra arguments passed on to vectorDistFun
){
    if (is.vector(x) && is.vector(y)) {
        z <- vectorDistFun(x, y, ...)
    }
    else {
        if(identical(class(x), "data.frame"))
            x <- as.matrix(x)
        if(identical(class(y), "data.frame"))
            y <- as.matrix(y)
        if(is.null(y)){
            z <- matrix(0, nrow = nrow(x), ncol = nrow(x))
            for (k in 1:(nrow(x) - 1)) {
                for (l in (k + 1):nrow(x)) {
                    z[k, l] <- vectorDistFun(x[k, ], x[l, ], ...)
                    z[l, k] <- z[k, l]
                }
            }
            dimnames(z) <- list(rownames(x), rownames(x))
        }else{
            z <- matrix(0, nrow = nrow(x), ncol = nrow(y))
            for (k in 1:(nrow(x))) {
                for (l in 1:nrow(y)) {
                    z[k, l] <- vectorDistFun(x[k, ], y[l, ], ...)
                }
            }
            dimnames(z) <- list(rownames(x), rownames(y))
        }
    }
    z
### a matrix of distances between pairs of rows of x (if y is
### unspecified), or between all pairs of rows between x and y (if
### both are provided).
}
