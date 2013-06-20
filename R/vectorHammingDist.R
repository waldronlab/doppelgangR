vectorHammingDist <- function ###Calculate Hamming Distance between two vectors, using pairwise complete observations.
(x,
### a vector (normally character)
y
### a vector the same length as x
){
    z <- x != y
    z <- sum(z, na.rm=TRUE) / length(na.omit(z))
    return(z)
### Returns a numeric value, the Hamming Distance (the number of
### non-equal values between x and y).
}
