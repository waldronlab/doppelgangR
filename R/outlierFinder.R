outlierFinder <- function ###Identifies outliers in a similarity matrix.
### By default uses the
### Fisher z-transform for Pearson correlation (atanh), and
### identifies outliers as those above the quantile of a normal
### distribution with mean and standard deviation estimated from the
### z-transformed matrix.  The quantile is calculated from the
### Bonferroni-corrected cumulative probability of the upper tail.
(similarity.mat,
### A matrix of similarities - larger values mean more similar.
bonf.prob=0.05,
### Bonferroni-corrected p-value.  A raw.prob is calculated by
### dividing this by the number of non-missing values in
### similarity.mat, and the rejection threshold is qnorm(1-raw.prob,
### mean, sd) where mean and sd are estimated from the
### transFun-transformed similarity.mat.
transFun=atanh,
### A function applied to the numeric values of similarity.mat, that
### should result in normally-distributed values.
normal.upper.thresh=NULL,
### Instead of specifying bonf.prob and transFun, an upper
### similarity threshold can be set, and values above this will be
### considered likely duplicates.
tail="upper",
### "upper" to look for samples with very high similarity values,
### "lower" to look for very low values, or "both" to look for both.
prune.output=TRUE
### If prune.output=TRUE, only return likely doppelgangers.
){
    if(!is.null(bonf.prob) & !is.null(normal.upper.thresh))
        stop("Specify only one of bonf.prob and normal.upper.thresh")
    if(is.null(normal.upper.thresh) & !is.null(bonf.prob) & !is.null(transFun)){
        zmat <- transFun(similarity.mat)
        raw.prob <- bonf.prob / sum(!is.na(zmat))
        if(identical(tail, "upper")){
            z.cutoff <- qnorm(1-raw.prob, mean=mean(as.numeric(zmat), na.rm=TRUE), sd=sd(as.numeric(zmat), na.rm=TRUE))
            outlier.mat <- zmat > z.cutoff
        }else if(identical(tail, "lower")){
            z.cutoff <- qnorm(raw.prob, mean=mean(as.numeric(zmat), na.rm=TRUE), sd=sd(as.numeric(zmat), na.rm=TRUE))
            outlier.mat <- zmat < z.cutoff
        }else if(identical(tail, "both")){
            z.cutoff <- qnorm(c(raw.prob, 1-raw.prob), mean=mean(as.numeric(zmat), na.rm=TRUE), sd=sd(as.numeric(zmat), na.rm=TRUE))
            outlier.mat <- (zmat < z.cutoff[1]) | (zmat > z.cutoff[2])
        }else{ stop("tail argument should be upper, lower, or both.") }
    }else if(!is.null(normal.upper.thresh)){
        outlier.mat <- similarity.mat > normal.upper.thresh
    }else{
        return(NULL)
    }
    outlier.mat[is.na(outlier.mat)] <- FALSE
    output <- .outer2df(rownames(outlier.mat), colnames(outlier.mat), bidirectional=TRUE, diag=TRUE)
    output$similarity <- .outer2df(similarity.mat, bidirectional=TRUE, diag=TRUE)
    output$doppel <- .outer2df(outlier.mat, bidirectional=TRUE, diag=TRUE)
    output <- output[!is.na(output$similarity), ]
    output <- output[output[, 1] != output[, 2], ]
    if(prune.output){
        if (!any(outlier.mat))
            return(NULL)
        output <- output[output$doppel, ]
    }
    colnames(output)[1:2] <- c("sample1", "sample2")
    return(output)
### Returns either NULL or a dataframe with three columns: sample1, sample2, and similarity.
}
