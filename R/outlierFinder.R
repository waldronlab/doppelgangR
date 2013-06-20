outlierFinder <- function ###Identifies outliers in a similarity matrix. 
### By default uses the
### Fisher z-transform for Pearson correlation (atanh), and
### identifies outliers as those above the quantile of a normal
### distribution with mean and standard deviation estimated from the
### z-transformed matrix.  The quantile is calculated from the
### Bonferroni-corrected cumulative probability of the upper tail.
(similarity.mat,
### A matrix of similarities - larger values mean more similar.
bonf.pvalue=0.05,
### Bonferroni-corrected p-value.  A raw.pvalue is calculated by
### dividing this by the number of non-missing values in
### similarity.mat, and the rejection threshold is qnorm(1-raw.pvalue,
### mean, sd) where mean and sd are estimated from the
### transFun-transformed similarity.mat.
transFun=atanh,
### A function applied to the numeric values of similarity.mat, that
### should result in normally-distributed values.
normal.upper.thresh=NULL,
### Instead of specifying bonf.pvalue and transFun, an upper
### similarity threshold can be set, and values above this will be
### considered likely duplicates.
prune.output=TRUE
### If prune.output=TRUE, only return likely doppelgangers.
){
    if(!is.null(bonf.pvalue) & !is.null(normal.upper.thresh))
        stop("Specify only one of bonf.pvalue and normal.upper.thresh")
    if(is.null(normal.upper.thresh) & !is.null(bonf.pvalue) & !is.null(transFun)){
        zmat <- transFun(similarity.mat)
        raw.pvalue <- bonf.pvalue / sum(!is.na(zmat))
        z.cutoff <- qnorm(1-raw.pvalue, mean=mean(as.numeric(zmat), na.rm=TRUE), sd=sd(as.numeric(zmat), na.rm=TRUE))
        outlier.mat <- zmat > z.cutoff
    }else if(!is.null(normal.upper.thresh)){
        outlier.mat <- similarity.mat > normal.upper.thresh
    }else{
        return(NULL)
    }
    outlier.mat[is.na(outlier.mat)] <- FALSE
##    output <- data.frame(sample1=rownames(outlier.mat)[row(outlier.mat)[which(outlier.mat)]],
##                         sample2=colnames(outlier.mat)[col(outlier.mat)[which(outlier.mat)]],
##                         similarity=similarity.mat[which(outlier.mat)])
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
