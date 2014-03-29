doppelgangR <- structure(function
### Identify samples with suspiciously high correlations and phenotype similarities
(esets,
### a list of ExpressionSets, containing the numeric and phenotypic data to be analyzed.
separator=":",
### a delimitor to use between dataset names and sample names
corFinder.args=list(separator=separator, use.ComBat=TRUE, method="pearson"),
### a list of arguments to be passed to the corFinder function.
phenoFinder.args=list(separator=separator),
### a list of arguments to be passed to the phenoFinder function.  If
### NULL, samples with similar phenotypes will not be searched for.
phenoDist.args=list(vectorDistFun=vectorHammingDist),
### a list of arguments to be passed to the phenoDist function
outlierFinder.expr.args=list(bonf.pvalue=0.005, transFun=atanh, tail="upper"),
### a list of arguments to be passed to outlierFinder when called for expression data
outlierFinder.pheno.args=list(normal.upper.thresh=0.99, bonf.pvalue=NULL, tail="upper"),
### a list of arguments to be passed to outlierFinder when called for phenotype data
smokingGunFinder.args=list(transFun=I),
### a list of arguments to be passed to smokingGunFinder
impute.knn.args=list(k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069),
### a list of arguments to be passed to impute::impute.knn.  Set to
### NULL to do no knn imputation.
manual.smokingguns=NULL,
### a character vector of phenoData columns that, if identical, will
### be considered evidence of duplication
automatic.smokingguns=TRUE,
### automatically look for "smoking guns."  If TRUE, look for
### phenotype variables that are unique to each patient in dataset 1,
### also unique to each patient in dataset 2, but contain exact
### matches between datasets 1 and 2.
within.datasets.only=FALSE,
### If TRUE, only search within each dataset for doppelgangers.
cache.dir="cache"
### The name of a directory in which to cache or look up results to save
### re-calculating correlations.  Set to NULL for no caching.
){
    if(!is.null(cache.dir))
        dir.create(cache.dir, showWarnings=FALSE)
    if (is.null(names(esets)))
        names(esets) <- make.names(1:length(esets))
    if(length(esets) > length(unique(names(esets))))
        names(esets) <- make.unique(names(esets))
    for (i in 1:length(esets)){
        if(!is.null(impute.knn.args) & any(!complete.cases(exprs(esets[[i]])))){
            ##KNN imputation
            message(paste("KNN imputation for", names(esets)[i]))
            impute.knn.args$data <- exprs(esets[[i]])
            impute.knn.output  <- do.call(impute::impute.knn, args=impute.knn.args)
            exprs(esets[[i]]) <- impute.knn.output$data
            .Random.seed <- impute.knn.output$rng.state  ##restore original RNG state
        }
    }
    output.full <- lapply(1:length(esets), function(i){
        if(within.datasets.only){
            jseq <- i
        }else{
            jseq <- i:length(esets)
        }
        output2 <- lapply(jseq, function(j, i){
            message(paste("Working on datasets", names(esets)[i], "and", names(esets)[j]))
            output3 <- list()
            ## calculate correlation matrix
            corFinder.args$eset.pair <- esets[c(i, j)]
            if(!is.null(cache.dir)){
                cache.file <- paste(cache.dir, "/", digest(corFinder.args), ".rda", sep="")
                if(file.exists(cache.file))
                    load(cache.file)
            }
            if(!exists("cor.sim"))
                cor.sim <- do.call(corFinder, corFinder.args)
            if(!is.null(cache.dir) && !file.exists(cache.file))
                save(cor.sim, file=cache.file)
            output3[["correlations"]] <- cor.sim
            ## find numeric (expression) doppelgangers
            outlierFinder.expr.args$similarity.mat <- cor.sim
            outlierFinder.expr.args$prune.output <- FALSE
            output3[["expr.doppels"]] <-
                do.call(outlierFinder, outlierFinder.expr.args)
            ## If there is no phenoData in one of the esets, do not
            ## search for phenotype doppelgangers:
            if(min(c(dim(pData(esets[[1]])), dim(pData(esets[[2]])))) == 0)
                phenoFinder.args <- NULL
            ## find potential "smoking gun" phenotypes
            if(automatic.smokingguns & !is.null(phenoFinder.args)){
                new.smokinggun.phenotypes <- unlist(sapply(colnames(pData(esets[[i]])), function(cname){
                    if(cname %in% colnames(pData(esets[[j]]))){
                        if((sum(!is.na(pData(esets[[i]])[, cname])) > 2 &
                            sum(!is.na(pData(esets[[j]])[, cname])) > 2) &
                           (identical(length(pData(esets[[i]])[, cname]), length(unique(pData(esets[[i]])[, cname]))) |
                            identical(length(pData(esets[[j]])[, cname]), length(unique(pData(esets[[j]])[, cname]))))
                           )
                            return(cname)}}))
                manual.smokingguns <- unique(c(manual.smokingguns, new.smokinggun.phenotypes))
            }
            output3[["smokingguns"]] <- manual.smokingguns
            if(!is.null(manual.smokingguns)){
                smokingGunFinder.args$eset.pair <- esets[c(i, j)]
                smokingGunFinder.args$smokingguns <- manual.smokingguns
                output3[["smokinggun.doppels"]] <-
                    do.call(smokingGunFinder, smokingGunFinder.args)
            }else{
                output3[["smokinggun.doppels"]] <- NULL
            }
            ## calculate phenotype similarity matrix
            if(!is.null(phenoFinder.args)){
                phenoFinder.args$eset.pair <- esets[c(i, j)]
                ## do not use "smoking guns" when calculating phenotype distances
                for (k in 1:2)
                    if(any(manual.smokingguns %in% colnames(pData(phenoFinder.args$eset.pair[[k]]))))
                        pData(phenoFinder.args$eset.pair[[k]]) <-
                            pData(phenoFinder.args$eset.pair[[k]])[, -na.omit(match(manual.smokingguns, colnames(pData(phenoFinder.args$eset.pair[[k]]))))]
                pheno.sim <- do.call(phenoFinder, phenoFinder.args)
                ## find phenotype doppelgangers
                outlierFinder.pheno.args$similarity.mat <- pheno.sim
                outlierFinder.pheno.args$prune.output <- FALSE
                output3[["pheno.doppels"]] <- do.call(outlierFinder, outlierFinder.pheno.args)
                if(nrow(output3[["pheno.doppels"]]) == 0){
                    output3[["pheno.doppels"]] <- output3[["expr.doppels"]][, 1:2]
                    output3[["pheno.doppels"]]$similarity <- NA
                    output3[["pheno.doppels"]]$doppel <- FALSE
                }
                ## This next block seems like a complicated way just to
                ## prune sample pairs that are not identified by any of
                ## the three searches.  Problem is that pheno.doppels and
                ## expr.doppels keep all pairs, but smokinggun.doppels
                ## shows only suspected doppels.
                keep.rows <- output3[["pheno.doppels"]]$doppel | output3[["expr.doppels"]]$doppel
                if(!is.null(output3[["smokinggun.doppels"]])){
                    sample1.pheno <- output3[["pheno.doppels"]][, 1]
                    sample2.pheno <- output3[["pheno.doppels"]][, 2]
                    sample1.smokinggun <- output3[["smokinggun.doppels"]][, 1]
                    sample2.smokinggun <- output3[["smokinggun.doppels"]][, 2]
                    keep.rows <- keep.rows | (paste(sample1.pheno, sample2.pheno, collapse=":::") %in% paste(sample1.smokinggun, sample2.smokinggun, collapse=":::"))
                }
                output3[["pheno.doppels"]] <- output3[["pheno.doppels"]][keep.rows, ]
                output3[["expr.doppels"]] <- output3[["expr.doppels"]][keep.rows, ]
            }else{
                ##Make pheno.doppels a zero-line data.frame with the same colnames as expr.doppels:
                output3[["pheno.doppels"]] <- output3[["expr.doppels"]][0, ]
                output3[["expr.doppels"]] <- output3[["expr.doppels"]][output3[["expr.doppels"]]$doppel, ]
            }
            return(output3)
        }, i=i)
        names(output2) <- names(esets)[jseq]
        return(output2)
    })
    names(output.full) <- names(esets)
    message("Finalizing...")
    wrapUp <- function(object, element){
        tmp <- lapply(object, function(x) lapply(x, function(y) y[[element]]))
        tmp <- lapply(object, function(x) do.call(rbind, lapply(x, function(y) y[[element]])))
        tmp <- tmp[!sapply(tmp, is.null)]
        do.call(rbind, tmp)
    }
    pheno.doppels <- wrapUp(output.full, "pheno.doppels")
    expr.doppels <- wrapUp(output.full, "expr.doppels")
    smokinggun.doppels <- wrapUp(output.full, "smokinggun.doppels")
    ## merge all doppelganger types
    all.doppels <- expr.doppels
    colnames(all.doppels)[3:4] <- paste("expr.", colnames(all.doppels)[3:4], sep="")
    addCols <- function(orig, add){
        newrows1 <- paste(add[, 1], add[, 2])
        newrows2 <- paste(add[, 2], add[, 1])
        oldrows <- paste(orig[, 2], orig[, 1])
        if(any(!(newrows1 %in% oldrows | newrows2 %in% oldrows))){
            tmp <- add[!(newrows1 %in% oldrows | newrows2 %in% oldrows), 1:2]
            tmp <- cbind(tmp, matrix(NA, nrow=nrow(tmp), ncol=(ncol(orig) - 2)))
            colnames(tmp) <- colnames(orig)
            orig <- rbind(orig, tmp)
        }
        match.rows <- match(paste(add[, 1], add[, 2]), paste(orig[, 1], orig[, 2]))
        if(nrow(add) > 0 & nrow(orig) == 0){
            ##adding to a zero-row dataframe
            newdf <- add[, 1:2]
            for (i in 3:ncol(orig))
                newdf[[colnames(orig)[i]]] <- NA
            orig <- newdf
        }
        if(nrow(add) == 0 & nrow(orig) == 0){
            return(cbind(orig, add[, -1:-2]))
        }
        orig[[colnames(add)[3]]] <- NA
        orig[[colnames(add)[3]]][match.rows] <- add[, 3]
        orig[[colnames(add)[4]]] <- FALSE
        orig[[colnames(add)[4]]][match.rows] <- add[, 4]
        orig
    }
    all.doppels <- addCols(all.doppels, pheno.doppels)
    colnames(all.doppels)[5:6] <- paste("pheno.", colnames(all.doppels)[5:6], sep="")
    if(is.null(smokinggun.doppels)){
        if(nrow(all.doppels) > 0){
            all.doppels$smokinggun.similarity <- 0
            all.doppels$smokinggun.doppel <- FALSE
        }else{
            all.doppels <- cbind(all.doppels, all.doppels[, 3:4])
            colnames(all.doppels)[7:8] <- c("smokinggun.similarity", "smokinggun.doppel")
        }
    }else{
        all.doppels <- addCols(all.doppels, smokinggun.doppels)
        colnames(all.doppels)[7:8] <- paste("smokinggun.", colnames(all.doppels)[7:8], sep="")
    }
    all.doppels <- all.doppels[all.doppels[, 4] | all.doppels[, 6] | all.doppels[, 8], ]
    rownames(all.doppels) <- NULL
    ##If all esets have the same colnames of pdata, add merged
    ##pairwise sample pdata to all.doppels:
    if( identical(nrow(all.doppels) > 0, TRUE) &&
       identical(colnames(pData(esets[[1]])), unique(unlist(lapply(esets, function(eset) colnames(pData(eset)))))) ){
        merged.pdat <- apply(as.matrix(all.doppels[, 1:2]), 1, function(x){
            pdat1 <- pData(esets[[strsplit(x[1], split=separator)[[1]][1]]])[strsplit(x[1], split=separator)[[1]][2], ]
            pdat2 <- pData(esets[[strsplit(x[2], split=separator)[[1]][1]]])[strsplit(x[2], split=separator)[[1]][2], ]
            paste(pdat1, pdat2, sep=separator)
        })
        merged.pdat <- t(merged.pdat)
        colnames(merged.pdat) <- colnames(pData(esets[[1]]))
        all.doppels <- cbind(all.doppels, merged.pdat)
    }
    new("DoppelGang", fullresults=output.full, summaryresults=all.doppels)
### Returns an object of S4-class "DoppelGang".  See ?DoppelGang-class.
}, ex=function(){
    library(doppelgangR)
    library(curatedOvarianData)

    data(GSE32062.GPL6480_eset)
    data(GSE32063_eset)
    data(GSE12470_eset)
    data(GSE17260_eset)

    testesets <- list(JapaneseA=GSE32062.GPL6480_eset,
                      JapaneseB=GSE32063_eset,
                      Yoshihara2009=GSE12470_eset,
                      Yoshihara2010=GSE17260_eset)

    testesets <- lapply(testesets, function(X){
        # standardize the sample ids to improve matching based on clinical annotation
        sampleNames(X) <- make.names(paste(X$sample_type,
            gsub("\\\\D","",X$alt_sample_name), sep="_"))
        X$alt_sample_name <- sampleNames(X)
        pData(X) <- pData(X)[, !grepl("uncurated_author_metadata", colnames(pData(X)))]
        X })

    results1 <- doppelgangR(testesets, corFinder.args=list(use.ComBat=TRUE), cache.dir=NULL)
    summary(results1)
    plot(results1)
    ## Set phenoFinder.args to ignore similar phenotypes:
    results2 <- doppelgangR(testesets, corFinder.args=list(use.ComBat=FALSE), phenoFinder.args=NULL, cache.dir=NULL)
    summary(results2)
})
