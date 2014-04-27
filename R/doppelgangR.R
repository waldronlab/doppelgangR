doppelgangR <- structure(function
### Identify samples with suspiciously high correlations and phenotype similarities
(esets,
### a list of ExpressionSets, containing the numeric and phenotypic data to be analyzed.
separator=":",
### a delimitor to use between dataset names and sample names
corFinder.args=list(separator=separator, use.ComBat=TRUE, method="pearson"),
### a list of arguments to be passed to the corFinder function.
phenoFinder.args=list(separator=separator, vectorDistFun=vectorWeightedDist),
### a list of arguments to be passed to the phenoFinder function.  If
### NULL, samples with similar phenotypes will not be searched for.
outlierFinder.expr.args=list(bonf.prob=0.5, transFun=atanh, tail="upper"),
### a list of arguments to be passed to outlierFinder when called for expression data
outlierFinder.pheno.args=list(normal.upper.thresh=0.99, bonf.prob=NULL, tail="upper"),
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
cache.dir="cache",
### The name of a directory in which to cache or look up results to save
### re-calculating correlations.  Set to NULL for no caching.
verbose=TRUE
### Print progress information
 ){
    ##Save input args except for esets:
    input.argnames <- ls()[-match("esets", ls())]
    input.args <- lapply(input.argnames, function(x) get(x))
    names(input.args) <- input.argnames
    input.args$esets.names <- names(esets)
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
    ds.combns <- lapply(1:length(esets), function(i) c(i, i))
    if(!within.datasets.only)
        ds.combns <- c(ds.combns, combn(1:length(esets), 2, simplify=FALSE))

    output.full <- lapply(ds.combns, function(ij){
        i <- ij[1]
        j <- ij[2]
        if (verbose) message(paste("Working on datasets", names(esets)[i], "and", names(esets)[j]))
        output3 <- list()
        ## calculate correlation matrix
        corFinder.args$eset.pair <- esets[c(i, j)]
        if(!is.null(cache.dir)){
            cache.file <- paste(cache.dir, "/", digest::digest(list(corFinder, corFinder.args)), ".rda", sep="")
            if(file.exists(cache.file)) {
		if (verbose) message("\tSkipping corFinder, loading cached results.")
                load(cache.file)
            }
        }
        if(!exists("cor.sim")){
            if(verbose) message("Calculating correlations...")
            cor.sim <- do.call(corFinder, corFinder.args)
        }
        if(!is.null(cache.dir) && !file.exists(cache.file))
            save(cor.sim, file=cache.file)
        output3[["correlations"]] <- cor.sim
        ## find numeric (expression) doppelgangers
        outlierFinder.expr.args$similarity.mat <- cor.sim
        outlierFinder.expr.args$prune.output <- FALSE
        if(verbose) message("Identifying correlation doppelgangers...")
        output3[["expr.doppels"]] <-
            do.call(outlierFinder, outlierFinder.expr.args)
        ## If there is no phenoData in one of the esets, do not
        ## search for phenotype doppelgangers:
        if(min(c(dim(pData(esets[[1]])), dim(pData(esets[[2]])))) == 0)
            phenoFinder.args <- NULL
        ## automatically find potential "smoking gun" phenotypes
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
        output3[["smokingguns"]] <- manual.smokingguns  ##FIXME: write full arguments somewhere else
        if(!is.null(manual.smokingguns)){
            ## find smokinggun doppelgangers
            if(verbose) message("Identifying smoking-gun doppelgangers...")
            smokingGunFinder.args$eset.pair <- esets[c(i, j)]
            smokingGunFinder.args$smokingguns <- manual.smokingguns
            outlierFinder.smokinggun.args <- list()
            outlierFinder.smokinggun.args$similarity.mat <- do.call(smokingGunFinder, smokingGunFinder.args)
            outlierFinder.smokinggun.args$prune.output <- FALSE
            outlierFinder.smokinggun.args$normal.upper.thresh <- 0.5
            output3[["smokinggun.doppels"]] <- do.call(outlierFinder, outlierFinder.smokinggun.args)
        }
        ## calculate phenotype similarity matrix
        if(!is.null(phenoFinder.args)){
            phenoFinder.args$eset.pair <- esets[c(i, j)]
            ##keep no rows because we only need pData(x):
            for (k in 1:2)
                phenoFinder.args$eset.pair[[k]] <- phenoFinder.args$eset.pair[[k]][0, ]
            ## do not use "smoking guns" when calculating phenotype distances
            for (k in 1:2)
                if(any(manual.smokingguns %in% colnames(pData(phenoFinder.args$eset.pair[[k]]))))
                    pData(phenoFinder.args$eset.pair[[k]]) <-
                        pData(phenoFinder.args$eset.pair[[k]])[, (!colnames(pData(phenoFinder.args$eset.pair[[k]])) %in% manual.smokingguns)]
            if(!is.null(cache.dir)){
                cache.file <- paste(cache.dir, "/", digest::digest(list(phenoFinder, phenoFinder.args)), ".rda", sep="")
                if(file.exists(cache.file)) {
		    if (verbose) message("\tSkipping phenoFinder, loading cached results.")
                    load(cache.file)
		}
            }
            if(!exists("pheno.sim")){
                if(verbose) message("Calculating phenotype similarities...")
                pheno.sim <- do.call(phenoFinder, phenoFinder.args)
            }
            if(!is.null(cache.dir) && !file.exists(cache.file))
                save(pheno.sim, file=cache.file)
            ## find phenotype doppelgangers
            outlierFinder.pheno.args$similarity.mat <- pheno.sim
            outlierFinder.pheno.args$prune.output <- FALSE
            if(verbose) message("Identifying phenotype doppelgangers...")
            output3[["pheno.doppels"]] <- do.call(outlierFinder, outlierFinder.pheno.args)
            ## if(nrow(output3[["pheno.doppels"]]) == 0){ #deprecated
            ##     output3[["pheno.doppels"]] <- output3[["expr.doppels"]][, 1:2]
            ##     output3[["pheno.doppels"]]$similarity <- NA
            ##     output3[["pheno.doppels"]]$doppel <- FALSE
            ## }
        }
        ##If a sample is identified as a doppelganger by any method,
        ##then keep that sample for all methods, so we don't lose the
        ##data when wrapping up:
        keep.rows <- output3[["pheno.doppels"]]$doppel | output3[["expr.doppels"]]$doppel | output3[["smokinggun.doppels"]]$doppel
        output3[["smokinggun.doppels"]] <- output3[["smokinggun.doppels"]][keep.rows, ]
        output3[["pheno.doppels"]] <- output3[["pheno.doppels"]][keep.rows, ]
        output3[["expr.doppels"]] <- output3[["expr.doppels"]][keep.rows, ]
        return(output3)
    })

    names(output.full) <- sapply(ds.combns, function(ij) paste(names(esets)[c(ij[1], ij[2])], collapse=separator))
    if (verbose) message("Finalizing...")
    wrapUp <- function(object, element){
        tmp <- lapply(object, function(x) x[[element]])
        tmp <- tmp[!sapply(tmp, is.null)]
    ##    tmp <- lapply(tmp, function(x) x[x[, 4], ])
        do.call(rbind, tmp)
    }
    pheno.doppels <- wrapUp(output.full, "pheno.doppels")
    expr.doppels <- wrapUp(output.full, "expr.doppels")
    smokinggun.doppels <- wrapUp(output.full, "smokinggun.doppels")
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
        match.rows2 <- match(paste(add[, 2], add[, 1]), paste(orig[, 1], orig[, 2]))
        match.rows[is.na(match.rows)] <- match.rows2[is.na(match.rows)]
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
    ## merge all doppelganger types
    all.doppels <- expr.doppels
    colnames(all.doppels)[3:4] <- paste("expr.", colnames(all.doppels)[3:4], sep="")
    all.doppels <- addCols(all.doppels, pheno.doppels)
    colnames(all.doppels)[5:6] <- paste("pheno.", colnames(all.doppels)[5:6], sep="")
    all.doppels <- addCols(all.doppels, smokinggun.doppels)
    colnames(all.doppels)[7:8] <- sub("pheno", "smokinggun", colnames(all.doppels)[5:6])
    for (k in c(4, 6, 8))
        all.doppels[, k][is.na(all.doppels[, k])] <- FALSE
    all.doppels <- all.doppels[all.doppels[, 4] | all.doppels[, 6] | all.doppels[, 8], ]
    rownames(all.doppels) <- NULL
    ##If all esets have the same colnames of pdata, add merged
    ##pairwise sample pdata to all.doppels:
    if( identical(nrow(all.doppels) > 0, TRUE) &&
       identical(colnames(pData(esets[[1]])), unique(unlist(lapply(esets, function(eset) colnames(pData(eset)))))) ){
        all.pdat <- lapply(esets, pData)
        for (k in 1:length(esets))
            all.pdat[[k]]$sampleid <- paste(names(esets)[k], rownames(all.pdat[[k]]), sep=separator)
        all.pdat <- do.call(rbind, all.pdat)
        pdat.sample1 <- all.pdat[match(all.doppels[, 1], all.pdat$sampleid), ]
        pdat.sample2 <- all.pdat[match(all.doppels[, 2], all.pdat$sampleid), ]
        merged.pdat <- lapply(1:ncol(pdat.sample1), function(k) paste(pdat.sample1[, k], pdat.sample2[, k], sep=separator))
        merged.pdat <- do.call(cbind, merged.pdat)
        merged.pdat <- merged.pdat[, -ncol(merged.pdat)]
        colnames(merged.pdat) <- colnames(pData(esets[[1]]))
        all.doppels <- cbind(all.doppels, merged.pdat)
    }
    new("DoppelGang", fullresults=output.full, summaryresults=all.doppels, inputargs=input.args)
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
        X$alt_sample_name <- paste(X$sample_type, gsub("[^0-9]", "", X$alt_sample_name), sep="_")
        pData(X) <- pData(X)[, !grepl("uncurated_author_metadata", colnames(pData(X)))]
        return(X) })

    results1 <- doppelgangR(testesets, corFinder.args=list(use.ComBat=TRUE), cache.dir=NULL)
    summary(results1)
    plot(results1)
    ## Set phenoFinder.args to ignore similar phenotypes:
    results2 <- doppelgangR(testesets, corFinder.args=list(use.ComBat=FALSE), phenoFinder.args=NULL, cache.dir=NULL)
    summary(results2)
})
