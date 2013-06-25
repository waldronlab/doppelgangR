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
outlierFinder.expr.args=list(bonf.pvalue=0.01, transFun=atanh, tail="upper"),
### a list of arguments to be passed to outlierFinder when called for expression data
outlierFinder.pheno.args=list(normal.upper.thresh=0.99, bonf.pvalue=NULL, tail="upper"),
### a list of arguments to be passed to outlierFinder when called for phenotype data
smokingGunFinder.args=list(transFun=I),
### a list of arguments to be passed to smokingGunFinder
manual.smokingguns=NULL,
### a character vector of phenoData columns that, if identical, will
### be considered evidence of duplication
automatic.smokingguns=TRUE,
### automatically look for "smoking guns."  If TRUE, look for
### phenotype variables that are unique to each patient in dataset 1,
### also unique to each patient in dataset 2, but contain exact
### matches between datasets 1 and 2.
verbose=TRUE
### print progress information.
){
    if (is.null(names(esets)))
        names(esets) <- make.names(1:length(esets))
    output.full <- lapply(1:length(esets), function(i){
        output2 <- lapply(i:length(esets), function(j, i){
            if (verbose)
                print(paste("Working on datasets", names(esets)[i], "and", names(esets)[j]))
            output3 <- list()
            ## calculate correlation matrix
            corFinder.args$eset.pair <- esets[c(i, j)]
            cor.sim <- do.call(corFinder, corFinder.args)
            output3[["correlations"]] <- cor.sim
            ## find numeric (expression) doppelgangers
            outlierFinder.expr.args$similarity.mat <- cor.sim
            outlierFinder.expr.args$prune.output <- FALSE
            output3[["expr.doppels"]] <-
                do.call(outlierFinder, outlierFinder.expr.args)
            ## find potential "smoking gun" phenotypes
            new.smokinggun.phenotypes <- unlist(sapply(colnames(pData(esets[[i]])), function(cname){
                if(cname %in% colnames(pData(esets[[j]]))){
                    if((sum(!is.na(pData(esets[[i]])[, cname])) > 2 &
                        sum(!is.na(pData(esets[[j]])[, cname])) > 2) &
                       (identical(length(pData(esets[[i]])[, cname]), length(unique(pData(esets[[i]])[, cname]))) |
                        identical(length(pData(esets[[j]])[, cname]), length(unique(pData(esets[[j]])[, cname]))))
                       )
                        return(cname)}}))
            if(automatic.smokingguns){
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
        names(output2) <- names(esets)[i:length(esets)]
        return(output2)
    })
    names(output.full) <- names(esets)
    pheno.doppels <- do.call(rbind, lapply(output.full, function(x) do.call(rbind, lapply(x, function(y) y$pheno.doppels))))
    expr.doppels <- do.call(rbind, lapply(output.full, function(x) do.call(rbind, lapply(x, function(y) y$expr.doppels))))
    smokinggun.doppels <- do.call(rbind, lapply(output.full, function(x) do.call(rbind, lapply(x, function(y) y$smokinggun.doppels))))
    ## merge all doppelganger types
    all.doppels <- expr.doppels
    colnames(all.doppels)[3:4] <- paste("expr.", colnames(all.doppels)[3:4], sep="")
    addCols <- function(orig, add){
        match.rows <- match(paste(add[, 1], add[, 2]), paste(orig[, 1], orig[, 2]))
        orig[[colnames(add)[3]]] <- NA
        orig[[colnames(add)[3]]][match.rows] <- add[, 3]
        orig[[colnames(add)[4]]] <- NA
        orig[[colnames(add)[4]]][match.rows] <- add[, 4]
        orig
    }
    all.doppels <- addCols(all.doppels, pheno.doppels)
    colnames(all.doppels)[5:6] <- paste("pheno.", colnames(all.doppels)[5:6], sep="")
    if(is.null(smokinggun.doppels)){
        all.doppels$smokinggun.similarity <- 0
        all.doppels$smokinggun.doppel <- FALSE
    }else{
        all.doppels <- addCols(all.doppels, smokinggun.doppels)
        colnames(all.doppels)[7:8] <- paste("smokinggun.", colnames(all.doppels)[7:8], sep="")
    }
    all.doppels <- all.doppels[all.doppels[, 4] | all.doppels[, 6] | all.doppels[, 8], ]
    new("DoppelGang", fullresults=output.full, summaryresults=all.doppels)
### Returns an object of S4-class "DoppelGang".  See ?DoppelGang-class.
}, ex=function(){
    library(curatedOvarianData)
    data(GSE32062.GPL6480_eset)
    data(GSE32063_eset)
    data(GSE12470_eset)
    data(GSE17260_eset)
    testesets <- list(JapaneseA=GSE32062.GPL6480_eset,
                      JapaneseB=GSE32063_eset, 
                      Yoshihara2009=GSE12470_eset, 
                      Yoshihara2010=GSE17260_eset)
    testesets <- lapply(testesets, function(X) { sampleNames(X) <-
                                                     X$alt_sample_name; X })
    library(doppelgangR)
    tmp <- doppelgangR(testesets, corFinder.args=list(use.ComBat=TRUE))
})
