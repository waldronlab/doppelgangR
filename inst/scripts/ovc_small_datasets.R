library(Biobase)

if(file.exists("ovarian_esets.rda")){
    load("ovarian_esets.rda")
}else{
    library(curatedOvarianData)
    library(logging)
    source(system.file("extdata",
                       "patientselection_all.config",package="curatedOvarianData"))
    min.number.of.genes <- 2000
    source(system.file("extdata", "createEsetList.R", package =
                       "curatedOvarianData"))
    rm(list=ls(pattern="_eset"))
    save(esets, file="ovarian_esets.rda")
}

esets <- esets[c("GSE14764_eset", "GSE2109_eset", "GSE32062.GPL6480_eset")]
table(table(unlist(lapply(esets, sampleNames))))

library(doppelgangR)
res.orig <- doppelgangR(esets, phenoFinder.args=NULL, smokingGunFinder.args=NULL, within.datasets.only=TRUE)
confirmed.dups <- summary(res.orig)[-7, 1:2]
dups <- c(confirmed.dups[, 1], confirmed.dups[, 2])
dups <- sub(".+:", "", dups)

reduceEset <- function(eset, dups, n){
    keeps <- which(sampleNames(eset) %in% dups)
    maybes <- which(!sampleNames(eset) %in% dups)
    n.extras <- n - length(keeps)
    keeps <- c(keeps, sample(maybes, n.extras))
    return(eset[, keeps])
}


if(file.exists("ovc_small_datasets.rds")){
    res <- readRDS("ovc_small_datasets.rds")
}else{
    res <- list()
    for (i in 1:length(esets)){
        res[[i]] <- list()
        nseq <- c(5, seq(10, ncol(esets[[i]]), by=10))
        min.keeps <- sum(sampleNames(esets[[i]]) %in% dups)
        nseq <- nseq[nseq > min.keeps]
        for (k in 1:length(nseq)){
            eset <- reduceEset(esets[[i]], dups, nseq[k])
            res[[i]][[k]] <- doppelgangR(eset, phenoFinder.args=NULL, smokingGunFinder.args=NULL, intermediate.pruning=FALSE)
        }
        names(res[[i]]) <- nseq
    }
    names(res) <- names(esets)
    saveRDS(res, file="ovc_small_datasets.rds")
}

for (i in 1:length(res)){
    print(" ")
    print(names(res)[i])
    print(sapply(res[[i]], function(x) nrow(summary(x))))
}

plot(res[["GSE32062.GPL6480_eset"]][[1]])
plot(res[["GSE32062.GPL6480_eset"]][[2]])
