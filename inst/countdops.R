library(gdata)
library(R.utils)

snames <- paste(c("bladder", "CRC", "breast", "ovarian"), "_dop.csv", sep="")

if(file.exists("manual_duplicates_curation.rda")){
    load("manual_duplicates_curation.rda")
}else{
    if(file.exists("manual_duplicates_curation.xls.bz2"))
        bunzip2("manual_duplicates_curation.xls.bz2", "manual_duplicates_curation.xls")
    manuals <- lapply(snames, function(x) read.xls("manual_duplicates_curation.xls", sheet=x, as.is=TRUE))
    names(manuals) <- snames
    save(manuals, file="manual_duplicates_curation.rda")
}

autos <- lapply(snames, function(x) read.csv(x, as.is=TRUE))
names(autos) <- snames

output <- sapply(1:length(autos), function(i){
    auto.ids <- with(autos[[i]], paste(sample1, sample2))
    if(grepl("breast", names(manuals)[i])){
        is.healthy <- rep(FALSE, nrow(manuals[[i]]))
    }else{
        is.healthy <- grepl("healthy|adjacentnormal", manuals[[i]]$sample_type)
    }
    manual.ids.healthy <- with(manuals[[i]][is.healthy & manuals[[i]]$manually.confirmed, ], paste(sample1, sample2))
    manual.ids.tumor <- with(manuals[[i]][!is.healthy & manuals[[i]]$manually.confirmed, ], paste(sample1, sample2))
    ##Dataframe for counting within and between-dataset:
    manual.ds <- strsplit(manual.ids.tumor, split=" ")
    manual.ds <- t(sapply(manual.ds, function(x) sub(":.*$", "", x)))
    auto.ds <- strsplit(auto.ids, split=" ")
    auto.ds <- t(sapply(auto.ds, function(x) sub(":.*$", "", x)))
    df <- data.frame(truepos=sum(auto.ids %in% manual.ids.tumor),
                     falsepos=sum(!auto.ids %in% manual.ids.tumor),
                     falsepos.is.healthy=sum(!auto.ids %in% manual.ids.healthy),
                     totalautopos=length(auto.ids),
                     within.ds.doppels=sum(manual.ds[, 1] == manual.ds[, 2]),
                     between.ds.doppels=sum(manual.ds[, 1] != manual.ds[, 2]),
                     total.dspairs.with.dops=length(unique(paste(manual.ds[, 1], manual.ds[, 2]))),
                     total.dspairs.found=sum(unique(paste(manual.ds[, 1], manual.ds[, 2])) %in% paste(auto.ds[, 1], auto.ds[, 2])),
                     total.manual.pos=length(manual.ids.tumor)
                     )
    rownames(df) <- names(autos)[i]
    df
})
colnames(output) <- names(autos)
output <- t(output)
output

write.csv(output, file="doppelgangers_summary.csv")

## summarize all samples
library(Biobase)
output <- list()
esetsfiles <- dir(pattern="_esets\\.rda$")
k <- 0
for (i in 1:length(esetsfiles)){
    load(esetsfiles[i])
    for (j in 1:length(esets)){
        k <- k+1
        output[[k]] <- c(sub("_esets.rda", "", esetsfiles[i]),
                         names(esets)[j],
                         ncol(esets[[j]]),
                         experimentData(esets[[j]])@lab,
                         esets[[j]]@annotation,
                         experimentData(esets[[j]])@pubMedIds,
                         experimentData(esets[[j]])@name)
    }
}

output <- do.call(rbind, output)
colnames(output) <- c("Cancer Type", "Dataset", "N.samples", "Lab", "Annotation", "PMID", "Name")
write.csv(output, file="AllStudies.csv")

## -----------------------------------------
## -----------------------------------------
## ROC plot from all correlations
## -----------------------------------------
## -----------------------------------------

## Curation of smoking gun doppelgangers

## Ovarian:
## GSE12470_eset GSE17260_eset
## GSE12470_eset GSE32062.GPL6480_eset
## GSE12470_eset GSE32063_eset
## PMID15897565_eset:X0074_1776_h133a_1784	PMID17290060_eset:X0074_1776_h133a_1784
## PMID17290060_eset:X0193_0000_h133a_D1805	PMID19318476_eset:D1805

library(Biobase)
load("ovarian_esets.rda")
esets <- esets[c("GSE12470_eset", "GSE17260_eset", "GSE32062.GPL6480_eset", "GSE32063_eset", "PMID15897565_eset", "PMID17290060_eset", "PMID19318476_eset")]
##
esets[["GSE12470_eset"]]$uniqueid <- paste("Yoshihara", esets[["GSE12470_eset"]]$sample_type, gsub("[^0-9]", "", esets[["GSE12470_eset"]]$alt_sample_name), sep="_")
esets[["GSE32062.GPL6480_eset"]]$uniqueid <- paste("Yoshihara", esets[["GSE32062.GPL6480_eset"]]$sample_type, gsub("[^0-9]", "", esets[["GSE32062.GPL6480_eset"]]$alt_sample_name), sep="_")
esets[["GSE32063_eset"]]$uniqueid <- paste("Yoshihara", esets[["GSE32063_eset"]]$sample_type, gsub("[^0-9]", "", esets[["GSE32063_eset"]]$alt_sample_name), sep="_")
esets[["GSE17260_eset"]]$uniqueid <- paste("Yoshihara", esets[["GSE17260_eset"]]$sample_type, gsub("[^0-9]", "", esets[["GSE17260_eset"]]$alt_sample_name), sep="_")
esets[["PMID15897565_eset"]]$uniqueid <- paste("Duke", esets[["PMID15897565_eset"]]$sample_type, sub(".+_", "", sampleNames(esets[["PMID15897565_eset"]])), sep="_")
esets[["PMID17290060_eset"]]$uniqueid <- paste("Duke", esets[["PMID17290060_eset"]]$sample_type, gsub("[^0-9]", "", sub(".+_", "", sampleNames(esets[["PMID17290060_eset"]]))), sep="_")
esets[["PMID17290060_eset"]]$uniqueid[esets[["PMID17290060_eset"]]$alt_sample_name=="M3142"] <- "M3142"
esets[["PMID19318476_eset"]]$uniqueid <- paste("Duke", esets[["PMID19318476_eset"]]$sample_type, gsub("[^0-9]", "", esets[["PMID19318476_eset"]]$alt_sample_name), sep="_")
library(doppelgangR)
tmp <- doppelgangR(esets, corFinder.args=NULL, phenoFinder.args=NULL, manual.smokingguns="uniqueid", cache.dir=NULL)
rm(esets)
goldstandard <- list()
goldstandard[["ovarian_dop.csv"]] <- rbind(paste(summary(tmp)[, 1], summary(tmp)[, 2]),
                                           paste(summary(tmp)[, 2], summary(tmp)[, 1]))

## ## Breast:
## TRANSBIG:VDXOXFU_104	UNT:OXFU_104
## TRANSBIG:VDXKIU_136B04	UPP:UPP_136B04
## UNT:KIU_101B88	UPP:UPP_101B88
load("breast_esets.rda")
esets <- esets[c("TRANSBIG", "UNT", "UPP")]
esets[["UNT"]]$uniqueid <- sub(".+_", "", sampleNames(esets[c("UNT")]))
esets[["UPP"]]$uniqueid <- sub(".+_", "", sampleNames(esets[c("UPP")]))
esets[["TRANSBIG"]]$uniqueid <- make.unique(sub(".+_", "", sampleNames(esets[c("TRANSBIG")])))
tmp <- doppelgangR(esets, corFinder.args=NULL, phenoFinder.args=NULL, manual.smokingguns="uniqueid", cache.dir=NULL)
goldstandard[["breast_dop.csv"]] <- rbind(paste(summary(tmp)[, 1], summary(tmp)[, 2]),
                                          paste(summary(tmp)[, 2], summary(tmp)[, 1]))
rm(esets)


## ## CRC:
## GSE13067_eset	GSE14333_eset:  41 manually
## GSE18105_eset	GSE21510_eset
load("CRC_esets.rda")

esets <- esets[c("GSE13067_eset", "GSE14333_eset", "GSE18105_eset", "GSE21510_eset")]
esets[["GSE13067_eset"]]$uniqueid <- paste("Jorissen", esets[["GSE13067_eset"]]$sample_type, sub(": .+", "", esets[["GSE13067_eset"]]$alt_sample_name), sep="_")
esets[["GSE14333_eset"]]$uniqueid <- paste("Jorissen", esets[["GSE14333_eset"]]$sample_type, sub("-[0-9][A-Z]$", "", esets[["GSE14333_eset"]]$alt_sample_name), sep="_")
##
celltype <- sub(".+, ", "", esets[["GSE18105_eset"]]$alt_sample_name)
esets[["GSE18105_eset"]]$uniqueid <- paste("Matsuyama", esets[["GSE18105_eset"]]$sample_type, celltype, gsub("[^0-9]", "", esets[["GSE18105_eset"]]$alt_sample_name), sep="_")
##
celltype <- sub(".+, ", "", esets[["GSE21510_eset"]]$alt_sample_name)
celltype <- sub(" .+", "", celltype)
esets[["GSE21510_eset"]]$uniqueid <- paste("Matsuyama", esets[["GSE21510_eset"]]$sample_type, celltype, gsub("[^0-9]", "", esets[["GSE21510_eset"]]$alt_sample_name), sep="_")

crc.dop <- doppelgangR(esets, corFinder.args=NULL, phenoFinder.args=NULL, manual.smokingguns="uniqueid", cache.dir=NULL)
goldstandard[["CRC_dop.csv"]] <- rbind(paste(summary(crc.dop)[, 1], summary(crc.dop)[, 2]),
                                       paste(summary(crc.dop)[, 2], summary(crc.dop)[, 1]))
rm(esets)



## ## Bladder:
## GSE5287_eset	GSE89_eset
## GSE19915.GPL3883_eset	GSE32894_eset
## GSE19915.GPL3883_eset GSE19915.GPL5186_eset

load("bladder_esets.rda")

esets <- esets[c("GSE5287_eset", "GSE89_eset", "GSE19915.GPL3883_eset", "GSE32894_eset", "GSE19915.GPL5186_eset")]
##
esets[["GSE5287_eset"]]$uniqueid <- paste("Als", esets[["GSE5287_eset"]]$sample_type, sub("Sample # ", "", esets[["GSE5287_eset"]]$alt_sample_name))
esets[["GSE89_eset"]]$uniqueid <- paste("Als", esets[["GSE89_eset"]]$sample_type, sub("Sample # ", "", esets[["GSE89_eset"]]$alt_sample_name))
##
esets[["GSE19915.GPL3883_eset"]]$uniqueid <- paste("Lindgren", esets[["GSE19915.GPL3883_eset"]]$sample_type, sub("([UCN]+_[0-9]+).+", "\\1", esets[["GSE19915.GPL3883_eset"]]$alt_sample_name))
esets[["GSE32894_eset"]]$uniqueid <- paste("Lindgren", esets[["GSE32894_eset"]]$sample_type, sub("([UCN]+_[0-9]+).+", "\\1", esets[["GSE32894_eset"]]$alt_sample_name))
esets[["GSE19915.GPL5186_eset"]]$uniqueid <- paste("Lindgren", esets[["GSE19915.GPL5186_eset"]]$sample_type, sub("([UCN]+_[0-9]+).+", "\\1", esets[["GSE19915.GPL5186_eset"]]$alt_sample_name))

summary(esets[["GSE19915.GPL5186_eset"]]$uniqueid %in% esets[["GSE19915.GPL3883_eset"]]$uniqueid)
summary(esets[["GSE32894_eset"]]$uniqueid %in% esets[["GSE19915.GPL3883_eset"]]$uniqueid)
summary(esets[["GSE32894_eset"]]$uniqueid %in% esets[["GSE19915.GPL5186_eset"]]$uniqueid)

bladder.dop <- doppelgangR(esets, corFinder.args=NULL, phenoFinder.args=NULL, manual.smokingguns="uniqueid", cache.dir=NULL)
goldstandard[["bladder_dop.csv"]] <- rbind(paste(summary(bladder.dop)[, 1], summary(bladder.dop)[, 2]),
                                           paste(summary(bladder.dop)[, 2], summary(bladder.dop)[, 1]))
rm(esets)


matrixToDf <- function(x){
    data.frame(sample=array(outer(rownames(x), colnames(x), FUN=paste)),
               correlation=as.numeric(x))
}

doppelmelt <- function(obj, ds){
    cormat <- obj@fullresults[[ds]]$correlations
    if(identical(rownames(cormat), colnames(cormat))){
        cormat[lower.tri(cormat)] <- cormat[upper.tri(cormat)]
    }
    if(nrow(cormat) < ncol(cormat)) cormat <- t(cormat)
    idx <- sapply(rownames(cormat), function(x) which.max(cormat[x, ]))
    corvec <- sapply(1:nrow(cormat), function(i) cormat[i, idx[i]])
    output <- data.frame(sample=paste(rownames(cormat), colnames(cormat)[idx]),
                         correlation=corvec, stringsAsFactors=FALSE)
    return(output)
}
plotROC <- function(pred, labels, plot = TRUE, na.rm = TRUE, colorize = FALSE, ...) {
    require(ROCR)
    require(pROC)
    if (na.rm) {
        idx <- !is.na(labels)
        pred <- pred[idx]
        labels <- labels[idx]
    }
    pred.rocr <- prediction(pred, labels)
    perf.rocr <- performance(pred.rocr, "tpr", "fpr")
    auc <- performance(pred.rocr, "auc")@y.values[[1]][[1]]
    roc.obj <- roc(labels, pred)
    auc.ci <- ci(roc.obj)
    significant <- ifelse(ci(roc.obj, conf.level=0.9)[1] > 0.5, "*", "")
    best <- coords(roc.obj,x="best", transpose = FALSE)
    if (plot) {
        plot(perf.rocr, colorize = colorize, cex.lab = 1.3, ...)
        text(0, 0.9, paste("AUC = ", round(auc, digits = 2), significant,
        sep=""), cex = 1.5, pos = 4)
        abline(a = 0, b = 1, lty = 2)
        text(1, 0.1, paste("n =", length(labels)), cex = 1.5, pos = 2)
        abline(a = 0, b = 1, lty = 2)
    }
    invisible(list(auc,auc.ci,best))
}



library(ROCR)
if(file.exists("pairwise.perf.rda")){
    load("pairwise.perf.rda")
}else{
    pairwise.perf <- list()
    pairwise.auc <- list()
    for (sname in snames){
        print(sname)
        load(sub("csv", "rda", sname))
##        pairwise.cors <- lapply(dop@fullresults, function(x) na.omit(matrixToDf(x$correlations)))
        pairwise.dsnames <- names(dop@fullresults)
        ## Look across studies only
        pairwise.dsnames <- pairwise.dsnames[sapply(strsplit(pairwise.dsnames, ":"), function(x) !identical(x[1], x[2]))]
        pairwise.cors <- lapply(pairwise.dsnames, function(ds){
            doppelmelt(dop, ds)
        })
        names(pairwise.cors) <- pairwise.dsnames
##        pairwise.cors <- do.call(rbind, pairwise.cors)
        ##
        pairwise.perf[[sname]] <- list()
        pairwise.auc[[sname]] <- list()
        for (i in 1:length(pairwise.cors)){
            pairwise.cors[[i]]$TP <- 0
            pairwise.cors[[i]]$TP[ pairwise.cors[[i]]$sample %in% goldstandard[[sname]] ] <- 1
            if(identical(all.equal(sum(pairwise.cors[[i]]$TP), 0), TRUE)) next
            pairwise.pred <- prediction(predictions=pairwise.cors[[i]]$correlation, labels=pairwise.cors[[i]]$TP)
            pairwise.perf[[sname]][[names(pairwise.cors)[i]]] <- performance(pairwise.pred, measure="tpr", x.measure="fpr")
            pairwise.auc[[sname]][[names(pairwise.cors)[i]]] <- performance(pairwise.pred, measure="auc")
        }
        ##
        ##        pairwise.perf[[sname]] <- roc(predictor=pairwise.cors$correlation, response=pairwise.cors$TP, plot=TRUE)
    }
    save(pairwise.perf, pairwise.auc, file="pairwise.perf.rda")
}

pdf("pairwise.perf.pdf")
par(mfrow=c(2, 2))
for (i in 1:length(pairwise.perf)){
    for (j in 1:length(pairwise.perf[[i]])){
        plot(pairwise.perf[[i]][[j]], add=(j>1), main=sub("_dop.csv", "", names(pairwise.perf)[i]), lty=j)
    }
    mycex <- ifelse(names(pairwise.perf)[i] %in% c("bladder_dop.csv", "ovarian_dop.csv"), 0.75, 1)
    legend("bottomright", legend=gsub("_eset", "", names(pairwise.perf[[i]])), lty=1:length(pairwise.perf[[i]]), bty='n', cex=mycex)
}
dev.off()

##Mean AUC is 0.97:
mean(unlist(sapply(pairwise.auc, function(auc1) sapply(auc1, function(auc2) auc2@y.values[[1]]))))

## pairwise.spline <- lapply(pairwise.perf, function(pairwise.cancer){
##     lapply(pairwise.cancer, function(perf1){
##         spline(perf1@x.values[[1]], perf1@y.values[[1]], xout=seq(0, 1, length.out=1001))
##     })
## })

## pairwise.splineout <- lapply(pairwise.spline, function(pairwise.1cancer){
##     data.frame(smoothx=rowMeans(sapply(pairwise.1cancer, function(x) x$x)),
##                smoothy=rowMeans(sapply(pairwise.1cancer, function(x) x$y)))
## })

## pdf("pairwise.perf.pdf")
## par(mfrow=c(2, 2))
## for (i in 1:4){
##     plot(pairwise.splineout[[i]], main=sub("_dop.csv", "", snames[i]), type="l", xlim=c(0, 1), ylim=c(0, 1))
##     abline(a=0, b=1, lty=2)
## }
## dev.off()


