## TCGA.R should be run first, as this script uses its product suitability.table.csv.


if( !require(RTCGAToolbox) ){
    library(devtools)
    install_github("LiNk-NY/RTCGAToolbox")
    library(RTCGAToolbox)
}

all.dates <- getFirehoseRunningDates()
all.datasets <- getFirehoseDatasets()

data.path <- "."
data.file <- file.path(data.path, "TCGA.rda")
if(file.exists(data.file)){
    load(data.file)
}else{
    tcga.res <- list()
    for (i in 1:length(all.datasets)){
        (ds.name <- all.datasets[i])
        if(!ds.name %in% names(tcga.res)){
            res <- try(getFirehoseData(ds.name, runDate=all.dates[1], RNAseq_Gene=TRUE, RNAseq2_Gene_Norm=TRUE, mRNA_Array=TRUE))
        }
        if(!is(res, "try-error")){
            tcga.res[[ds.name]] <- res
        }
    }
    save(tcga.res, file=data.file)
}


library(Biobase)
if(file.exists(file.path(data.path, "tcga.microarray_RNAseq.rda"))){
  load(file.path(data.path, "tcga.microarray_RNAseq.rda"))
}else{
  eset.list <- list()
    for (i in 1:length(tcga.res)){
      print(names(tcga.res)[i])
      RNAseq <- list()
      arraydat <- try(extract(tcga.res[[i]], "mRNA_Array"))
      if(is(arraydat, "try-error")) next
      RNAseq[["rnaseq"]] <- try(extract(tcga.res[[i]], "RNAseq_Gene"))
      RNAseq[["rnaseq2"]] <- try(extract(tcga.res[[i]], "RNAseq2_Gene_Norm"))
      if(is(RNAseq[["rnaseq"]], "try-error") && is(RNAseq[["rnaseq2"]], "try-error")) next
      RNAseq <- RNAseq[sapply(RNAseq, function(x) !is(x, "try-error"))]
      pickplat <- which.max(sapply(RNAseq, ncol))
      print(paste(names(tcga.res)[i], ":", names(RNAseq)[pickplat]))
      eset.list[[names(tcga.res)[i]]] = list()
      ## use only primary tumors
      arraydat = arraydat[, substr(sampleNames(arraydat),14,15) %in% c( "01", "03", "09" )]
      seqdat = RNAseq[[pickplat]][ ,substr(sampleNames(RNAseq[[pickplat]]),14,15) %in% c( "01", "03", "09" )]
      eset.list[[ names(tcga.res)[i] ]][["microarray"]] = arraydat
      eset.list[[ names(tcga.res)[i] ]][[paste(names(tcga.res)[i], names(RNAseq)[pickplat])]] <- seqdat
    }
    save(eset.list, file=file.path(data.path, "tcga.microarray_RNAseq.rda"))
}



tcga.remove <- read.delim(system.file("extdata", "TCGA_remove.txt",
                                      package="doppelgangR"), as.is=TRUE)[, 1]
tcga.remove <- tolower(tcga.remove)
eset.list$OV$microarray <- eset.list$OV$microarray[, !make.names(substr(sampleNames(eset.list$OV$microarray), 1, 12)) %in% tcga.remove]

library(doppelgangR)
if(file.exists("doppelgangR.microarray_RNAseq_nolog2.rda")){
  load("doppelgangR.microarray_RNAseq_nolog2.rda")
}else{
  doppelgangR.microarray_RNAseq <- list()
  for (i in 1:length(eset.list)){
    print(names(eset.list)[i])
    eset.pair = eset.list[[i]]
#    exprs(eset.pair[[2]]) = log(exprs(eset.pair[[2]]) + 1)
    doppelgangR.microarray_RNAseq[[names(eset.list)[i]]] = doppelgangR(esets=eset.pair, phenoFinder.args = NULL, smokingGunFinder.args = NULL)
  }
  save( doppelgangR.microarray_RNAseq, file="doppelgangR.microarray_RNAseq_nolog2.rda")
}

doppelmelt <- function(obj, ds1, ds2){
  if(paste(ds1, ds2, sep=":") %in% names(obj@fullresults)){
    ds <- paste(ds1, ds2, sep=":")
  }else if(paste(ds2, ds1, sep=":") %in% names(obj@fullresults)){
    ds <- paste(ds2, ds1, sep=":")
  }else{
    return(NULL)
  }
  cormat <- obj@fullresults[[ds]]$correlations
  if(nrow(cormat) < ncol(cormat)) cormat <- t(cormat)
  idx <- sapply(rownames(cormat), function(x) which.max(cormat[x, ]))
  corvec <- sapply(1:nrow(cormat), function(i) cormat[i, idx[i]])
  output <- data.frame(sample1=rownames(cormat),
                       sample2=colnames(cormat)[idx],
                       cor=corvec, stringsAsFactors=FALSE)
  output$truepos <- sub(".+:", "", output[, 1]) == sub(".+:", "", output[, 2])
  return(output)
}
plotROC <- function(pred, labels, plot = TRUE, na.rm = TRUE, colorize = FALSE, addtext=TRUE, ...) {
  require(ROCR)
  require(pROC)
  if (na.rm) {
    idx <- !is.na(labels)
    pred <- pred[idx]
    labels <- labels[idx]
  }
  pred.rocr <- ROCR::prediction(pred, labels)
  perf.rocr <- ROCR::performance(pred.rocr, "tpr", "fpr")
  auc <- performance(pred.rocr, "auc")@y.values[[1]][[1]]
  roc.obj <- roc(labels, pred)
  auc.ci <- ci(roc.obj)
  significant <- ifelse(ci(roc.obj, conf.level=0.9)[1] > 0.5, "*", "")
  best <- coords(roc.obj,x="best")
  if (plot) {
    plot(perf.rocr, colorize = colorize, cex.lab = 1.3, bty="n", lty=1:length(perf.rocr),...)
    abline(a = 0, b = 1, lty = 2)
    if(addtext){
      text(0, 0.9, paste("AUC = ", round(auc, digits = 2), significant,
                         sep=""), cex = 1.5, pos = 4)
      text(1, 0.1, paste("n =", length(labels)), cex = 1.5, pos = 2)
    }
  }
  invisible(list(auc,auc.ci,best))
}

pdf("microarray_RNAseq_nolog2.pdf", width=8, height=4)
par(mfrow=c(1, 2))
res <- list()
for (i in 1:length(doppelgangR.microarray_RNAseq)){
  exptnames <- names(doppelgangR.microarray_RNAseq[[i]]@fullresults)[[3]]
  exptnames <- strsplit(exptnames, ":")[[1]]
  roc1 <- doppelmelt(doppelgangR.microarray_RNAseq[[i]], exptnames[1], exptnames[2])
  res[[i]] <- plotROC(roc1$cor, roc1$truepos, main=names(doppelgangR.microarray_RNAseq)[[i]])
  plot(doppelgangR.microarray_RNAseq[[i]], plot.pair=exptnames)
}
names(res) <- names(doppelgangR.microarray_RNAseq)
dev.off()

pdf("microarray_RNAseq_doppelgangerout_nolog2.pdf", width=9, height=3)
par(mfrow=c(1, 3))
for (i in 1:length(doppelgangR.microarray_RNAseq)){
  plot(doppelgangR.microarray_RNAseq[[i]])
}
dev.off()

suitability.table <- read.csv("suitability.table.csv", row.names=1, as.is=TRUE)
suitability.table$AUC = NA
rocvec <- sapply(res, function(x) x[[1]])
suitability.table$AUC[match(names(rocvec), suitability.table$cancertype)] <- rocvec
write.csv(suitability.table[, c("cancertype", "Study.Name", "quantile999", "AUC")], file="TCGA_microarray_RNAseq_nolog2.csv")

