if( !require(RTCGAToolbox) ){
    library(devtools)
    install_github("LiNk-NY/RTCGAToolbox")
    library(RTCGAToolbox)
}

all.dates <- getFirehoseRunningDates()
all.datasets <- getFirehoseDatasets()

data.path <- "."
# data.file <- file.path(data.path, "TCGA.rda")
# if(file.exists(data.file)){
#     load(data.file)
# }else{
#     tcga.res <- list()
#     for (i in 1:length(all.datasets)){
#         (ds.name <- all.datasets[i])
#         if(!ds.name %in% names(tcga.res)){
#             res <- try(getFirehoseData(ds.name, runDate=all.dates[1], RNAseq_Gene=TRUE, RNAseq2_Gene_Norm=TRUE, mRNA_Array=TRUE))
#         }
#         if(!is(res, "try-error")){
#             tcga.res[[ds.name]] <- res
#         }
#     }
#     save(tcga.res, file=data.file)
# }


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
      arraydat = arraydat[,substr(sampleNames(arraydat),14,15) %in% c( "01", "03", "09" )] )
      seqdat = [RNAseq[[pickplat]],substr(sampleNames(RNAseq[[pickplat]]),14,15) %in% c( "01", "03", "09" )] )
      eset.list[[ names(tcga.res)[i] ]][["microarray"]] = arraydat
      eset.list[[ names(tcga.res)[i] ]][[paste(names(tcga.res)[i], names(RNAseq)[pickplat])]] <- seqdat
    }
    save(eset.list, file=file.path(data.path, "tcga.microarray_RNAseq.rda"))
}


library(doppelgangR)
if(file.exists("doppelgangR.microarray_RNAseq.rda")){
  load("doppelgangR.microarray_RNAseq.rda")
}else{
  doppelgangR.microarray_RNAseq <- list()
  for (i in 1:length(eset.list)){
    print(names(eset.list)[i])
    eset.pair = eset.list[[i]]
    exprs(eset.pair[[2]]) = log(exprs(eset.pair[[2]]) + 1)
    doppelgangR.microarray_RNAseq[[names(eset.list)[i]]] = doppelgangR(esets=eset.pair, phenoFinder.args = NULL, smokingGunFinder.args = NULL)
  }
  save( doppelgangR.microarray_RNAseq, file="doppelgangR.microarray_RNAseq.rda")
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

pdf("microarray_RNAseq.pdf", width=8, height=4)
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

tcgacodes <-
structure(list(Study.Abbreviation = c("GBM", "OV", "LUSC", "LAML",
"COAD", "KIRC", "LUAD", "READ", "BRCA", "STAD", "UCEC", "KIRP",
"HNSC", "LIHC", "LGG", "BLCA", "THCA", "CESC", "PRAD", "PAAD",
"DLBC", "SKCM", "SARC", "KICH", "ESCA", "UCS", "ACC", "MESO",
"PCPG", "UVM", "CHOL", "TGCT", "THYM"), Study.Name = c("Glioblastoma multiforme",
"Ovarian serous cystadenocarcinoma", "Lung squamous cell carcinoma",
"Acute Myeloid Leukemia", "Colon adenocarcinoma", "Kidney renal clear cell carcinoma",
"Lung adenocarcinoma", "Rectum adenocarcinoma", "Breast invasive carcinoma",
"Stomach adenocarcinoma", "Uterine Corpus Endometrial Carcinoma",
"Kidney renal papillary cell carcinoma", "Head and Neck squamous cell carcinoma",
"Liver hepatocellular carcinoma", "Brain Lower Grade Glioma",
"Bladder Urothelial Carcinoma", "Thyroid carcinoma", "Cervical SCC and endocervical AC",
"Prostate adenocarcinoma", "Pancreatic adenocarcinoma", "Diffuse Large B-cell Lymphoma",
"Skin Cutaneous Melanoma", "Sarcoma", "Kidney Chromophobe", "Esophageal carcinoma ",
"Uterine Carcinosarcoma", "Adrenocortical carcinoma", "Mesothelioma",
"Pheochromocytoma and Paraganglioma", "Uveal Melanoma", "Cholangiocarcinoma",
"Testicular Germ Cell Tumors", "Thymoma")), .Names = c("Study.Abbreviation",
"Study.Name"), row.names = c(2L, 10L, 24L, 26L, 29L, 33L, 35L,
43L, 48L, 49L, 50L, 52L, 55L, 56L, 79L, 87L, 88L, 89L, 92L, 107L,
136L, 180L, 218L, 226L, 254L, 302L, 304L, 353L, 366L, 416L, 427L,
429L, 430L), class = "data.frame")
