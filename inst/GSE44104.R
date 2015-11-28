library(GEOquery)
library(doppelgangR)
eset <- getGEO("GSE44104")[[1]]
wuetal <- doppelgangR(eset, outlierFinder.expr.args = list(bonf.prob = 40, transFun = atanh, tail = "upper"))
wuetaltable <- summary(wuetal); rownames(wuetaltable) <- NULL
wuetaltable[, 1] <- paste0("GSE44104:", wuetaltable[, 1])
wuetaltable[, 2] <- paste0("GSE44104:", wuetaltable[, 2])
write.csv(wuetaltable, file="wuetal.csv")


library(pheatmap)
cormat <- wuetal@fullresults[[1]][[1]]
rownames(cormat) <- sub("ExpressionSet1:", "", rownames(cormat))
colnames(cormat) <- sub("ExpressionSet2:", "", colnames(cormat))

pdf("GSE44104.pdf", width=6.5, height=6.5)
plot(wuetal)
pheatmap(cormat, fontsize_row=6, fontsize_col=6)
dev.off()

if(!file.exists("GSE44104_esets.rda")){
  (fnames=dir(pattern="^.*\\.CEL\\.gz$"))
  batch.var <-
    sapply(fnames,function(fname)
    {library(affyio)
      tempDate1 <- read.celfile.header(fname,info="full")$ScanDate
      output <- as.character(as.Date(tempDate1, format = "%m/%d/%y %H:%M:%S"))
      if(length(tempDate1) > 0 && is.na(output))
        output <- as.character(as.Date(tempDate1))
      if(length(output)==0) output <- NA
      return(output)
    })
  batch.var=factor(batch.var)
  library(affy)
  affybatch <- read.affybatch(filenames=fnames)
  ##library(arrayQualityMetrics)
  ##aqm <- arrayQualityMetrics(affybatch, do.logtransform=TRUE)
  eset2 = mas5(affybatch)
  exprs(eset2) = log2(exprs(eset2))
  rm(affybatch); gc()
  library(sva)
  big.matrix.combat <- sva::ComBat(exprs(eset2), mod=model.matrix(~(rep(1, length(batch.var)))), batch=batch.var)
  eset3 <- eset2
  exprs(eset3) = big.matrix.combat
  save(eset2, eset3, file="GSE44104_esets.rda")
}else{
  load("GSE44104_esets.rda")
}

wuetal2 <- doppelgangR(eset2, outlierFinder.expr.args = list(bonf.prob = 13, transFun = atanh, tail = "upper"), 
                       phenoFinder.args = NULL, smokingGunFinder.args = NULL, intermediate.pruning = TRUE)

## eset3 has ComBat batch correction, does not help.
wuetal3 <- doppelgangR(eset3, outlierFinder.expr.args = list(bonf.prob = 5, transFun = atanh, tail = "upper"), 
                       phenoFinder.args = NULL, smokingGunFinder.args = NULL, intermediate.pruning = TRUE)
plot(wuetal3)

library(pheatmap)
pheatmap(cor(exprs(eset2)))
pheatmap(cor(exprs(eset3)))
