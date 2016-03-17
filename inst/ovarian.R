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
table(table(unlist(lapply(esets, sampleNames))))

library(doppelgangR)

if(!file.exists("ovarian_dop.rda")){
  dop <- doppelgangR(esets, phenoFinder.args=NULL, smokingGunFinder.args=NULL,
                   outlierFinder.expr.args=list(bonf.prob = 1.0, transFun = atanh, tail = "upper"))
  warnings()
  save(dop, file="ovarian_dop.rda")
}else{
  load("ovarian_dop.rda")
}

##Look for smoking guns only:
dop.gun <- doppelgangR(esets, manual.smokingguns="alt_sample_name", phenoFinder.args=NULL,
                       corFinder.args=NULL, impute.knn.args=NULL, cache.dir = NULL)
save(dop.gun, file="ovarian_dopgun.rda")
write.csv(summary(dop.gun), file="ovarian_dopgun.csv")

##Look for smoking guns only:
for (i in 1:length(esets))
    esets[[i]]$samplenames <- sampleNames(esets[[i]])
dop.gun <- doppelgangR(esets, manual.smokingguns=c("alt_sample_name", "unique_patient_ID", "samplenames"),
                       phenoFinder.args=NULL, corFinder.args=NULL, impute.knn.args=NULL)
save(dop.gun, file="ovarian_dopgun.rda")
write.csv(summary(dop.gun), file="ovarian_dopgun.csv")

load("ovarian_dop.rda")
write.csv(dop@summaryresults, file="ovarian_dop.csv")
pdf("ovarian_dop.pdf")
plot(dop, skip.no.doppels=TRUE)
dev.off()

wuetal <- doppelgangR(esets[["GSE44104_eset"]], outlierFinder.expr.args = list(bonf.prob = 40, transFun = atanh, tail = "upper"))
wuetaltable <- summary(wuetal); rownames(wuetaltable) <- NULL
wuetaltable[, 1] <- paste0("GSE41044_eset:", wuetaltable[, 1])
wuetaltable[, 2] <- paste0("GSE41044_eset:", wuetaltable[, 2])
write.csv(wuetaltable, file="wuetal.csv")


library(pheatmap)
cormat <- wuetal@fullresults[[1]][[1]]
rownames(cormat) <- sub("ExpressionSet1:", "", rownames(cormat))
colnames(cormat) <- sub("ExpressionSet2:", "", colnames(cormat))

pdf("GSE44104.pdf", width=6.5, height=6.5)
plot(wuetal)
pheatmap(cormat, fontsize_row=6, fontsize_col=6)
dev.off()

library(GEOquery)
eset=getGEO("GSE44104")[[1]]
wuetal2 <- doppelgangR(eset, outlierFinder.expr.args = list(bonf.prob = 40, transFun = atanh, tail = "upper"))
plot(wuetal2)
