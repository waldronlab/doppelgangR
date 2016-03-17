library(curatedCRCData)
library(Biobase)
library(BiocParallel)
multicoreParam <- MulticoreParam()

source(system.file("extdata",
    "patientselection_all.config",package="curatedCRCData"))
min.number.of.genes <- 2000
source(system.file("extdata", "createEsetList.R", package =
    "curatedCRCData"))
rm(list=ls(pattern="_eset"))

library(doppelgangR)
##esets <- esets[c("GSE11237_eset", "GSE14095_eset")]
##esets <- esets[c("GSE11237_eset", "GSE3294_eset")]
esets <- esets[!names(esets) %in% c("GSE3294_eset","TCGA.RNASeqV2_eset", "TCGA.RNASeqV2.READ_eset")]  ##uses Genbank IDs instead of gene symbols.
esets <- esets[!names(esets) %in% c("GSE17536_eset", "GSE17537_eset")] #subseries of GSE17538.GPL570_eset

## Double-check there are no duplicated sample IDs:
table(table(unlist(lapply(esets, sampleNames))))

save(esets, file="CRC_esets.rda", compress="bzip2")

if(file.exists("crc_dop.rda")){
  load("crc_dop.rda")
}else{
  dop <- doppelgangR(esets, phenoFinder.args=NULL, smokingGunFinder.args=NULL, outlierFinder.expr.args=list(bonf.prob = 1.0, transFun = atanh, tail = "upper"))
  warnings()
  save(dop, file="crc_dop.rda")
}

##Look for smoking guns only:
for (i in 1:length(esets))
    esets[[i]]$samplenames <- sampleNames(esets[[i]])
dop.gun <- doppelgangR(esets, manual.smokingguns=c("alt_sample_name", "unique_patient_ID", "samplenames"),
                       phenoFinder.args=NULL, corFinder.args=NULL, impute.knn.args=NULL)
save(dop.gun, file="crc_dopgun.rda")
write.csv(summary(dop.gun), file="crc_dopgun.csv")

load("crc_dop.rda")
write.csv(dop@summaryresults, file="crc_dop.csv")
pdf("CRC_dop.pdf")
plot(dop, skip.no.doppels=TRUE)
dev.off()
