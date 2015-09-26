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
save(esets, file="CRC_esets.rda", compress="bzip2")
dop <- doppelgangR(esets, phenoFinder.args=NULL, smokingGunFinder.args=NULL, outlierFinder.expr.args=list(bonf.prob = 1.0, transFun = atanh, tail = "upper"))
warnings()
save(dop, file="crc_dop_1.0.rda")

##Look for smoking guns only:
dop.gun <- doppelgangR(esets, manual.smokingguns="alt_sample_name", phenoFinder.args=NULL, 
                       corFinder.args=NULL, impute.knn.args=NULL)
save(dop.gun, file="crc_dopgun_1.0.rda")
write.csv(summary(dop.gun), file="crc_dopgun_1.0.csv")

load("crc_dop_1.0.rda")
write.csv(dop@summaryresults, file="crc_dop_1.0.csv")
pdf("CRC_dop_1.0.pdf")
plot(dop, skip.no.doppels=TRUE)
dev.off()
