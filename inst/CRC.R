library(curatedCRCData)
library(affy)

source(system.file("extdata", 
    "patientselection_all.config",package="curatedCRCData"))
min.number.of.genes <- 2000
source(system.file("extdata", "createEsetList.R", package =
    "curatedCRCData"))
rm(list=ls(pattern="_eset"))

library(doppelgangR)
##esets <- esets[c("GSE11237_eset", "GSE14095_eset")]
##esets <- esets[c("GSE11237_eset", "GSE3294_eset")]
esets <- esets[-match(c("GSE3294_eset","TCGA.RNASeqV2_eset", "TCGA.RNASeqV2.READ_eset"), names(esets))]  ##uses Genbank IDs instead of gene symbols.
dop <- doppelgangR(esets, phenoFinder.args=NULL, smokingGunFinder.args=NULL)
warnings()
#dop <- doppelgangR(esets)
save(dop, file="crc_dop.rda")