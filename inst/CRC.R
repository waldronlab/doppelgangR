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
dop <- doppelgangR(esets, phenoFinder.args=NULL, smokingGunFinder.args=NULL)
#dop <- doppelgangR(esets)
save(dop, file="crc_dop.rda")
