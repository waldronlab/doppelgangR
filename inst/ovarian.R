library(curatedOvarianData)
library(affy)
source(system.file("extdata",
    "patientselection_all.config",package="curatedOvarianData"))
min.number.of.genes <- 2000
source(system.file("extdata", "createEsetList.R", package =
    "curatedOvarianData"))
rm(list=ls(pattern="_eset"))

library(doppelgangR)

dop <- doppelgangR(esets, phenoFinder.args=NULL, smokingGunFinder.args=NULL)
#dop <- doppelgangR(esets)
save(dop, file="ovarian_dop.rda")
load("ovarian_dop.rda")
write.csv(dop@summaryresults, file="ovarian_dop.csv")

pdf("ovarian_dop.pdf")
plot(dop, skip.no.doppels=TRUE)
dev.off()
