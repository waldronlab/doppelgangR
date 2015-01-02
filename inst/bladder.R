library(curatedBladderData)
library(affy)
source(system.file("extdata",
    "patientselection_all.config",package="curatedBladderData"))
min.number.of.genes <- 2000
source(system.file("extdata", "createEsetList.R", package =
    "curatedBladderData"))
rm(list=ls(pattern="_eset"))

library(doppelgangR)

dop <- doppelgangR(esets, phenoFinder.args=NULL, smokingGunFinder.args=NULL)
#dop <- doppelgangR(esets)
save(dop, file="bladder_dop.rda")

load("bladder_dop.rda")
write.csv(dop@summaryresults, file="bladder_dop.csv")

pdf("bladder_dop.pdf")
plot(dop, skip.no.doppels=TRUE)
dev.off()
