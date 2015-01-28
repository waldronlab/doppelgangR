library(curatedBladderData)
library(affy)
library(BiocParallel)
multicoreParam <- MulticoreParam()

source(system.file("extdata",
    "patientselection_all.config",package="curatedBladderData"))
min.number.of.genes <- 2000
source(system.file("extdata", "createEsetList.R", package =
    "curatedBladderData"))
rm(list=ls(pattern="_eset"))

library(doppelgangR)

dop <- doppelgangR(esets, phenoFinder.args=NULL, smokingGunFinder.args=NULL,
                   outlierFinder.expr.args=list(bonf.prob = 1.0, transFun = atanh, tail = "upper"))
#dop <- doppelgangR(esets)
save(dop, file="bladder_dop.rda")

load("bladder_dop.rda")
write.csv(dop@summaryresults, file="bladder_dop.csv")

pdf("bladder_dop.pdf")
plot(dop, skip.no.doppels=TRUE)
dev.off()
