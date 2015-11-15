library(curatedBladderData)
library(Biobase)
library(BiocParallel)
multicoreParam <- MulticoreParam()

source(system.file("extdata",
    "patientselection_all.config",package="curatedBladderData"))
min.number.of.genes <- 2000
source(system.file("extdata", "createEsetList.R", package =
    "curatedBladderData"))
rm(list=ls(pattern="_eset"))

table(table(unlist(lapply(esets, sampleNames))))

library(doppelgangR)
save(esets, file="bladder_esets.rda", compress="bzip2")

dop <- doppelgangR(esets, phenoFinder.args=NULL, smokingGunFinder.args=NULL,
                   outlierFinder.expr.args=list(bonf.prob = 1.0, transFun = atanh, tail = "upper"))
#dop <- doppelgangR(esets)
save(dop, file="bladder_dop_1.0.rda")

##Look for smoking guns only:
for (i in 1:length(esets))
    esets[[i]]$samplenames <- sampleNames(esets[[i]])
dop.gun <- doppelgangR(esets, manual.smokingguns=c("alt_sample_name", "unique_patient_ID", "samplenames"),
                       phenoFinder.args=NULL, corFinder.args=NULL, impute.knn.args=NULL)
save(dop.gun, file="bladder_dopgun_1.0.rda")
write.csv(summary(dop.gun), file="bladder_dopgun_1.0.csv")

load("bladder_dop_1.0.rda")
write.csv(dop@summaryresults, file="bladder_dop_1.0.csv")

pdf("bladder_dop_1.0.pdf")
plot(dop, skip.no.doppels=TRUE)
dev.off()
