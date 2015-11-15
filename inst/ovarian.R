library(curatedOvarianData)
library(Biobase)
library(logging)
library(BiocParallel)
multicoreParam <- MulticoreParam()

source(system.file("extdata",
    "patientselection_all.config",package="curatedOvarianData"))
min.number.of.genes <- 2000
source(system.file("extdata", "createEsetList.R", package =
    "curatedOvarianData"))
rm(list=ls(pattern="_eset"))

table(table(unlist(lapply(esets, sampleNames))))

library(doppelgangR)
save(esets, file="ovarian_esets.rda")

dop <- doppelgangR(esets, phenoFinder.args=NULL, smokingGunFinder.args=NULL,
                   outlierFinder.expr.args=list(bonf.prob = 1.0, transFun = atanh, tail = "upper"))
#dop <- doppelgangR(esets)
save(dop, file="ovarian_dop_1.0.rda")

##Look for smoking guns only:
dop.gun <- doppelgangR(esets, manual.smokingguns="alt_sample_name", phenoFinder.args=NULL,
                       corFinder.args=NULL, impute.knn.args=NULL)
save(dop.gun, file="ovarian_dopgun_1.0.rda")
write.csv(summary(dop.gun), file="ovarian_dopgun_1.0.csv")

##Look for smoking guns only:
for (i in 1:length(esets))
    esets[[i]]$samplenames <- sampleNames(esets[[i]])
dop.gun <- doppelgangR(esets, manual.smokingguns=c("alt_sample_name", "unique_patient_ID", "samplenames"),
                       phenoFinder.args=NULL, corFinder.args=NULL, impute.knn.args=NULL)
save(dop.gun, file="ovarian_dopgun_1.0.rda")
write.csv(summary(dop.gun), file="ovarian_dopgun_1.0.csv")

load("ovarian_dop_1.0.rda")
write.csv(dop@summaryresults, file="ovarian_dop_1.0.csv")
pdf("ovarian_dop_1.0.pdf")
plot(dop, skip.no.doppels=TRUE)
dev.off()
