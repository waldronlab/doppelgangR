library(BiocParallel)
multicoreParam <- MulticoreParam()

breast.packages <- c("breastCancerMAINZ", "breastCancerNKI", "breastCancerTRANSBIG", "breastCancerUNT", "breastCancerUPP", "breastCancerVDX")

other.packages <- "WGCNA"

if (!require(BiocInstaller))
    stop("You need to install Bioconductor, which includes BiocInstaller.")


for (pkg in breast.packages){
    if(!require(package=pkg, character.only=TRUE)){
        print(paste("Need to install", pkg))
        biocLite(pkg, suppressUpdates=TRUE, suppressAutoUpdate=TRUE, ask=FALSE)
    }
}


for (pkg in other.packages){
    if(!require(package=pkg, character.only=TRUE)){
        print(paste("Need to install", pkg))
        biocLite(pkg, suppressUpdates=TRUE, suppressAutoUpdate=TRUE, ask=FALSE)
    }
}

if(file.exists("breast_esets.rda")) {
  load("breast_esets.rda")
}else{
  esets <- bplapply(breast.packages, function(pkg) {
    print(pkg)
    library(Biobase)
    esetname <- tolower(sub("breastCancer", "", pkg))
    data(list = esetname)
    output <- get(esetname)
    output <- output[!is.na(featureData(output)$EntrezGene.ID),]
    merge.probeset <- WGCNA::collapseRows(
      datET = exprs(output),
      rowGroup = featureData(output)$EntrezGene.ID,
      rowID = featureNames(output)
    )
    output <- output[merge.probeset$selectedRow,]
    featureNames(output) <- featureData(output)$EntrezGene.ID
    return(output)
  })
  names(esets) <- sub("breastCancer", "", breast.packages)
  save(esets, file = "breast_esets.rda", compress = "bzip2")
}

library(Biobase)
table(table(unlist(lapply(esets, sampleNames))))

#set.seed(1)
#eset2 <- lapply(esets, function(x) x[sample(1:nrow(x), 300), ])

library(doppelgangR)

if(file.exists("breast_dop.rda")){
  load("breast_dop.rda")
}else{
  dop <- doppelgangR(esets, phenoFinder.args=NULL, smokingGunFinder.args=NULL,
                     outlierFinder.expr.args=list(bonf.prob = 1.0, transFun = atanh, tail = "upper"))
  warnings()
  save(dop, file="breast_dop.rda")
}
##Look for smoking guns only:
for (i in 1:length(esets))
    esets[[i]]$samplenames <- sampleNames(esets[[i]])
dop.gun <- doppelgangR(esets, manual.smokingguns=c("alt_sample_name", "unique_patient_ID", "samplenames"),
                       phenoFinder.args=NULL, corFinder.args=NULL, impute.knn.args=NULL)
save(dop.gun, file="breast_dopgun.rda")
write.csv(summary(dop.gun), file="breast_dopgun.csv")

write.csv(dop@summaryresults, file="breast_dop.csv")
pdf("breast_dop.pdf")
plot(dop)
dev.off()

