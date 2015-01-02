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



esets <- lapply(breast.packages, function(pkg){
    print(pkg)
    library(affy)
    esetname <- tolower(sub("breastCancer", "", pkg))
    data(list=esetname)
    output <- get(esetname)
    output <- output[!is.na(featureData(output)$EntrezGene.ID), ]
    merge.probeset <- WGCNA::collapseRows(datET=exprs(output),
                                          rowGroup=featureData(output)$EntrezGene.ID,
                                          rowID=featureNames(output))
    output <- output[merge.probeset$selectedRow, ]
    featureNames(output) <- featureData(output)$EntrezGene.ID
    return(output)
})
names(esets) <- sub("breastCancer", "", breast.packages)

save(esets, file="esets.rda")
load("esets.rda")

#set.seed(1)
#eset2 <- lapply(esets, function(x) x[sample(1:nrow(x), 300), ])

library(doppelgangR)

dop <- doppelgangR(esets)
save(dop, file="breast_dop.rda")
load("breast_dop.rda")



pdf("~/Dropbox/tmp/breastdoppel.pdf")
plot(dop)
dev.off()
system("evince ~/Dropbox/breastdoppel.pdf &")
