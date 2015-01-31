if( !require(RTCGAToolbox) ){
    library(devtools)
    install_github("mksamur/RTCGAToolbox")
    library(RTCGAToolbox)
}

all.dates <- getFirehoseRunningDates()
all.datasets <- getFirehoseDatasets()

tcga.res <- list()
for (i in 1:length(all.datasets)){
    (ds.name <- all.datasets[i])
    if(!ds.name %in% names(tcga.res))
        tcga.res[[ds.name]] <- getFirehoseData(ds.name, runDate=all.dates[1], RNAseq_Gene=TRUE, RNAseq2_Gene_Norm=TRUE, mRNA_Array=TRUE)
}

save(tcga.res, file="TCGA.rda")
