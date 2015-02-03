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
    if(!ds.name %in% names(tcga.res)){
        res <- try(getFirehoseData(ds.name, runDate=all.dates[1], RNAseq_Gene=TRUE, RNAseq2_Gene_Norm=TRUE, mRNA_Array=TRUE))
    }
    if(!is(res, "try-error")){
        tcga.res[[ds.name]] <- res
    }
}

save(tcga.res, file="TCGA.rda")

extractRTCGA <- function(object, type){
    typematch <- match.arg(type,
          choices=c("RNAseq_Gene", "Clinic", "miRNASeq_Gene",
              "RNAseq2_Gene_Norm", "CNA_SNP", "CNV_SNP", "CNA_Seq",
              "CNA_CGH", "Methylation", "Mutation", "mRNA_Array",
              "miRNA_Array", "RPPA"))
    if(identical(typematch, "RNAseq_Gene")){
        output <- object@RNASeqGene
    }else if(identical(typematch, "RNAseq2_Gene_Norm")){
        output <- object@RNASeq2GeneNorm
    }else if(identical(typematch, "mRNA_Array")){
        if(is(object@mRNAArray, "FirehosemRNAArray")){
            output <- object@mRNAArray@Datamatrix
        }else if(is(object@mRNAArray, "list")){
            output <- lapply(object@mRNAArray, function(tmp){
                tmp@DataMatrix
            })
            warning(paste("Taking the mRNA_Array platform with the greatest number of samples:", which.max(sapply(output, ncol))))
            ## just silently take the platform with the greatest
            ## number of samples:
            output <- output[[which.max(sapply(output, ncol))]]
        }
    }else{
        stop(paste("Type", typematch, "not yet supported."))
    }
    return(output)
}

load("ov.array.rda")
dim(extractRTCGA(ov.array, "RNAseq_Gene"))
dim(extractRTCGA(ov.array, "RNAseq2_Gene_Norm"))
dim(extractRTCGA(ov.array, "mRNA_Array"))

load("TCGA.rda")

