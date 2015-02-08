if( !require(RTCGAToolbox) ){
    library(devtools)
    install_github("mksamur/RTCGAToolbox")
    library(RTCGAToolbox)
}

all.dates <- getFirehoseRunningDates()
all.datasets <- getFirehoseDatasets()

if(file.exists("tcga.res")){
    load("/scratch/lw391/doppelgangR/inst")
}else{
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
    save(tcga.res, file="/scratch/lw391/doppelgangR/inst/TCGA.rda")
}else{
    load("/scratch/lw391/doppelgangR/inst")
}

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
            if(length(output) == 0){
                output <- matrix(NA, nrow=0, ncol=0)
            }else if(length(output) == 1){
                output <- output[[1]]
            }else{
                ## just silently take the platform with the greatest
                ## number of samples:
                keeplist <- which.max(sapply(output, ncol))
                output <- output[[keeplist]]
                warning(paste("Taking the mRNA_Array platform with the greatest number of samples:", keeplist))
            }
        }
    }else{
        stop(paste("Type", typematch, "not yet supported."))
    }
    return(output)
}

library(affy)
if(file.exists("/scratch/lw391/doppelgangR/inst/tcga.esets.rda")){
    load("/scratch/lw391/doppelgangR/inst/tcga.esets.rda")
}else{
    tcga.esets <- list()
    for (i in 1:length(tcga.res)){
        print(names(tcga.res)[i])
        tmp <- list()
        tmp[["mrna"]] <- extractRTCGA(tcga.res[[i]], "mRNA_Array")
        tmp[["rnaseq"]] <- extractRTCGA(tcga.res[[i]], "RNAseq_Gene")
        tmp[["rnaseq2"]] <- extractRTCGA(tcga.res[[i]], "RNAseq2_Gene_Norm")
        pickplat <- which.max(sapply(tmp, ncol))
        tcga.esets[[paste(names(tcga.res)[i], names(tmp)[pickplat])]] <- ExpressionSet(tmp[[pickplat]])
    }
    save(tcga.esets, file="/scratch/lw391/doppelgangR/inst/tcga.esets.rda")
}

cor.list <- lapply(tcga.esets, function(eset){
    output <- cor(exprs(eset))
    output[upper.tri(output)]
})
names(cor.list) <- names(tcga.esets)

save(cor.list, file="/scratch/lw391/doppelgangR/inst/cor.list.rda")

load("/scratch/lw391/doppelgangR/inst/cor.list.rda")

ztrans.list <- lapply(cor.list, atanh)

suitability.table <- data.frame(perc.gt.95=signif(100*sapply(cor.list, function(x) sum(x > 0.95) / length(x)), 1),
                                quantile999=round(sapply(cor.list, quantile, 0.999), 2),
                                nsamples=sapply(cor.list, function(x) (1 + sqrt(1 + 4*2*length(x))) / 2))
rownames(suitability.table) <- sub("mrna", "microarray", rownames(suitability.table))
suitability.table$cancertype <- sub(" .+", "", rownames(suitability.table))
suitability.table$assaytype <- sub(".+ ", "", rownames(suitability.table))
rownames(suitability.table) <- NULL
##
tcgacodes <-
structure(list(Study.Abbreviation = c("GBM", "OV", "LUSC", "LAML",
"COAD", "KIRC", "LUAD", "READ", "BRCA", "STAD", "UCEC", "KIRP",
"HNSC", "LIHC", "LGG", "BLCA", "THCA", "CESC", "PRAD", "PAAD",
"DLBC", "SKCM", "SARC", "KICH", "ESCA", "UCS", "ACC", "MESO",
"PCPG", "UVM", "CHOL", "TGCT", "THYM"), Study.Name = c("Glioblastoma multiforme",
"Ovarian serous cystadenocarcinoma", "Lung squamous cell carcinoma",
"Acute Myeloid Leukemia", "Colon adenocarcinoma", "Kidney renal clear cell carcinoma",
"Lung adenocarcinoma", "Rectum adenocarcinoma", "Breast invasive carcinoma",
"Stomach adenocarcinoma", "Uterine Corpus Endometrial Carcinoma",
"Kidney renal papillary cell carcinoma", "Head and Neck squamous cell carcinoma",
"Liver hepatocellular carcinoma", "Brain Lower Grade Glioma",
"Bladder Urothelial Carcinoma", "Thyroid carcinoma", "Cervical squamous cell carcinoma and endocervical adenocarcinoma",
"Prostate adenocarcinoma", "Pancreatic adenocarcinoma", "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma",
"Skin Cutaneous Melanoma", "Sarcoma", "Kidney Chromophobe", "Esophageal carcinoma ",
"Uterine Carcinosarcoma", "Adrenocortical carcinoma", "Mesothelioma",
"Pheochromocytoma and Paraganglioma", "Uveal Melanoma", "Cholangiocarcinoma",
"Testicular Germ Cell Tumors", "Thymoma")), .Names = c("Study.Abbreviation",
"Study.Name"), row.names = c(2L, 10L, 24L, 26L, 29L, 33L, 35L,
43L, 48L, 49L, 50L, 52L, 55L, 56L, 79L, 87L, 88L, 89L, 92L, 107L,
136L, 180L, 218L, 226L, 254L, 302L, 304L, 353L, 366L, 416L, 427L,
429L, 430L), class = "data.frame")
##
suitability.table$Study.Name <- tcgacodes[match(suitability.table$cancertype, tcgacodes$Study.Abbreviation), "Study.Name"]
suitability.table <- suitability.table[, c(4, 6, 5, 3, 2:1)]
##
suitability.table <- suitability.table[order(suitability.table$quantile99, suitability.table$perc.gt.95), ]
rownames(suitability.table) <- 1:nrow(suitability.table)
suitability.table$Study.Name <- tolower(suitability.table$Study.Name)


library(xtable)
sink("suitability.table.html")
print(xtable(suitability.table), type="html")
sink()


pdf("TCGA_PairwisePearson.pdf")
j <- 0
for (i in match(suitability.table$cancertype, sub(" .+", "", names(cor.list)))){
    j <- j+1
    hist(cor.list[[i]], breaks="FD", xlab="PCC", xlim=c(0, 1),
         main=paste(names(cor.list)[i], "\n", suitability.table$Study.Name[j]))
    quant999 <- quantile(cor.list[[i]], 0.999)
    abline(v=quant999, col="red", lw=2); abline(v=0.95, col="black", lw=2)
    legend("topleft", pch=-1, col=c("black", "red"), bty="n",
           lw=2, legend=c("PCC=0.95", "99.9 percentile"))
}
dev.off()


##bimodal: KICH (RNAseq2), maybe KIRC and KIRP, LUSC, PAAD, PCPG, THCA, GBM

##very high correlations normally: THCA, HNSC, LIHC, maybe KIRP, LAML, PCPG, PRAD, STAD, THYM, ESCA

##looks good: ACC, BLCA, BRCA, CESC, COAD, COADREAD, DLBC, LGG,
