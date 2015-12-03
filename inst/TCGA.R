if( !require(RTCGAToolbox) ){
    library(devtools)
    install_github("LiNk-NY/RTCGAToolbox")
    library(RTCGAToolbox)
}

all.dates <- getFirehoseRunningDates()
all.datasets <- getFirehoseDatasets()

data.path <- "/scratch/lw391/doppelgangR/inst"
data.path <- "."
data.file <- file.path(data.path, "TCGA.rda")


if(file.exists(data.file)){
    load(data.file)
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
    save(tcga.res, file=data.file)
}

library(Biobase)
if(file.exists(file.path(data.path, "tcga.esets.rda"))){
    load(file.path(data.path, "tcga.esets.rda"))
}else{
    tcga.esets <- list()
    for (i in 1:length(tcga.res)){
        tmp <- list()
##        tmp[["mrna"]] <- try(extract(tcga.res[[i]], "mRNA_Array"))  #RNAseq and RNAseq2 only
        tmp[["rnaseq"]] <- try(extract(tcga.res[[i]], "RNAseq_Gene"))
        tmp[["rnaseq2"]] <- try(extract(tcga.res[[i]], "RNAseq2_Gene_Norm"))
        tmp <- tmp[sapply(tmp, function(x) !is(x, "try-error"))]
        pickplat <- which.max(sapply(tmp, ncol))
        print(paste(names(tcga.res)[i], ":", names(tmp)[pickplat]))
        tcga.esets[[paste(names(tcga.res)[i], names(tmp)[pickplat])]] <- tmp[[pickplat]]
    }
    save(tcga.esets, file=file.path(data.path, "tcga.esets.rda"))
}

# use only primary tumors
tcga.esets <- lapply(tcga.esets, function(x) x[,substr(sampleNames(x),14,15) %in% c( "01", "03", "09" )] )

cor.list <- lapply(tcga.esets, function(eset){
    output <- cor(exprs(eset))
    output[upper.tri(output)]
})
names(cor.list) <- names(tcga.esets)

save(cor.list, file=file.path(data.path, "cor.list.rda"))

load(file.path(data.path, "cor.list.rda"))

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
"Bladder Urothelial Carcinoma", "Thyroid carcinoma", "Cervical SCC and endocervical AC",
"Prostate adenocarcinoma", "Pancreatic adenocarcinoma", "Diffuse Large B-cell Lymphoma",
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
suitability.table$Study.Name[is.na(suitability.table$Study.Name )] <- "Colorectal adenocarcinoma"
suitability.table$Study.Name <- tolower(suitability.table$Study.Name)

##Last publication restrictions in place until 12/18/2015
# library(pipeR)
# library(XML)
# library(stringr)
# url1 <- "http://cancergenome.nih.gov/publications/publicationguidelines"
# tabb <- readHTMLTable(url1, stringsAsFactors=FALSE)
# tt <- tabb[1][[1]]
# names(tt) <- c("disease", "restrict")
# tt$disease %>>% strsplit("\\(") %>>% sapply("[", 2) %>>%
#   gsub(pattern="\\)", replacement="", x=.) %>>% tolower -> tt$disease
# tt <- tt[complete.cases(tt),]
# tt[tt[,"disease"]=="coad, read",1] <- c("coad")
# tt <- rbind(tt[1:7,],c("read", NA),tt[-(1:7),])
# imin <- which.max(which(is.na(tt[, 2])))
# #for (i in imin:nrow(tt))
# tt$restrict <- str_extract_all(tt$restrict, "[0-9]{2}/[0-9]{2}/[0-9]{4}")
# tt$restrict <- ifelse(str_detect(tt$restrict, "No restrictions"), NA, tt$restrict)
# 
# tt2 <- tt[match(suitability.table$cancertype, toupper(tt[, 1])), ]
# tt2$restrict[is.na(tt2$disease)] <- "unknown"
# tt2$restrict[is.na(tt2$restrict)] <- "unrestricted"
# suitability.table$embargoed <- tt2$restrict
# 
# write.csv(suitability.table, file="suitability.table.csv")

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

bimods <- c("PAAD rnaseq2", "COADREAD rnaseq2", "UCEC rnaseq2", "KICH rnaseq2", "TGCT rnaseq2", "ESCA rnaseq", "HNSC rnaseq2", "LGG rnaseq2", "THYM rnaseq2", "STAD rnaseq", "GBM rnaseq2")
pdf("TCGA_PairwisePearson_bimodal.pdf", width=7.5, height=9)
par(mfrow=c(4,3))
j <- 0
bimods2 <- sub(" .+", "", bimods)
st <- suitability.table[match(bimods2, suitability.table$cancertype), ]
for (i in match(bimods, names(cor.list))){
  j <- j+1
  hist(cor.list[[i]], breaks="FD", xlab="PCC", xlim=c(0, 1),
       main=paste(names(cor.list)[i], "\n", st$Study.Name[j]))
}
dev.off()

# remove restricted 
#cor.list <- cor.list[-unlist(sapply(toupper(tt[!is.na(tt$restrict),1]), grep, names(cor.list)))]
#cor.list <- cor.list[!grepl("mrna$", names(cor.list))]
names(cor.list) <- gsub(" rnaseq2","", names(cor.list))
names(cor.list) <- gsub(" rnaseq","", names(cor.list))
names(cor.list) <- suitability.table$Study.Name[match(names(cor.list), suitability.table$cancertype)]
#names(cor.list) <- sub("cervical squamous cell carcinoma and endocervical adenocarcinoma", "cervical squamous cell carcinoma \n endocervical adenocarcinoma", names(cor.list))

d.f <- do.call(rbind, lapply(names(cor.list), function(i) data.frame(Cancer=i, COR=cor.list[[i]])))
d.f$Cancer <- factor(d.f$Cancer, levels=names(cor.list)[order(sapply(cor.list, quantile, p=0.99),decreasing=TRUE)])

library(ggplot2)
pdf("tcgacor.pdf", width=12,height=6)
ggplot(d.f, aes(Cancer,COR))+geom_violin()+theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+xlab("")+ylab("Pairwise Pearson Correlation")
dev.off()

##bimodal: KICH (RNAseq2), maybe KIRC and KIRP, LUSC, PAAD, PCPG, THCA, GBM

##very high correlations normally: THCA, HNSC, LIHC, maybe KIRP, LAML, PCPG, PRAD, STAD, THYM, ESCA

##looks good: ACC, BLCA, BRCA, CESC, COAD, COADREAD, DLBC, LGG,
