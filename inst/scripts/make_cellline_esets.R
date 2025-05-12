datadir <- "."
load("mRNAexpression/NCI60.dataset.rda")
library(org.Hs.eg.db)
xx <- as.list(org.Hs.egALIAS2EG)
# Remove pathway identifiers that do not map to any entrez gene id
xx <- xx[!is.na(xx)]

rna <- rna[!is.na(rna$X.2),]
rna.entrez <- sapply(xx[rna$X.2], function(x) x[[1]])
idx <- !is.na(names(rna.entrez))
rna <- rna[idx,]
rna.entrez <- unlist(rna.entrez[idx])

rna.eset <- rna[!duplicated(rna.entrez), -(1:8)]
rownames(rna.eset) <- paste(rna.entrez[!duplicated(rna.entrez)], "at", sep="_")
nci60.eset <- ExpressionSet(as.matrix(rna.eset))

saveRDS(nci60.eset, file="mRNAexpression/NCI60_mrna.rds")

#ccle

expr.dat <- read.table(file.path(datadir,"mRNAexpression/CCLE_Expression_Entrez_2012-09-29.gct"), skip=2, sep="\t", as.is=TRUE,header = TRUE, check.names=FALSE)
cor(expr.dat[, colnames(expr.dat) %in% "NCIH292_LUNG"])

expr.names <- colnames(expr.dat)[-1:-2]
expr.matrix <- as.matrix(expr.dat[,-(1:2)])
rownames(expr.matrix) <- expr.dat$Name

ccle.eset <- ExpressionSet(expr.matrix)

saveRDS(ccle.eset, file="mRNAexpression/CCLE_mrna.rds")
