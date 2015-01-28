library(RUnit)
library(affy)
library(doppelgangR)

ncor <- 0; npheno <- 0; nsmoking <- 0
options(stringsAsFactors=FALSE)
set.seed(1)
m1 <- matrix(rnorm(1100), ncol=11)
colnames(m1) <- paste("m", 1:11, sep="")
rownames(m1) <- make.names(1:nrow(m1))
n1 <- matrix(rnorm(1000), ncol=10)
colnames(n1) <- paste("n", 1:10, sep="")
rownames(n1) <- make.names(1:nrow(n1))
##m:1 & n:1 are expression doppelgangers:
m1[, 1] <- n1[, 1] + rnorm(100, sd=0.25); ncor <- ncor+1
##m:2 & m:3 are expression doppelgangers:
m1[, 2] <- m1[, 3] + rnorm(100, sd=0.25); ncor <- ncor+1
##n:2 & n:3 are expression doppelgangers:
n1[, 2] <- n1[, 3] + rnorm(100, sd=0.25); ncor <- ncor+1
#n:4 & m:6 are expression doppelgangers:
n1[, 4] <- m1[, 6] + rnorm(100, sd=0.25); ncor <- ncor+1
#n:5 & m:4 are expression doppelgangers:
n1[, 5] <- m1[, 4] + rnorm(100, sd=0.25); ncor <- ncor+1

##
##m:10 and n:10 are phenotype doppelgangers:
m.pdata <- matrix(letters[sample(1:26, size=110, replace=TRUE)], ncol=10)
n.pdata <- matrix(letters[sample(1:26, size=100, replace=TRUE)], ncol=10)
n.pdata[10, ] <- m.pdata[10, ]; npheno <- npheno+1
##Create ExpressionSets
m.eset <- ExpressionSet(assayData=m1)
m.eset$id <- toupper(colnames(m1))
pData(m.eset) <- data.frame(c(pData(m.eset), data.frame(m.pdata)))
##
n.eset <- ExpressionSet(assayData=n1)
n.eset$id <- toupper(colnames(n1))
##m5 and n4 are "smoking gun" doppelgangers:
n.eset$id[4] <- "gotcha"
m.eset$id[5] <- "gotcha"
pData(n.eset) <- data.frame(c(pData(n.eset), data.frame(n.pdata)))
nsmoking <- nsmoking+1
##
esets <- list(m=m.eset, n=n.eset)

##------------------------------------------
##Check of all three types of doppelgangers:
##------------------------------------------
res1 <- doppelgangR(esets, manual.smokingguns="id", automatic.smokingguns=FALSE, cache.dir=NULL)
df1 <- summary(res1)

checkIdentical(df1[df1$sample1=="m:1" & df1$sample2=="n:1", "expr.doppel"], TRUE)
checkIdentical(df1[df1$sample1=="m:2" & df1$sample2=="m:3", "expr.doppel"], TRUE)
checkIdentical(df1[df1$sample1=="m:6" & df1$sample2=="n:4", "expr.doppel"], TRUE)
checkIdentical(df1[df1$sample1=="n:2" & df1$sample2=="n:3", "expr.doppel"], TRUE)
checkIdentical(df1[df1$sample1=="m:6" & df1$sample2=="n:4", "expr.doppel"], TRUE)
checkIdentical(df1[df1$sample1=="m:10" & df1$sample2=="n:10", "pheno.doppel"], TRUE)
checkIdentical(df1[df1$sample1=="m:5" & df1$sample2=="n:4", "smokinggun.doppel"], TRUE)
checkEquals(nrow(df1), ncor+npheno+nsmoking)
checkEquals(sum(df1$expr.doppel), ncor)
checkEquals(sum(df1$pheno.doppel), npheno)
checkEquals(sum(df1$smokinggun.doppel), nsmoking)
checkEquals(sum(is.na(df1$expr.similarity)), 0)
checkEquals(sum(is.na(df1$pheno.similarity)), 0)
checkEquals(sum(is.na(df1$smokinggun.similarity)), 0)
checkEquals(df1$id, c("M2:M3", "N2:N3", "M1:N1", "gotcha:gotcha", "M6:gotcha", "M4:N5", "M10:N10"))
for (i in match(paste("X", 1:10, sep=""), colnames(df1))){
    cat(paste("Checking column", i, "\n"))
    checkEquals(all(grepl("[a-z]:[a-z]", df1[[i]])), TRUE)
}

##------------------------------------------
cat("\n")
cat("Check without smoking guns: \n")
##------------------------------------------
res2 <- doppelgangR(esets, smokingGunFinder.args=NULL, cache.dir=NULL)
df2 <- summary(res2)
for (i in grep("pheno.similarity|smokinggun.similarity", colnames(df1), invert=TRUE)){
    cat(paste("Checking column", i, "\n"))
    checkEquals(df2[, i], df1[!df1$smokinggun.doppel, i])
}


##------------------------------------------
cat("\n")
cat("Check without phenotype: \n")
##------------------------------------------
res3 <- doppelgangR(esets, phenoFinder.args=NULL, manual.smokingguns="id", automatic.smokingguns=FALSE, cache.dir=NULL)
df3 <- summary(res3)
for (i in grep("pheno.similarity", colnames(df1), invert=TRUE)){
    cat(paste("Checking column", i, "\n"))
    checkEquals(df3[, i], df1[!df1$pheno.doppel, i])
}

##------------------------------------------
cat("\n")
cat("Check without expression: \n")
##------------------------------------------
res4 <- doppelgangR(esets, corFinder.args=NULL, manual.smokingguns="id", automatic.smokingguns=FALSE, cache.dir=NULL)
df4 <- summary(res4)
for (i in grep("expr.similarity", colnames(df1), invert=TRUE)){
    cat(paste("Checking column", i, "\n"))
    checkEquals(df4[, i], df1[!df1$expr.doppel, i])
}


##------------------------------------------
cat("\n")
cat("Check pruning: \n")
##------------------------------------------
res5 <- doppelgangR(esets, manual.smokingguns="id", automatic.smokingguns=FALSE, intermediate.pruning=TRUE, cache.dir=NULL)
df5 <- summary(res5)
checkEquals(df1, df5)


##------------------------------------------
cat("\n")
cat("Check corFinder function: \n")
##------------------------------------------
cor1 <- corFinder(eset.pair=esets)
cor2 <- corFinder(eset.pair=esets[c(2, 1)])
checkEquals(cor1, t(cor2))

cor1 <- corFinder(eset.pair=esets, use.ComBat=FALSE)
cor2 <- corFinder(eset.pair=esets[c(2, 1)], use.ComBat=FALSE)
checkEquals(cor1, t(cor2))

cor1 <- corFinder(eset.pair=esets[c(1, 1)])
cor2 <- corFinder(eset.pair=esets[c(1, 1)], use.ComBat=FALSE)
checkEquals(cor1, cor2)
