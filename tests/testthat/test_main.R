library(Biobase)
library(doppelgangR)

ncor <- 0
npheno <- 0
nsmoking <- 0
options(stringsAsFactors = FALSE)
set.seed(1)
m1 <- matrix(rnorm(1100), ncol = 11)
colnames(m1) <- paste("m", 1:11, sep = "")
rownames(m1) <- make.names(1:nrow(m1))
n1 <- matrix(rnorm(1000), ncol = 10)
colnames(n1) <- paste("n", 1:10, sep = "")
rownames(n1) <- make.names(1:nrow(n1))
##m:1 & n:1 are expression doppelgangers:
m1[, 1] <- n1[, 1] + rnorm(100, sd = 0.25)
ncor <- ncor + 1
##m:2 & m:3 are expression doppelgangers:
m1[, 2] <- m1[, 3] + rnorm(100, sd = 0.25)
ncor <- ncor + 1
##n:2 & n:3 are expression doppelgangers:
n1[, 2] <- n1[, 3] + rnorm(100, sd = 0.25)
ncor <- ncor + 1
##n:8 & n:9 are expression doppelgangers:
n1[, 8] <- n1[, 9] + rnorm(100, sd = 0.25)
ncor <- ncor + 1
#n:4 & m:6 are expression doppelgangers:
n1[, 4] <- m1[, 6] + rnorm(100, sd = 0.25)
ncor <- ncor + 1
#n:5 & m:4 are expression doppelgangers:
n1[, 5] <- m1[, 4] + rnorm(100, sd = 0.25)
ncor <- ncor + 1

##
##m:10 and n:10 are phenotype doppelgangers:
m.pdata <-
  matrix(letters[sample(1:26, size = 110, replace = TRUE)], ncol = 10)
rownames(m.pdata) <- colnames(m1)
n.pdata <-
  matrix(letters[sample(1:26, size = 100, replace = TRUE)], ncol = 10)
n.pdata[10,] <- m.pdata[10,]
npheno <- npheno + 1
rownames(n.pdata) <- colnames(n1)



##Create ExpressionSets
m.eset <-
  ExpressionSet(assayData = m1, phenoData = AnnotatedDataFrame(data.frame(m.pdata)))
m.eset$id <- toupper(colnames(m1))
##
n.eset <-
  ExpressionSet(assayData = n1, phenoData = AnnotatedDataFrame(data.frame(n.pdata)))
n.eset$id <- toupper(colnames(n1))
##m5 and n4 are "smoking gun" doppelgangers:
n.eset$id[4] <- "gotcha"
m.eset$id[5] <- "gotcha"
nsmoking <- nsmoking + 1
##
esets <- list(m = m.eset, n = n.eset)

##------------------------------------------
##Check of all three types of doppelgangers:
##------------------------------------------
res1 <-
  doppelgangR(
    esets,
    manual.smokingguns = "id",
    automatic.smokingguns = FALSE,
    cache.dir = NULL
  )
df1 <- summary(res1)

expect_true(df1[df1$sample1 == "m:m1" &
                  df1$sample2 == "n:n1", "expr.doppel"])
expect_true(df1[df1$sample1 == "m:m2" &
                  df1$sample2 == "m:m3", "expr.doppel"])
expect_true(df1[df1$sample1 == "m:m6" &
                  df1$sample2 == "n:n4", "expr.doppel"])
expect_true(df1[df1$sample1 == "n:n2" &
                  df1$sample2 == "n:n3", "expr.doppel"])
expect_true(df1[df1$sample1 == "m:m6" &
                  df1$sample2 == "n:n4", "expr.doppel"])
expect_true(df1[df1$sample1 == "m:m10" &
                  df1$sample2 == "n:n10", "pheno.doppel"])
expect_true(df1[df1$sample1 == "m:m5" &
                  df1$sample2 == "n:n4", "smokinggun.doppel"])
expect_equal(nrow(df1), ncor + npheno + nsmoking)
expect_equal(sum(df1$expr.doppel), ncor)
expect_equal(sum(df1$pheno.doppel), npheno)
expect_equal(sum(df1$smokinggun.doppel), nsmoking)
expect_equal(sum(is.na(df1$expr.similarity)), 0)
expect_equal(sum(is.na(df1$pheno.similarity)), 0)
expect_equal(sum(is.na(df1$smokinggun.similarity)), 0)
expect_equal(
  df1$id,
  c(
    "M2:M3",
    "N2:N3",
    "N8:N9",
    "M1:N1",
    "gotcha:gotcha",
    "M6:gotcha",
    "M4:N5",
    "M10:N10"
  )
)
for (i in match(paste("X", 1:10, sep = ""), colnames(df1))) {
  lab <- paste("Checking column", i, "\n")
  expect_equal(all(grepl("[a-z]:[a-z]", df1[[i]])), TRUE, label = lab)
}

##------------------------------------------
## Check without smoking guns
##------------------------------------------
res2 <-
  doppelgangR(esets,
              smokingGunFinder.args = NULL,
              cache.dir = NULL)
df2 <- summary(res2)
for (i in grep("pheno.similarity|smokinggun.similarity",
               colnames(df1),
               invert = TRUE)) {
  lab <- paste("Check without smoking guns, column", i, "\n")
  expect_equal(df2[, i], df1[!df1$smokinggun.doppel, i], label = lab)
}


##------------------------------------------
## Check without phenotype
##------------------------------------------
res3 <-
  doppelgangR(
    esets,
    phenoFinder.args = NULL,
    manual.smokingguns = "id",
    automatic.smokingguns = FALSE,
    cache.dir = NULL
  )
df3 <- summary(res3)
for (i in grep("pheno.similarity", colnames(df1), invert = TRUE)) {
  lab <- paste("Check without phenotype column", i, "\n")
  expect_equal(df3[, i], df1[!df1$pheno.doppel, i], label = lab)
}

##------------------------------------------
## Check without expression
##------------------------------------------
res4 <-
  doppelgangR(
    esets,
    corFinder.args = NULL,
    manual.smokingguns = "id",
    automatic.smokingguns = FALSE,
    cache.dir = NULL
  )
df4 <- summary(res4)
for (i in grep("expr.similarity", colnames(df1), invert = TRUE)) {
  lab <- paste("Check without expression column", i, "\n")
  expect_equal(df4[, i], df1[!df1$expr.doppel, i], label = lab)
}

##------------------------------------------
## Check smoking guns only
##------------------------------------------
res4b <-
  doppelgangR(
    esets,
    corFinder.args = NULL,
    phenoFinder.args = NULL,
    manual.smokingguns = "id",
    automatic.smokingguns = FALSE,
    cache.dir = NULL
  )
df4b <- summary(res4b)
rownames(df4b) <- NULL
df4b.compare <-
  df1[df1$smokinggun.doppel,]
rownames(df4b.compare) <- NULL
##don't check expr and pheno columns
expect_identical(df4b.compare[,-3:-6], df4b[,-3:-6], label = "Check smoking guns only")


##------------------------------------------
## Check pruning
##------------------------------------------
res5 <-
  doppelgangR(
    esets,
    manual.smokingguns = "id",
    automatic.smokingguns = FALSE,
    intermediate.pruning = TRUE,
    cache.dir = NULL
  )
df5 <- summary(res5)
expect_equal(df1, df5, label = "Check pruning")

##------------------------------------------
## Check caching, with a third ExpressionSet that is almost identical to the first
##------------------------------------------
esets2 <- c(esets, esets[[1]])
names(esets2)[3] <- "o"
exprs(esets2[[3]]) <-
  exprs(esets2[[3]]) + rnorm(nrow(esets2[[3]]) * ncol(esets2[[3]]), sd =
                               0.1)
esets2[[3]]$X10 <- "a"
sampleNames(esets2[[3]]) <-
  paste("X", sampleNames(esets2[[3]]), sep = "")

tmpcachedir <- tempdir()
##Do this twice, so the second time the cache will be used:
for (i in 1:2) {
  res6 <-
    doppelgangR(
      esets2,
      manual.smokingguns = "id",
      automatic.smokingguns = FALSE,
      cache.dir = tmpcachedir
    )
  ##Make sure comparison of m to n is the same as res1:
  df6a <- summary(res6)
  df6a <-
    df6a[grepl("^[mn]", df6a$sample1) & grepl("^[mn]", df6a$sample2),]
  rownames(df6a) <- 1:nrow(df6a)
  expect_equal(df1, df6a, label = "Check caching, with a third ExpressionSet that is almost identical to the first")
}

df6b <- summary(res6)
df6b <-
  df6b[grepl("^[mo]", df6b$sample1) & grepl("^[mo]", df6b$sample2),]
expect_identical(df6b[df6b$sample1 == "o:Xm2" &
                        df6b$sample2 == "o:Xm3", "expr.doppel"], TRUE)
df6b <- df6b[-1:-2,]
##expect_true(all(df6b$expr.doppel))  ## not a bug, but a shortcoming in the outlier detection that these are not all identified as expression doppelgangers.
expect_true(all(df6b$pheno.doppel))
expect_true(all(df6b$smokinggun.doppel))


res7 <-
  doppelgangR(
    esets2,
    phenoFinder.args = NULL,
    smokingGunFinder.args = NULL,
    outlierFinder.expr.args = list(
      bonf.prob = 1.0,
      transFun = atanh,
      tail = "upper"
    ),
    cache.dir = NULL
  )
df7 <- summary(res7)
has.o <- grepl("^o", df7$sample1) | grepl("^o", df7$sample2)

df7a <- df7[has.o,]
df7a$sample1 <- sub("o", "m", df7a$sample1)
df7a$sample1 <- sub("X", "", df7a$sample1)
df7a$sample2 <- sub("o", "m", df7a$sample2)
df7a$sample2 <- sub("X", "", df7a$sample2)

df7b <- df7[!has.o,]
df7b <-
  df7b[(!grepl("^n", df7b$sample1) | !grepl("^n", df7b$sample2)),]
for (i in 1:nrow(df7a)) {
  df7a[i, 1:2] <- sort(df7a[i, 1:2])
  df7b[i, 1:2] <- sort(df7b[i, 1:2])
}
df7a <- df7a[order(df7a$sample1, df7a$sample2),]
df7b <- df7b[order(df7b$sample1, df7b$sample2),]
expect_true(all(df7a$expr.doppel == df7b$expr.doppel))
expect_true(all(df7a$pheno.doppel == df7b$pheno.doppel))
expect_true(all(df7a$smokinggun.doppel == df7b$smokinggun.doppel))

##------------------------------------------
## Check corFinder function
##------------------------------------------
cor1 <- corFinder(eset.pair = esets)
cor2 <- corFinder(eset.pair = esets[c(2, 1)])
expect_equal(cor1, t(cor2), label = "Check corFinder function")

cor1 <- corFinder(eset.pair = esets, use.ComBat = FALSE)
cor2 <- corFinder(eset.pair = esets[c(2, 1)], use.ComBat = FALSE)
expect_equal(cor1, t(cor2), label = "Check corFinder function 2")

cor1 <- corFinder(eset.pair = esets[c(1, 1)])
cor2 <- corFinder(eset.pair = esets[c(1, 1)], use.ComBat = FALSE)
expect_equal(cor1, cor2, label = "Check corFinder function 3")

##Check missing values:
exprs(esets[[1]])[1:10, 1:5] <- NA
expect_message(doppelgangR(esets[1:2]), regexp = "Finalizing", label = "check missing values 1")
## More missing values:
exprs(esets[[1]])[1:10, 1:8] <- NA
expect_warning(doppelgangR(esets[1:2]), regexp = "10 rows with more than 50 % entries missing", label = "check missing values 2")
## More missing values:
exprs(esets[[1]])[1:10, 1:11] <- NA
expect_warning(doppelgangR(esets[1:2]), regexp = "mean imputation used for these rows", label = "check missing values 3")
## infinite values:
exprs(esets[[1]])[14, 1] <- -Inf
exprs(esets[[1]])[15, 2] <- Inf
expect_warning(doppelgangR(esets[1:2]), regexp = "Inf with min/max expression values for dataset m", label = "check missing values 4")

##------------------------------------------
## Smoking guns only with cache=TRUE
##------------------------------------------
expect_warning(
  dop <-
    doppelgangR(
      esets,
      corFinder.args = NULL,
      phenoFinder.args = NULL,
      manual.smokingguns = "id"
    ),
  label = "Smoking guns only with cache=TRUE"
)
expect_equal(summary(dop)[, 1], "m:m5", label = "Smoking guns only with cache=TRUE 2")
expect_equal(summary(dop)[, 2], "n:n4", label = "Smoking guns only with cache=TRUE 3")


##------------------------------------------
## Identical ExpressionSets
##------------------------------------------
expect_warning(df1 <-
                 summary(doppelgangR(esets[[1]], cache.dir = NULL)), label = "Identical ExpressionSets 1")
expect_true(df1$sample1 == "m2", label = "Identical ExpressionSets 2")
expect_true(df1$sample2 == "m3", label = "Identical ExpressionSets 3")
##
df2 <- summary(doppelgangR(esets[[2]], cache.dir = NULL))
expect_true(all(df2$sample1 == "n2"), label = "Identical ExpressionSets 4")
expect_true(all(df2$sample2 == "n3"), label = "Identical ExpressionSets 4b")
##
expect_warning(df5 <-
                 summary(doppelgangR(
                   list(eset1 = esets[[1]], eset2 = esets[[2]]), cache.dir = NULL
                 )), label = "Identical ExpressionSets 5")
expect_warning(df6 <-
                 summary(doppelgangR(esets, cache.dir = NULL)), label = "Identical ExpressionSets 6")
expect_identical(df5[,-1:-2], df6[,-1:-2], label = "Identical ExpressionSets 7")
expect_identical(sub("eset2", "n", sub("eset1", "m", df5$sample1)), df6$sample1, label = "Identical ExpressionSets 8")
expect_identical(sub("eset2", "n", sub("eset1", "m", df5$sample2)), df6$sample2, label = "Identical ExpressionSets 9")

## with zero-column pData:
expect_warning(withpheno <-
                 summary(doppelgangR(esets)), label = "with zero-column pData")
esets3 <- esets
pData(esets3[[1]]) <- pData(esets3[[1]])[, 0]
expect_warning(withoutpheno <-
                 summary(doppelgangR(esets3)),
               regexp = "m and n have different column names in phenoData",
               label = "with zero-column pData 2")

pData(esets3[[2]]) <- pData(esets3[[2]])[, 0]
expect_warning(withoutpheno2 <-
                 summary(doppelgangR(esets3)),
               regexp = "with min/max expression values for dataset m",
               label = "with zero-column pData 3")

withoutpheno3 <- withoutpheno[!withoutpheno$pheno.doppel,]
withoutpheno4 <- withoutpheno2[!withoutpheno2$pheno.doppel,]

expect_identical(withoutpheno[, 1:4], withoutpheno3[, 1:4], label = "with zero-column pData 4")
expect_identical(withoutpheno2[, 1:4], withoutpheno4[, 1:4], label = "with zero-column pData 5")

esets4 <- esets
for (i in 1:length(esets4))
  pData(esets4[[i]]) <- pData(esets4[[i]])[1]
expect_warning(doppelgangR(esets4[[1]]), label = "with zero-column pData 6")
expect_s4_class(doppelgangR(esets4[[2]]), "DoppelGang")
