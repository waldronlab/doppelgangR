
#' Identify likely duplicate samples from genomic or meta-data
#'
#' The main function is doppelgangR(), which takes as minimal input a list of
#' ExpressionSet object, and searches all list pairs for duplicated samples.
#' The search is based on the genomic data (exprs(eset)), phenotype/clinical
#' data (pData(eset)), and "smoking guns" - supposedly unique identifiers found
#' in pData(eset).
#'
#' @importFrom Biobase exprs featureNames sampleNames pData "sampleNames<-"
#' "exprs<-" "pData<-"
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom sva ComBat
#' @importFrom digest digest
#' @importFrom impute impute.knn
#' @importFrom mnormt pmt biv.nt.prob pmnorm
#' @import methods
#' @importFrom grDevices dev.interactive
#' @importFrom graphics abline contour curve hist legend lines pairs par points title
#' @importFrom stats coef complete.cases cor cov2cor
#'           dchisq dnorm dt integrate lm.fit
#'           model.matrix na.omit nlminb optim pchisq
#'           pf pnorm pt qchisq qf qnorm rchisq
#'           resid rnorm runif uniroot var
#'           weighted.mean
#' @importFrom utils combn
#' @importFrom SummarizedExperiment assay colData
#'
#' @name doppelgangR-package
#' @aliases NULL doppelgangR-package
#' @docType package
#' @author Levi Waldron, Markus Riester, Marcel Ramos
#'
#' @keywords package
"_PACKAGE"
