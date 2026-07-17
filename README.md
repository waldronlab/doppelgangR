
# doppelgängR

## Introduction

*[doppelgangR](https://bioconductor.org/packages/3.23/doppelgangR)* is a
package for identifying duplicate samples within or between datasets of
transcriptome profiles. It is intended for microarray and RNA-seq gene
expression profiles where biological replicates are ordinarily more
distinct than technical replicates, as is the case for cancer types with
“noisy” genomes. It is intended for cases where per-gene summaries are
available but full genotypes are not, which is typical of public
databases such as the Gene Expression Omnibus.

Results from running `doppelgangR` on CRC, bladder, and ovarian [on
Dropbox](httpss://www.dropbox.com/sh/or3l7oo96x1ilm8/AAD1MVQVnoS6H7OExoYwOovSa?dl=0).

For the manuscript vignette, visit
<https://github.com/waldronlab/doppelgangR_paper>.

The `doppelgangR()` function identifies duplicates in three different
ways:

- **“expression”** doppelgängers have highly similar expression
  profiles, which are identified by default by having higher Pearson
  correlation than expected based on an empirical distribution of
  Pearson correlations between biological replicates. The type of
  correlation, and default use of ComBat batch correction, can be
  changed using the “corFinder.args” argument.

- **“phenotype”** doppelgängers have highly similar clinical or
  phenotype data, as contained in the phenoData slot of the
  `ExpressionSet`. In order to identify duplicates this way, it is
  required to curate the phenoData of each ExpressionSet they have
  identical column names, and encode phenotypes in the same way. For
  example, if each dataset provides information on age, this column of
  the phenoData could be called “age” in every dataset, and encoded as
  an integer number of years. If the phenoData slots are NULL then this
  type of checking will automatically be turned off. If they are not
  NULL but are also not curated, you should turn off phenotype checking
  by setting `phenoFinder.args=NULL`.

- **“smoking gun”** doppelgängers have the same value for an identifier
  that should be unique. You can enable this type of check by setting
  the argument “manual.smokingguns” to the names of columns containing
  supposedly unique identifiers, or setting “automatic.smokingguns” to
  TRUE, and the function will assume any column containing unique values
  within the column should also be unique across datasets.

This vignette focuses on the “expression” type of doppelgänger.

# Data types

Identification of doppelgängers is effective for both microarray and
**log-transformed** RNA-seq data, and even for matching samples that
have been profiled by microarray and RNA-seq.

# Case Study: Batch correction in Japanese datasets

We load for datasets by Yoshihara **et al.** that have been curated in
*[curatedOvarianData](https://bioconductor.org/packages/3.23/curatedOvarianData)*.
These are objects of class `ExpressionSet`.

``` r
library(curatedOvarianData)
data(GSE32062.GPL6480_eset)
data(GSE17260_eset)
```

The `doppelgangR` function requires a list of `ExpressionSet` objects as
input, which we create here:

``` r
testesets <- list(JapaneseA=GSE32062.GPL6480_eset,
                  Yoshihara2010=GSE17260_eset)
```

Now run `doppelgangR` with default arguments, except for setting
`phenoFinder.args=NULL`, which turns off checking for similar clinical
data in the `phenoData` slot of the ExpressionSet objects:

``` r
results1 <- doppelgangR(testesets, phenoFinder.args=NULL)
```

This creates an object of class `DoppelGang`, which has print, summary,
and plot methods. Summary method output not shown here due to voluminous
output:

``` r
summary(results1)
```

Plot creates a histogram of sample pairwise correlations within and
between each study:

``` r
par(mfrow=c(2,2), las=1)
plot(results1)
```

<div class="figure">

<img src="https://raw.githubusercontent.com/waldronlab/figures/refs/heads/devel/plotdop-1.png?token=GHSAT0AAAAAADOAOGYE6J2QRP54HL4XFORY2I3VLIA" alt="Doppelgängers identified on the basis of similar expression profiles. The vertical red lines indicate samples that were flagged." width="100%" />
<p class="caption">

Doppelgängers identified on the basis of similar expression profiles.
The vertical red lines indicate samples that were flagged.
</p>

</div>

One of these histograms can be drawn using the plot.pair argument:

``` r
plot(results1, plot.pair=c("JapaneseA", "JapaneseA"))
```

<div class="figure">

<img src="https://raw.githubusercontent.com/waldronlab/figures/refs/heads/devel/plotdop2-1.png?token=GHSAT0AAAAAADOAOGYFTY77UUIBHH744C7A2I3VNFA" alt="Pair plot of JapaneseA:JapaneseA Doppelgängers identified. The vertical red lines indicate samples that were flagged." width="100%" />
<p class="caption">

Pair plot of JapaneseA:JapaneseA Doppelgängers identified. The vertical
red lines indicate samples that were flagged.
</p>

</div>

# Important options

## Changing sensitivity

If after inspecting the histograms, you see that some visible outliers
were not caught, or non-outliers exceeded the sensitivity threshold, you
can change the default sensitivity using the argument:

``` r
outlierFinder.expr.args =
    list(bonf.prob = 0.5, transFun = atanh, tail = "upper")
```

The default 0.5 is a reasonable but arbitrary trade-off between
sensitivity and specificity which we have found to often select dataset
pairs containing duplicates, but to often not find *all* the duplicate
samples. Sensitivity can be increased by changing the bonf.prob
argument, *i.e.*:

``` r
results1 <- doppelgangR(testesets, 
        outlierFinder.expr.args = list(bonf.prob = 1.0, transFun = atanh, 
                                       tail = "upper"))
```

## Use of the ExpressionSet

The `doppelgangR()` function takes as its main argument a list of
`ExpressionSet` objects. If you just have matrices, you can easily
convert these to the `ExpressionSet` objects, for example:

``` r
mat <- matrix(1:4, ncol=2)
library(Biobase)
eset <- ExpressionSet(mat)
class(eset)
#> [1] "ExpressionSet"
#> attr(,"package")
#> [1] "Biobase"
```

## Parallelizing

The `doppelgangR()` function checks all pairwise combinations of
datasets in a list of `ExpressionSet` objects, and these dataset pairs
can be checked in parallel using multiple processing cores. This
functionality is imported from the "future" package. Please see
the "future" package documentation to configure your parallel backend.

``` r
library(future)
plan(multisession, workers = 8)
results2 <- doppelgangR(testesets)
```

## Caching

By default, the `doppelgangR()` function caches intermediate results to
make re-running with different arguments faster. Turn caching off by
setting the argument `cache.dir=NULL`.
