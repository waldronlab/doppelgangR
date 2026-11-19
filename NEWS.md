## Changes in version 1.42.0

### Significant user-visible changes

* Parallel computation now uses the CRAN `future` framework instead of
  `BiocParallel`. The `BPPARAM` argument of `doppelgangR()` is removed;
  choose a backend with `future::plan()`. Supplying `BPPARAM` gives a
  deprecation warning (#29, @lwaldron).

### New features

* `doppelgangR()` accepts `SummarizedExperiment` objects, and lists of
  them, converting them to `ExpressionSet` internally (#30, @lwaldron).

### Bug fixes and minor improvements

* `doppelgangR()` is faster and uses less memory on large datasets
  (#30, @lwaldron).
* Fix a recycling warning when smoking guns are enabled (#29, @lwaldron).
* Fix examples, build the vignette when suggested packages are missing,
  and expand test coverage (#26, #27, #28, @lwaldron).

## Changes in version 0.5.0

* Add knn imputation by default.

## Changes in version 0.1.1

* remove "smoking gun" columns from calculation of phenotype similarity.

## Changes in version 0.1.0

* initial commit
