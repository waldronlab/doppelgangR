## Changes in version 1.41.1

### Migration to future Parallel Framework
* Replaced `BiocParallel` package usage with the CRAN `future` and `future.apply` framework for parallel computing.
* Removed `BPPARAM` parameter from `doppelgangR()` to align with `future` package best practices (enabling global backend configuration via `future::plan()`).
* Added a dummy parameter catcher in `doppelgangR(...)` that throws a helpful deprecation warning if legacy `BPPARAM` is used.

### Performance and Memory Optimizations
* **Cartesian Product Vector Optimization**: Added a vector fast path in `.outer2df` to construct coordinates directly using vector replication (`rep`), yielding a **46x speedup** on 2,000 samples and eliminating massive memory spikes from legacy string splits.
* **Phenotype Distance Lookups**: Precalculated value frequencies once per column in `phenoDist()`, allowing `vectorWeightedDist()` to perform $O(1)$ lookups instead of scanning columns dynamically. This reduced complexity from cubic to quadratic and achieved a **2.3x speedup** on 500 samples.
* **Hamming Distance Helper Optimization**: Replaced `length(na.omit(z))` with `sum(!is.na(z))` inside `vectorHammingDist()` for an **8x speedup** per call.

### Bug Fixes & Compliance
* **Recycling Warning Fixed**: Resolved a vector length mismatch recycling warning in `doppelgangR()` intermediate pruning when smoking guns are enabled. Replaced the element-wise logical OR on vectors of unequal sizes with a robust, key-based merging and subsetting strategy.
* **BiocCheck & Code Style Compliance**:
  * Migrated unsafe `sapply()` calls to type-safe `vapply()`.
  * Replaced manual integer ranges (`1:...`) with `seq_len()`.
  * Removed `paste()` and `paste0()` from inside error and warning condition signals.
  * Replaced raw `cat` and `print` calls in trace logging blocks with `message()`.

## Changes in version 0.5.0

* Add knn imputation by default.

## Changes in version 0.1.1

* remove "smoking gun" columns from calculation of phenotype similarity.

## Changes in version 0.1.0

* initial commit
