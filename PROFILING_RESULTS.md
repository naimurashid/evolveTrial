# evolveTrial Profiling Results

**Date:** 2025-12-03 **Profiler:** R Rprof (sequential execution for
visibility)

## Executive Summary

Profiling reveals that the main bottleneck in evolveTrial is **PSOCK
cluster spawn/teardown overhead**, not the simulation logic itself. Each
call to [`run_simulation_pure()`](reference/run_simulation_pure.md) with
`parallel_replicates = TRUE` creates and destroys a parallel cluster,
which adds significant overhead when called repeatedly (e.g., in
Bayesian optimization loops).

## Timing Breakdown

| Component                        | Time      | Notes                |
|----------------------------------|-----------|----------------------|
| Single simulator call (200 reps) | 40.59 sec | ~202 ms/rep          |
| Sequential same config           | 0.06 sec  | ~1.2 ms/rep          |
| Overhead ratio                   | ~170x     | Due to PSOCK cluster |

## Identified Bottlenecks

### 1. PSOCK Cluster Spawn/Teardown (CRITICAL)

**Location:** `R/simulation_driver.R` lines 689-705

``` r
cl <- parallel::makeCluster(workers, type = cluster_type)
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterCall(cl, function(pkg) {
  suppressPackageStartupMessages(require(pkg, character.only = TRUE))
}, pkg_name)
```

**Impact:** - Each BO evaluation creates/destroys a cluster - Worker
processes must load evolveTrial package each time - PSOCK requires full
R process spawn + socket communication

### 2. JIT Compilation Overhead

**Location:** Throughout R code

**Profile evidence:** - `h`, `tryInline`, `cmpCall`, `cmp`, `genCode`,
`cmpfun` in top functions - `findCenvVar`, `cb$makelabel`, `putconst`
show compiler activity

**Impact:** - Functions recompiled in each worker process - No caching
of compiled bytecode across calls

### 3. Invisible Worker Time

**Issue:** Rprof runs in parent process; PSOCK workers are separate R
processes

**Result:** - 40.59 sec total, but only 0.08 sec visible to profiler -
Cannot directly profile what happens inside `simulate_chunk()`

## Optimization Recommendations

### High Priority (Large Impact)

#### 1. Implement Cluster Pooling

Add option to reuse a persistent cluster across multiple BO evaluations:

``` r
# New API option
run_simulation_pure(..., cluster = NULL)  # If NULL, use internal pooling

# Usage pattern
cl <- evolveTrial::get_cluster(workers = 8)  # Create once
for (i in 1:100) {
  run_simulation_pure(..., cluster = cl)  # Reuse
}
evolveTrial::release_cluster(cl)  # Cleanup
```

**Expected improvement:** 10-50x faster for repeated calls

#### 2. Use FORK Instead of PSOCK on Linux

FORK clusters share memory and don’t require serialization:

``` r
cluster_type <- if (.Platform$OS.type == "unix") "FORK" else "PSOCK"
```

**Benefits:** - No package loading in workers (inherited) - No data
serialization (copy-on-write) - Much faster startup

**Expected improvement:** 2-5x on Linux/macOS

#### 3. Threshold for Sequential Execution

For small rep counts, sequential is faster than parallel overhead:

``` r
# Rule of thumb: parallel_replicates makes sense only if:
# num_simulations * time_per_rep > cluster_overhead (typically 5-10 sec)
if (num_simulations < 50) {
  parallel_replicates <- FALSE
}
```

### Medium Priority

#### 4. Pre-compile Critical Functions

Add to package `.onLoad()`:

``` r
.onLoad <- function(libname, pkgname) {
  # Pre-compile hot path functions
  assign("simulate_chunk", compiler::cmpfun(simulate_chunk), envir = parent.env(environment()))
  assign("interim_check", compiler::cmpfun(interim_check), envir = parent.env(environment()))
  # etc.
}
```

#### 5. Export Cluster for External Use

Allow BO packages to manage cluster lifecycle:

``` r
#' Create a reusable evolveTrial cluster
#' @export
create_simulation_cluster <- function(workers = NULL, cluster_type = "auto") {
  workers <- workers %||% max(1L, parallel::detectCores() - 1L)
  if (cluster_type == "auto") {
    cluster_type <- if (.Platform$OS.type == "unix") "FORK" else "PSOCK"
  }
  cl <- parallel::makeCluster(workers, type = cluster_type)
  parallel::clusterEvalQ(cl, library(evolveTrial))
  class(cl) <- c("evolveTrial_cluster", class(cl))
  cl
}
```

### Low Priority

#### 6. Internal Simulation Optimization

The `simulate_chunk()` function appears efficient (~1.2 ms/rep when run
sequentially). C++ implementations for posterior sampling are already in
place.

Potential micro-optimizations: - Vectorize patient enrollment loop -
Pre-allocate result matrices - Use integer arithmetic where possible

------------------------------------------------------------------------

## Implemented Optimizations (2025-12-03)

All three high-priority optimizations have been implemented in
`R/simulation_driver.R`:

### 1. Cluster Pooling ✅

New parameter `cluster` added to
[`run_simulation_pure()`](reference/run_simulation_pure.md) to accept an
external cluster:

``` r
# Create cluster once for all BO evaluations
cl <- evolveTrial::create_simulation_cluster(workers = 8)

# Use in repeated BO evaluations (100+ calls)
for (i in 1:100) {
  result <- run_simulation_pure(
    num_simulations = 500,
    ...,
    parallel_replicates = TRUE,
    cluster = cl  # Reuses cluster, no spawn/teardown overhead
  )
}

# Clean up when done
evolveTrial::release_cluster(cl)
```

**New exported functions:** -
`create_simulation_cluster(workers, cluster_type)` - Create reusable
cluster - `release_cluster(cluster)` - Clean up cluster resources

### 2. FORK Auto-Detection ✅

Default `cluster_type` changed from `"PSOCK"` to `"auto"`:

``` r
# Auto-detection logic (in simulation_driver.R)
if (cluster_type == "auto") {
  cluster_type <- if (.Platform$OS.type == "unix") "FORK" else "PSOCK"
}
```

**Benefits on Linux/macOS:** - FORK inherits parent environment (no
package loading in workers) - Copy-on-write memory (no serialization
overhead) - ~2-5x faster cluster startup

### 3. Sequential Threshold ✅

Increased threshold from `workers * 2` to `100` (configurable):

``` r
# Threshold is configurable via option
seq_threshold <- getOption("evolveTrial.sequential_threshold", 100L)

if (num_simulations < seq_threshold) {
  # Run sequentially - faster than parallel overhead
  chunk_results <- list(simulate_chunk(...))
}
```

**Rationale:** - Cluster spawn/teardown: ~5-10 sec - Sequential
simulation: ~1-2 ms/rep - Breakeven at ~100 reps with typical config

**Configuration:**

``` r
# Lower threshold for fast machines
options(evolveTrial.sequential_threshold = 50)

# Higher threshold if cluster overhead is high
options(evolveTrial.sequential_threshold = 200)
```

### Expected Performance Improvements

| Scenario                         | Before    | After    | Improvement |
|----------------------------------|-----------|----------|-------------|
| BO with 100 evals, 200 reps each | ~4000 sec | ~400 sec | **10x**     |
| Single call, 50 reps             | ~10 sec   | ~0.1 sec | **100x**    |
| Single call, 500 reps (Linux)    | ~40 sec   | ~15 sec  | **2-3x**    |

### API Changes Summary

**New parameters in
[`run_simulation_pure()`](reference/run_simulation_pure.md):** -
`cluster_type = c("auto", "PSOCK", "FORK")` - Default changed to
“auto” - `cluster = NULL` - Pass external cluster for pooling

**New exported functions:** -
[`create_simulation_cluster()`](reference/create_simulation_cluster.md) -
Create persistent cluster -
[`release_cluster()`](reference/release_cluster.md) - Clean up cluster

**New option:** - `evolveTrial.sequential_threshold` - Override default
threshold (100)

## Files Created

- `tests/profile_sequential.R` - Sequential profiling script
- `PROFILING_RESULTS.md` - This document

## Verification

To reproduce profiling:

``` bash
cd /home/naimrashid/Downloads/evolveTrial
Rscript tests/profile_sequential.R
```

To compare with parallel (from paper directory):

``` bash
cd /home/naimrashid/Downloads/adaptive-trial-bo-paper
Rscript scripts/profile_comprehensive.R
```
