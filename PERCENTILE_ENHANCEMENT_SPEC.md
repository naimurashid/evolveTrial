# evolveTrial Percentile Enhancement Specification

**STATUS: ✅ IMPLEMENTED** (2025-12-07)

## Overview

This document specifies how to add exact sample size percentile calculation to `run_simulation_pure()` with minimal computational overhead.

The implementation was completed following this specification. Key functions modified:
- `run_simulation_pure()` - Added `return_percentiles` and `percentile_probs` parameters
- `run_scenarios()` - Added pass-through support for percentile parameters

## Problem Statement

Currently, evolveTrial returns only aggregate statistics (mean EN, power, type I error). For statistical SAP sections, we need the **distribution** of sample sizes across simulations, specifically:
- Min, 25th percentile, Median, 75th percentile, 90th percentile, Max

## Design Decision: Accumulate Raw Values vs Online Quantile Estimation

### Option A: Accumulate Raw Values (Recommended)
Store all per-replicate sample sizes in a vector, compute percentiles at the end.

**Memory cost:** ~8 bytes × n_reps × n_arms = ~80KB for 10,000 reps × 1 arm

**Pros:**
- Exact percentiles
- Simple implementation
- Memory cost is trivial for typical use cases (10K-50K reps)
- Supports arbitrary percentile queries post-hoc

**Cons:**
- Memory scales with n_reps (but this is negligible)

### Option B: Online Quantile Estimation (P² Algorithm)
Estimate quantiles incrementally without storing all values.

**Memory cost:** O(1) per tracked percentile

**Pros:**
- Constant memory regardless of n_reps

**Cons:**
- Approximate, not exact
- More complex implementation
- Fixed percentiles must be chosen upfront
- P² can have accuracy issues for heavy-tailed distributions

### Recommendation
**Use Option A (accumulate raw values).** The memory overhead is negligible (<1MB even for 100K reps), and we get exact results with simple implementation.

## Implementation Plan

### 1. Add New Parameter to `run_simulation_pure()`

```r
run_simulation_pure <- function(
    num_simulations,
    # ... existing params ...
    return_percentiles = FALSE,  # NEW: Return sample size percentiles
    percentile_probs = c(0, 0.25, 0.5, 0.75, 0.9, 1.0),  # NEW: Which percentiles to compute
    # ...
)
```

### 2. Initialize Storage (only when return_percentiles = TRUE)

After the existing initialization (around line 430):

```r
# Percentile tracking (when enabled)
if (isTRUE(return_percentiles)) {
  # Pre-allocate vectors for all replicates
  all_final_n <- vector("list", length(arm_names))
  names(all_final_n) <- arm_names
  for (arm in arm_names) {
    all_final_n[[arm]] <- numeric(num_simulations)
  }
}
```

### 3. Modify `simulate_chunk()` to Return Raw Values

In `simulate_chunk()`, add collection of per-replicate values when enabled:

```r
simulate_chunk <- function(sim_indices, seed = NULL, tick = function() {}) {
  # ... existing code ...

  # NEW: Track per-replicate sample sizes for percentiles
  chunk_final_n_raw <- if (collect_raw) {
    lapply(arm_names, function(arm) numeric(length(sim_indices)))
  } else NULL
  if (!is.null(chunk_final_n_raw)) names(chunk_final_n_raw) <- arm_names

  for (i in seq_along(sim_indices)) {
    # ... existing simulation code ...

    # After computing final_n_vec (around line 626):
    chunk_sum_final_n <- chunk_sum_final_n + final_n_vec

    # NEW: Store raw values for percentile calculation
    if (!is.null(chunk_final_n_raw)) {
      for (arm in arm_names) {
        chunk_final_n_raw[[arm]][i] <- final_n_vec[arm]
      }
    }

    # ... rest of loop ...
  }

  # Return structure (around line 637):
  list(
    sum_final_n = chunk_sum_final_n,
    # ... existing fields ...
    final_n_raw = chunk_final_n_raw  # NEW (NULL if not collecting)
  )
}
```

### 4. Modify `aggregate_results()` to Combine Raw Vectors

```r
aggregate_results <- function(partials) {
  # ... existing aggregation ...

  # NEW: Combine raw vectors across chunks
  if (!is.null(partials[[1]]$final_n_raw)) {
    for (arm in arm_names) {
      all_final_n[[arm]] <<- unlist(lapply(partials, function(p) p$final_n_raw[[arm]]))
    }
  }
}
```

### 5. Compute and Return Percentiles

After `aggregate_results()` call (around line 775), add percentile computation:

```r
# Compute percentiles if requested
percentiles_result <- NULL
if (isTRUE(return_percentiles) && !is.null(all_final_n)) {
  percentiles_result <- lapply(arm_names, function(arm) {
    quantile(all_final_n[[arm]], probs = percentile_probs, na.rm = TRUE)
  })
  names(percentiles_result) <- arm_names
}

# Modify return (around line 798)
if (isTRUE(return_percentiles)) {
  return(list(
    summary = results_data,
    percentiles = list(
      N = percentiles_result,
      probs = percentile_probs
    )
  ))
} else {
  return(results_data)
}
```

### 6. Update `run_scenarios()` to Support Percentiles

In `run_scenarios()`, pass through the percentile parameters:

```r
run_scenarios <- function(base_args, scens, parallel = FALSE, seed = NULL,
                          return_percentiles = FALSE,
                          percentile_probs = c(0, 0.25, 0.5, 0.75, 0.9, 1.0)) {
  # ... existing code ...

  for (i in seq_along(scens)) {
    args_i <- modifyList(base_args, scens[[i]])
    args_i$return_percentiles <- return_percentiles
    args_i$percentile_probs <- percentile_probs

    res <- do.call(run_simulation_pure, args_i)

    if (return_percentiles) {
      all_results[[i]] <- res$summary
      all_percentiles[[i]] <- res$percentiles
    } else {
      all_results[[i]] <- res
    }
  }

  # Return
  if (return_percentiles) {
    list(
      summary = do.call(rbind, all_results),
      percentiles = all_percentiles
    )
  } else {
    do.call(rbind, all_results)
  }
}
```

## Usage Example

```r
# Run with percentile collection
result <- run_simulation_pure(
  num_simulations = 10000,
  # ... other params ...
  return_percentiles = TRUE,
  percentile_probs = c(0, 0.25, 0.5, 0.75, 0.9, 1.0)
)

# Access results
summary_df <- result$summary  # Same as before
percentiles <- result$percentiles$N  # List by arm

# For single-arm:
# percentiles$Experimental gives: c(min, 25th, median, 75th, 90th, max)
```

## Performance Impact

### Memory
- 10,000 reps × 1 arm: ~80 KB
- 10,000 reps × 3 arms: ~240 KB
- 100,000 reps × 3 arms: ~2.4 MB

This is negligible compared to the simulation state memory.

### CPU
- Storage: One assignment per replicate per arm (~negligible)
- Percentile computation: O(n log n) for sorting, done once at end
- For 10,000 reps: <10ms additional time

### Backward Compatibility
- When `return_percentiles = FALSE` (default), behavior is identical to current version
- Return type changes only when `return_percentiles = TRUE`

## Testing

Add tests to verify:
1. Percentiles are correct (compare to known distributions)
2. Means computed from raw values match Exp_N
3. Parallel and sequential modes give consistent results
4. Default behavior unchanged

## Files to Modify

1. `R/simulation_driver.R` - Main changes in `run_simulation_pure()`
2. `R/design_analysis.R` - Update `run_scenarios()` wrapper
3. `man/run_simulation_pure.Rd` - Document new parameters
4. `tests/testthat/test_percentiles.R` - New test file

## Estimated Effort

- Implementation: 2-3 hours
- Testing: 1-2 hours
- Documentation: 30 minutes

Total: ~4-5 hours
