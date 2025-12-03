# Pre-compute control arm posteriors for caching in multi-arm comparisons

PERFORMANCE: When comparing multiple experimental arms against the same
control, this function computes control posteriors once to avoid
redundant computation.

## Usage

``` r
precompute_ctrl_posteriors(slCtrl, args, num_samples)
```

## Arguments

- slCtrl:

  Control arm slice from slice_arm_data_at_time()

- args:

  Trial arguments containing prior parameters

- num_samples:

  Number of posterior samples

## Value

List with lamC (hazard samples matrix) and medCtrl (median vector)
