# Calibrate historical-control thresholds over a grid

Evaluates a grid of interim/final thresholds and superiority margins
under null and alternative scenarios, returning full operating
characteristics for each combination.

## Usage

``` r
grid_calibrate(
  base_args,
  null_med = 6,
  alt_med = 9,
  margins_abs = c(1, 2, 3),
  interim_thr_grid = c(0.9, 0.95),
  final_thr_grid = c(0.95, 0.975, 0.99),
  sims = 400,
  target_alpha = 0.1,
  seed = 123,
  parallel = TRUE
)
```

## Arguments

- base_args:

  Baseline argument list passed to
  [`run_scenarios()`](run_scenarios.md).

- null_med:

  Control-arm median under the null hypothesis.

- alt_med:

  Experimental median under the alternative.

- margins_abs:

  Numeric vector of absolute superiority margins (months).

- interim_thr_grid:

  Interim success probabilities to evaluate.

- final_thr_grid:

  Final posterior success probabilities to evaluate.

- sims:

  Number of simulations per grid point.

- target_alpha:

  Target type I error used when ranking feasible designs.

- seed:

  RNG seed.

- parallel:

  Logical; run scenarios in parallel.

## Value

A list with components `all`, `feasible`, and `top` summarising the grid
results.
