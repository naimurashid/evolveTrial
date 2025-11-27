# Calibrate interim and final thresholds for single-arm designs

Sweeps candidate interim and final posterior probability thresholds and
returns the best combination achieving the desired type I error under
the null scenario(s).

## Usage

``` r
calibrate_alpha(
  base_args,
  scens_null,
  thr_grid_interim = c(0.9, 0.95, 0.975),
  thr_grid_final = c(0.95, 0.975, 0.99),
  sims = 300
)
```

## Arguments

- base_args:

  Baseline argument list passed to
  [`run_scenarios()`](run_scenarios.md).

- scens_null:

  Scenario list representing null hypotheses.

- thr_grid_interim:

  Numeric vector of interim success thresholds to try.

- thr_grid_final:

  Numeric vector of final success thresholds to try.

- sims:

  Number of simulations per candidate setting.

## Value

A list containing the chosen thresholds and corresponding estimated type
I error.
