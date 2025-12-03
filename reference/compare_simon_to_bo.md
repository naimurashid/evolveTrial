# Validate BO calibration against Simon enumeration

Runs BO calibration on binary simulator and compares to exact Simon
design.

## Usage

``` r
compare_simon_to_bo(
  p0,
  p1,
  alpha = 0.1,
  beta = 0.2,
  n_max_search = 100,
  bo_fit = NULL,
  num_sims = 10000
)
```

## Arguments

- p0:

  Null response rate

- p1:

  Alternative response rate

- alpha:

  Type I error constraint

- beta:

  Type II error constraint (1 - power)

- n_max_search:

  Maximum N for Simon search

- bo_fit:

  Optional: pre-computed BO fit object

- num_sims:

  Number of simulations for Monte Carlo validation

## Value

Data frame comparing Simon (exact) vs BO (calibrated) designs
