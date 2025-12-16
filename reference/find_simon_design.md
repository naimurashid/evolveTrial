# Find optimal Simon design

Searches for the Simon optimal or minimax design meeting constraints.
Uses optimized search with early termination.

## Usage

``` r
find_simon_design(
  p0,
  p1,
  alpha = 0.1,
  beta = 0.2,
  n_max = 100,
  criterion = c("optimal", "minimax")
)
```

## Arguments

- p0:

  Null response rate

- p1:

  Alternative response rate

- alpha:

  Maximum type I error

- beta:

  Maximum type II error (1 - power)

- n_max:

  Maximum total sample size to search

- criterion:

  "optimal" (minimize E\[N\] under null) or "minimax" (minimize max N)

## Value

Data frame with design parameters and operating characteristics
