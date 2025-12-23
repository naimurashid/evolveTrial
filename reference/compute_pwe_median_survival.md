# Compute median survival from PWE hazards

Wrapper that delegates to compute_pwe_median() in hybrid_trial.R to
avoid code duplication.

## Usage

``` r
compute_pwe_median_survival(lambda, interval_cutpoints)
```

## Arguments

- lambda:

  Vector of hazard rates per interval

- interval_cutpoints:

  Interval boundaries

## Value

Median survival time
