# Simulate survival time from PWE model

Wrapper that delegates to simulate_pwe_survival() in hybrid_trial.R to
avoid code duplication.

## Usage

``` r
simulate_pwe_time(lambda, interval_cutpoints)
```

## Arguments

- lambda:

  Hazard rates per interval

- interval_cutpoints:

  Interval boundaries

## Value

Survival time (Inf if no event possible due to zero hazards)
