# Simulate survival time from PWE model

Uses inverse CDF method with piecewise constant hazard.

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

Survival time
