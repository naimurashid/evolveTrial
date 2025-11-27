# Calculate median survival for piecewise exponential model (C++ implementation)

Computes the median survival time for a piecewise exponential model. If
the 0.5 survival is not reached by the end of the last interval,
continues past the last cutpoint with the last interval's hazard
(open-ended tail). Only returns Inf if the last hazard is exactly zero.

## Usage

``` r
calculate_median_survival_piecewise_cpp(hazard_rates, interval_lengths)
```

## Arguments

- hazard_rates:

  Numeric vector of hazard rates for each interval

- interval_lengths:

  Numeric vector of interval lengths (durations)

## Value

Median survival time (numeric scalar, possibly Inf)
