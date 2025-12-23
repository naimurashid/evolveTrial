# Simulate PWE survival times (C++ vectorized)

Simulate PWE survival times (C++ vectorized)

## Usage

``` r
simulate_pwe_survival_batch_cpp(n, lambda, interval_cutpoints)
```

## Arguments

- n:

  Number of survival times to simulate

- lambda:

  Hazard rates per interval

- interval_cutpoints:

  Interval boundaries (length = n_intervals + 1)

## Value

Vector of n survival times
