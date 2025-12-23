# Simulate multiple survival times from PWE model (VECTORIZED)

Generates n survival times in a single vectorized operation. Much faster
than calling simulate_pwe_survival() n times.

## Usage

``` r
simulate_pwe_survival_batch(n, lambda, interval_cutpoints)
```

## Arguments

- n:

  Number of survival times to simulate

- lambda:

  Hazard rates per interval

- interval_cutpoints:

  Interval boundaries

## Value

Vector of n survival times
