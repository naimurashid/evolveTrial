# Calculate median survival for multiple hazard samples (C++ implementation)

Vectorized version that processes an entire matrix of posterior hazard
samples. Each row is a posterior sample; each column is an interval.

## Usage

``` r
calculate_median_survival_matrix_cpp(hazard_samples, interval_lengths)
```

## Arguments

- hazard_samples:

  Matrix with num_samples rows and num_intervals columns

- interval_lengths:

  Numeric vector of interval lengths (durations)

## Value

Numeric vector of median survival times (length num_samples)
