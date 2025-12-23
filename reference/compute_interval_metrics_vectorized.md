# Compute interval metrics (VECTORIZED)

Vectorized computation of events and exposure per interval. Much faster
than the row-by-row loop version.

## Usage

``` r
compute_interval_metrics_vectorized(
  registry,
  analysis_time,
  interval_cutpoints
)
```

## Arguments

- registry:

  Patient registry data frame

- analysis_time:

  Current analysis time

- interval_cutpoints:

  PWE interval boundaries

## Value

List with events_per_interval and exposure_per_interval
