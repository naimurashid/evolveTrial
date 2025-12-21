# Compute interval-specific metrics from registry

Compute interval-specific metrics from registry

## Usage

``` r
compute_interval_metrics(registry, analysis_time, interval_cutpoints)
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
