# Calculate interval-specific metrics from patient data

Recalculates events and person-time using a more efficient data.table
approach.

## Usage

``` r
calculate_interval_metrics_fast(patient_data, interval_cutpoints)
```

## Arguments

- patient_data:

  Data frame with columns `observed_time` and `event_status`.

- interval_cutpoints:

  Numeric vector of interval boundaries.

## Value

A list with `events_per_interval` and `person_time_per_interval`
vectors.
