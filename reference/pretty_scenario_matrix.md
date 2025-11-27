# Summarise simulation output by scenario and arm

Aggregates the per-arm results returned by
[`run_scenarios()`](run_scenarios.md) and pivots them into a wide table
(one row per scenario) for downstream reporting.

## Usage

``` r
pretty_scenario_matrix(results_df)
```

## Arguments

- results_df:

  Data frame/data.table produced by
  [`run_scenarios()`](run_scenarios.md) containing at least `scenario`
  and `Arm_Name`.

## Value

A data.frame with one row per scenario and columns for each arm's key
operating characteristics.
