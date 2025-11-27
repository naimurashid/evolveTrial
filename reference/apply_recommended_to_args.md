# Apply a recommended early-stopping configuration to the argument list

Apply a recommended early-stopping configuration to the argument list

## Usage

``` r
apply_recommended_to_args(args_star, rec_row)
```

## Arguments

- args_star:

  Baseline argument list (typically from
  [`adopt_calibration()`](adopt_calibration.md)).

- rec_row:

  Single-row data.table produced by
  [`recommend_design_from_early()`](recommend_design_from_early.md).

## Value

Modified argument list with the recommended early-stopping settings.
