# Filter early-stopping designs by operating targets

Filter early-stopping designs by operating targets

## Usage

``` r
filter_early_grid(early_df, alpha_cap = 0.1, power_floor = 0.7)
```

## Arguments

- early_df:

  data.table/data.frame produced by
  [`explore_early_stopping_from_cal()`](explore_early_stopping_from_cal.md).

- alpha_cap:

  Maximum acceptable type I error.

- power_floor:

  Minimum acceptable power.

## Value

data.table containing the subset that meets the supplied criteria.
