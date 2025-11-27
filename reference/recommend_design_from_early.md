# Recommend a single early-stopping design

Selects the top-performing row from an early-stopping grid subject to
alpha/power (and optionally PET) constraints.

## Usage

``` r
recommend_design_from_early(
  df,
  alpha_cap = 0.1,
  power_floor = 0.8,
  pet_fut_cap = NULL
)
```

## Arguments

- df:

  data.table/data.frame from
  [`explore_early_stopping_from_cal()`](explore_early_stopping_from_cal.md).

- alpha_cap:

  Maximum acceptable type I error.

- power_floor:

  Minimum acceptable power.

- pet_fut_cap:

  Optional cap on alternative PET for futility.

## Value

data.table row describing the recommended design.
