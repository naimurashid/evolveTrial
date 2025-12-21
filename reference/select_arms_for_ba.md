# Select arms for between-arm comparison

Determines which arms should proceed to BA phase.

## Usage

``` r
select_arms_for_ba(
  sa_efficacy_reached,
  active_arms,
  reference_arm,
  selection_strategy = "all_successful"
)
```

## Arguments

- sa_efficacy_reached:

  Named logical vector of SA efficacy status

- active_arms:

  Vector of active arm names

- reference_arm:

  Name of reference arm

- selection_strategy:

  Strategy: "all_successful", "best", or "reference_plus_best"

## Value

Vector of arm names for BA phase
