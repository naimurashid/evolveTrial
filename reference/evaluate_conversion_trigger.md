# Evaluate conversion trigger

Checks if conditions are met to transition from SA to BA phase.

## Usage

``` r
evaluate_conversion_trigger(
  sa_efficacy_reached,
  active_arms,
  trigger = "any_single_success",
  k_required = 1
)
```

## Arguments

- sa_efficacy_reached:

  Named logical vector of SA efficacy status per arm

- active_arms:

  Vector of active arm names

- trigger:

  Trigger type: "any_single_success", "all_single_success", or "k_of_K"

- k_required:

  Number required for k_of_K trigger

## Value

List with triggered status and reason
