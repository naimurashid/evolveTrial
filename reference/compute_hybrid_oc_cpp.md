# Compute hybrid trial operating characteristics (C++ implementation)

Compute hybrid trial operating characteristics (C++ implementation)

## Usage

``` r
compute_hybrid_oc_cpp(
  n_sim,
  theta_list,
  base_args_list,
  scenario_params_list,
  null_scenario_list
)
```

## Arguments

- n_sim:

  Number of simulations

- theta_list:

  Design parameters

- base_args_list:

  Base configuration

- scenario_params_list:

  Scenario parameters (alternative hypothesis)

- null_scenario_list:

  Scenario parameters (null hypothesis)

## Value

List with power, type1, EN, conversion rate (including per-arm metrics)
