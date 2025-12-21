# Compute predictive probability with full PWE model

Uses Monte Carlo sampling of median survival times for more accurate PWE
comparison.

## Usage

``` r
compute_pp_predictive_pwe(
  state,
  n_add,
  theta,
  base_args,
  scenario_params,
  n_outer = 1000,
  n_inner = 250
)
```

## Arguments

- state:

  Current hybrid trial state

- n_add:

  Additional patients per arm

- theta:

  Hybrid design parameters

- base_args:

  evolveTrial base configuration

- scenario_params:

  Scenario parameters

- n_outer:

  Number of outer samples

- n_inner:

  Number of inner samples for posterior

## Value

Predictive probability
