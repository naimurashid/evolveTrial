# Compute predictive probability for a single N_add value

This is the core Monte Carlo algorithm for predictive probability.
Implements early stopping for efficiency when outcome is clear.

## Usage

``` r
compute_pp_predictive_full(
  state,
  n_add,
  theta,
  base_args,
  scenario_params,
  n_outer = 1000
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

  Number of outer Monte Carlo samples

## Value

Predictive probability
