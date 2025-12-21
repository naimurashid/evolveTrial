# Compute PP curve using predictive probability method

Full Monte Carlo integration over posterior uncertainty. For each
candidate N_add:

1.  Draw "true" parameters from current posterior

2.  Simulate future patients under those parameters

3.  Update posterior with simulated data

4.  Compute P(BA success \| final data)

5.  Average over outer samples

## Usage

``` r
compute_pp_curve_predictive(
  state,
  n_candidates,
  theta,
  base_args,
  scenario_params = NULL
)
```

## Arguments

- state:

  Current hybrid trial state

- n_candidates:

  Vector of candidate N_add values

- theta:

  Hybrid design parameters

- base_args:

  evolveTrial base configuration

- scenario_params:

  Scenario parameters

## Value

Data frame with n_add and pp columns
