# Compute PP curve using posterior method

Fast analytical approximation using current posterior point estimate.
Less accurate but much faster than predictive method.

## Usage

``` r
compute_pp_curve_posterior(state, n_candidates, theta, base_args)
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

## Value

Data frame with n_add and pp columns
