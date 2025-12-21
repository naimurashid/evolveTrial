# Compute exact PP for exponential model (for validation)

Uses closed-form solutions based on F-distribution. This is used to
validate Monte Carlo methods.

## Usage

``` r
compute_ba_posterior_exponential(a_exp, b_exp, a_ref, b_ref)
```

## Arguments

- a_exp:

  Posterior shape for experimental arm

- b_exp:

  Posterior rate for experimental arm

- a_ref:

  Posterior shape for reference arm

- b_ref:

  Posterior rate for reference arm

## Value

P(lambda_exp \< lambda_ref \| data)
