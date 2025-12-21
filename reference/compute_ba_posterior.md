# Compute between-arm posterior probability

Uses closed-form F-distribution for exponential model, Monte Carlo
sampling for PWE.

## Usage

``` r
compute_ba_posterior(a_exp, b_exp, a_ref, b_ref, interval_cutpoints, base_args)
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

- interval_cutpoints:

  Interval boundaries

- base_args:

  evolveTrial base configuration

## Value

P(HR \< 1 \| data)
