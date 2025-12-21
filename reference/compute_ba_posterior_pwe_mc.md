# Monte Carlo comparison for PWE posteriors

Samples from posteriors and compares median survival times.

## Usage

``` r
compute_ba_posterior_pwe_mc(
  a_exp,
  b_exp,
  a_ref,
  b_ref,
  interval_cutpoints,
  n_samples = 1000
)
```

## Arguments

- a_exp:

  Posterior shapes for experimental arm

- b_exp:

  Posterior rates for experimental arm

- a_ref:

  Posterior shapes for reference arm

- b_ref:

  Posterior rates for reference arm

- interval_cutpoints:

  Interval boundaries

- n_samples:

  Number of Monte Carlo samples

## Value

P(median_exp \> median_ref)
