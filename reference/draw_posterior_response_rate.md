# Draw posterior samples for binary response rate

Uses Beta-Binomial conjugacy to sample from the posterior distribution
of the response rate.

## Usage

``` r
draw_posterior_response_rate(
  n_responses,
  n_total,
  alpha_prior = 1,
  beta_prior = 1,
  num_samples = 1000
)
```

## Arguments

- n_responses:

  Number of responders

- n_total:

  Total number of patients

- alpha_prior:

  Beta prior shape1 parameter (default 1 for uniform)

- beta_prior:

  Beta prior shape2 parameter (default 1 for uniform)

- num_samples:

  Number of posterior samples to draw

## Value

Numeric vector of posterior samples for response rate

## Details

Prior: Beta(alpha_prior, beta_prior) Likelihood: Binomial(n, p)
Posterior: Beta(alpha_prior + successes, beta_prior + failures)
