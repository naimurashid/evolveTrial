# Calculate posterior probability that response rate is below threshold

P(p \< threshold \| data) using Beta posterior

## Usage

``` r
prob_response_below(
  n_responses,
  n_total,
  threshold,
  alpha_prior = 1,
  beta_prior = 1
)
```

## Arguments

- n_responses:

  Number of responders

- n_total:

  Total number of patients

- threshold:

  Response rate threshold

- alpha_prior:

  Beta prior shape1 parameter (default 1)

- beta_prior:

  Beta prior shape2 parameter (default 1)

## Value

Posterior probability P(p \< threshold)
