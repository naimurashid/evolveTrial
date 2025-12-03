# Calculate binary endpoint interim probabilities

Computes posterior probabilities for efficacy and futility decisions in
binary endpoint trials.

## Usage

``` r
calculate_binary_probs(
  n_responses,
  n_total,
  p0,
  p1,
  alpha_prior = 1,
  beta_prior = 1
)
```

## Arguments

- n_responses:

  Number of responders observed

- n_total:

  Total patients enrolled

- p0:

  Null hypothesis response rate (for futility)

- p1:

  Alternative hypothesis response rate (for efficacy target)

- alpha_prior:

  Beta prior shape1 (default 1)

- beta_prior:

  Beta prior shape2 (default 1)

## Value

List with pr_eff (P(p \> p0)) and pr_fut (P(p \< p1))
