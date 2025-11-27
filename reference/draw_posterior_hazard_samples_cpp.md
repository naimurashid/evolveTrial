# Draw posterior hazard samples (C++ implementation)

Draws samples from the posterior distribution of hazard rates for each
interval in a Bayesian piecewise exponential model with Gamma priors.

## Usage

``` r
draw_posterior_hazard_samples_cpp(
  num_intervals,
  events_per_interval,
  person_time_per_interval,
  prior_alpha_params,
  prior_beta_params,
  num_samples
)
```

## Arguments

- num_intervals:

  Number of intervals in piecewise model

- events_per_interval:

  Integer vector of observed events per interval

- person_time_per_interval:

  Numeric vector of person-time at risk per interval

- prior_alpha_params:

  Numeric vector of Gamma prior shape parameters

- prior_beta_params:

  Numeric vector of Gamma prior rate parameters

- num_samples:

  Number of posterior samples to draw

## Value

Matrix with num_samples rows and num_intervals columns of hazard samples
