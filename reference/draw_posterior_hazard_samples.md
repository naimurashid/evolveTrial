# Draws samples from the posterior distribution of hazard rates for each interval in a Bayesian piecewise exponential model with Gamma priors.

Draws samples from the posterior distribution of hazard rates for each
interval in a Bayesian piecewise exponential model with Gamma priors.

## Usage

``` r
draw_posterior_hazard_samples(
  num_intervals,
  events_per_interval,
  person_time_per_interval,
  prior_alpha_params,
  prior_beta_params,
  num_samples = 1000
)
```
