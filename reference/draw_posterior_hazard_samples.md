# Draw posterior hazard samples (DISPATCHER)

Routes to C++ or R implementation based on EVOLVETRIAL_USE_CPP. Default:
C++ for better performance.

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

## Arguments

- num_intervals:

  Number of intervals

- events_per_interval:

  Events per interval

- person_time_per_interval:

  Person-time per interval

- prior_alpha_params:

  Gamma prior shapes

- prior_beta_params:

  Gamma prior rates

- num_samples:

  Number of samples to draw

## Value

Matrix of posterior hazard samples
