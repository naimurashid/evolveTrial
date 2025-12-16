# Sample medians for vs-ref PH model (DISPATCHER)

Routes to C++ or R implementation based on EVOLVETRIAL_USE_CPP. Uses
proportional hazards model for borrowing strength across arms.

## Usage

``` r
sample_vs_ref_medians_ph(slCtrl, slTrt, args, num_samples)
```

## Arguments

- slCtrl:

  Control arm slice (list with metrics\$events_per_interval,
  metrics\$person_time_per_interval)

- slTrt:

  Treatment arm slice (list with metrics\$events_per_interval,
  metrics\$person_time_per_interval)

- args:

  List of simulation args including prior_alpha_params_model,
  prior_beta_params_model, interval_cutpoints_sim, ph_loghr_prior_mean,
  ph_loghr_prior_sd, num_posterior_draws

- num_samples:

  Number of posterior samples to draw

## Value

List with medCtrl, medTrt (posterior median samples), and logHR (log
hazard ratio samples)
