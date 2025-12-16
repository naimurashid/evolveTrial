# Sample medians for vs-ref independent model (DISPATCHER)

Routes to C++ or R implementation based on EVOLVETRIAL_USE_CPP. Samples
hazards independently for each arm (no borrowing).

## Usage

``` r
sample_vs_ref_medians_independent(
  slCtrl,
  slTrt,
  args,
  num_samples,
  ctrl_cache = NULL
)
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
  prior_beta_params_model, interval_cutpoints_sim

- num_samples:

  Number of posterior samples to draw

- ctrl_cache:

  Optional cached control arm posteriors for multi-arm optimization
  (list with lamC, medCtrl, interval_lengths)

## Value

List with medCtrl, medTrt (posterior median samples), and logHR (NULL
for independent model)
