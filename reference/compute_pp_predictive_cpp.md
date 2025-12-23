# Compute predictive probability for hybrid trial conversion (C++)

Full Monte Carlo computation with optional antithetic variates and early
stopping.

## Usage

``` r
compute_pp_predictive_cpp(
  a_exp,
  b_exp,
  a_ref,
  b_ref,
  n_add,
  interval_cutpoints,
  accrual_rate,
  followup,
  eff_ba,
  pp_go,
  pp_nogo,
  n_outer = 500L,
  use_antithetic = TRUE,
  use_early_stop = TRUE
)
```

## Arguments

- a_exp:

  Current posterior shape for experimental arm (per interval)

- b_exp:

  Current posterior rate for experimental arm (per interval)

- a_ref:

  Current posterior shape for reference arm (per interval)

- b_ref:

  Current posterior rate for reference arm (per interval)

- n_add:

  Number of additional patients to simulate per arm

- interval_cutpoints:

  PWE interval boundaries

- accrual_rate:

  Patients per time unit

- followup:

  Follow-up time after enrollment

- eff_ba:

  Efficacy threshold for BA decision

- pp_go:

  PP threshold for GO decision

- pp_nogo:

  PP threshold for NO-GO decision

- n_outer:

  Number of outer MC samples

- use_antithetic:

  Use antithetic variates for variance reduction

- use_early_stop:

  Enable early stopping based on CI

## Value

Predictive probability of success
