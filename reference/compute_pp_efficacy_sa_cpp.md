# Compute predictive probability for single-arm efficacy (C++)

Monte Carlo computation of PP that single arm will reach efficacy
threshold vs historical control at the final analysis.

## Usage

``` r
compute_pp_efficacy_sa_cpp(
  a_arm,
  b_arm,
  hist_hazard,
  hr_threshold,
  n_add,
  interval_cutpoints,
  accrual_rate,
  followup,
  eff_threshold,
  n_outer = 200L,
  use_antithetic = TRUE
)
```

## Arguments

- a_arm:

  Current posterior shape for arm (per interval)

- b_arm:

  Current posterior rate for arm (per interval)

- hist_hazard:

  Historical hazard rates (per interval)

- hr_threshold:

  Hazard ratio threshold for efficacy

- n_add:

  Number of additional patients to simulate

- interval_cutpoints:

  PWE interval boundaries

- accrual_rate:

  Patients per time unit

- followup:

  Follow-up time after enrollment

- eff_threshold:

  Efficacy posterior threshold (e.g., 0.90)

- n_outer:

  Number of outer MC samples

- use_antithetic:

  Use antithetic variates for variance reduction

## Value

Predictive probability of reaching efficacy
