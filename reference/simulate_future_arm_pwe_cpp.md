# Simulate future arm data under PWE model (C++)

Simulate future arm data under PWE model (C++)

## Usage

``` r
simulate_future_arm_pwe_cpp(
  n_patients,
  lambda,
  interval_cutpoints,
  accrual_rate,
  followup
)
```

## Arguments

- n_patients:

  Number of future patients

- lambda:

  Hazard rates per interval

- interval_cutpoints:

  Interval boundaries

- accrual_rate:

  Patients per time unit

- followup:

  Follow-up time after enrollment complete

## Value

List with events and exposure vectors
