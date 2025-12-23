# Simulate future arm data under PWE model (VECTORIZED)

Vectorized version that simulates all patients at once. Much faster than
the patient-by-patient loop.

## Usage

``` r
simulate_future_arm_pwe_vectorized(
  n_patients,
  lambda,
  interval_cutpoints,
  accrual_rate,
  followup
)
```

## Arguments

- n_patients:

  Number of patients to simulate

- lambda:

  True hazard rates per interval

- interval_cutpoints:

  Interval boundaries

- accrual_rate:

  Patients per month

- followup:

  Follow-up time in months

## Value

List with events and exposure per interval
