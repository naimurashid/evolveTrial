# Simulate future arm data under PWE model

Simulate future arm data under PWE model

## Usage

``` r
simulate_future_arm_data(
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
