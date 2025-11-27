# Sample medians for vs-ref PH model (C++ implementation)

End-to-end posterior sampler for proportional hazards model. Samples
log-HR from posterior, then baseline hazards, then computes medians.

## Usage

``` r
sample_vs_ref_medians_ph_cpp(
  E_C,
  PT_C,
  E_T,
  PT_T,
  alpha_prior,
  beta_prior,
  interval_lengths,
  mu,
  sigma,
  num_samples
)
```

## Arguments

- E_C:

  Integer vector of control events per interval

- PT_C:

  Numeric vector of control person-time per interval

- E_T:

  Integer vector of treatment events per interval

- PT_T:

  Numeric vector of treatment person-time per interval

- alpha_prior:

  Numeric vector of Gamma prior shape parameters

- beta_prior:

  Numeric vector of Gamma prior rate parameters

- interval_lengths:

  Numeric vector of interval durations

- mu:

  Log-HR prior mean

- sigma:

  Log-HR prior SD

- num_samples:

  Number of posterior samples

## Value

List with "medCtrl", "medTrt", and "logHR" vectors
