# Sample medians for vs-ref independent model (C++ implementation)

End-to-end posterior sampler for independent hazards model. Used for
single-arm trials and multi-arm trials with independent hazards.

## Usage

``` r
sample_vs_ref_medians_independent_cpp(
  E_C,
  PT_C,
  E_T,
  PT_T,
  alpha_prior,
  beta_prior,
  interval_lengths,
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

- num_samples:

  Number of posterior samples

## Value

List with "medCtrl" and "medTrt" vectors (no logHR for independent
model)
