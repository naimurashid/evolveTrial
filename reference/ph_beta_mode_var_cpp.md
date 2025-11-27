# Compute mode and variance of log-HR posterior for PH model (C++ implementation)

Uses Newton-Raphson iteration to find the mode of the log-HR posterior
under a proportional hazards model with Gamma prior on baseline hazards
and normal prior on log-HR.

## Usage

``` r
ph_beta_mode_var_cpp(
  E_C,
  PT_C,
  E_T,
  PT_T,
  alpha_prior,
  beta_prior,
  mu,
  sigma,
  tol = 1e-06,
  max_iter = 50L
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

- mu:

  Log-HR prior mean

- sigma:

  Log-HR prior SD

- tol:

  Convergence tolerance (default 1e-6)

- max_iter:

  Maximum Newton-Raphson iterations (default 50)

## Value

List with "mean" and "sd" of log-HR posterior
