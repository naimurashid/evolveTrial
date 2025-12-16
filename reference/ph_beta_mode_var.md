# Compute mode and variance of log-HR posterior for PH model (DISPATCHER)

Compute mode and variance of log-HR posterior for PH model (DISPATCHER)

## Usage

``` r
ph_beta_mode_var(
  E_C,
  PT_C,
  E_T,
  PT_T,
  alpha_prior,
  beta_prior,
  mu,
  sigma,
  tol = 1e-06,
  max_iter = 50
)
```

## Arguments

- E_C:

  Events in control arm (vector by interval)

- PT_C:

  Person-time in control arm (vector by interval)

- E_T:

  Events in treatment arm (vector by interval)

- PT_T:

  Person-time in treatment arm (vector by interval)

- alpha_prior:

  Prior shape parameter for hazard (vector by interval)

- beta_prior:

  Prior rate parameter for hazard (vector by interval)

- mu:

  Prior mean for log-HR

- sigma:

  Prior SD for log-HR

- tol:

  Tolerance for Newton-Raphson convergence (default 1e-6)

- max_iter:

  Maximum iterations for Newton-Raphson (default 50)

## Value

Named list with mode and variance of log-HR posterior
