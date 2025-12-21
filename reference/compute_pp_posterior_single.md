# Compute posterior-based PP for single N_add

Uses normal approximation to log(HR) to compute conditional power.

## Usage

``` r
compute_pp_posterior_single(hr_hat, current_n, n_add, total_events, eff_ba)
```

## Arguments

- hr_hat:

  Current HR point estimate

- current_n:

  Current total sample size

- n_add:

  Additional patients per arm

- total_events:

  Current total events

- eff_ba:

  BA efficacy threshold

## Value

Conditional power (PP approximation)
