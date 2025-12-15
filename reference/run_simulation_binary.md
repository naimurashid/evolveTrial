# Run binary endpoint trial simulation

Simulates a two-stage binary endpoint trial with Bayesian decision
rules. Supports Simon-style designs with cohort-based enrollment.

## Usage

``` r
run_simulation_binary(
  num_simulations,
  arm_names,
  true_response_prob,
  n1_per_arm,
  n_total_per_arm,
  p0,
  p1 = NULL,
  efficacy_threshold_binary_prob = 0.95,
  futility_threshold_binary_prob = 0.95,
  r1_per_arm = NULL,
  r_per_arm = NULL,
  use_simon_rules = FALSE,
  alpha_prior = 1,
  beta_prior = 1,
  disable_interim_eff_stop = FALSE,
  diagnostics = FALSE,
  progress = interactive()
)
```

## Arguments

- num_simulations:

  Number of Monte Carlo replicates

- arm_names:

  Character vector of arm names (typically single arm)

- true_response_prob:

  Named numeric vector of true response probabilities

- n1_per_arm:

  Named integer vector of stage 1 sample sizes

- n_total_per_arm:

  Named integer vector of total sample sizes

- p0:

  Null hypothesis response rate

- p1:

  Alternative hypothesis response rate

- efficacy_threshold_binary_prob:

  Posterior probability threshold for efficacy

- futility_threshold_binary_prob:

  Posterior probability threshold for futility

- r1_per_arm:

  Optional: Simon-style stage 1 boundaries (stop if X1 \<= r1)

- r_per_arm:

  Optional: Simon-style total boundaries (success if X \> r)

- use_simon_rules:

  Use exact Simon counting rules instead of Bayesian

- alpha_prior:

  Beta prior shape1 (default 1)

- beta_prior:

  Beta prior shape2 (default 1)

- disable_interim_eff_stop:

  If TRUE, do not stop for efficacy at interim (only futility stopping
  at interim, efficacy only at final). This makes the Bayesian design
  comparable to Simon's two-stage design which only stops for futility
  at the interim.

- diagnostics:

  Print diagnostic messages

## Value

Data frame with operating characteristics per arm
