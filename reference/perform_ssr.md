# hybrid_ssr.R Sample Size Re-Estimation methods for hybrid trials

Two methods are provided:

1.  Predictive Probability: Full Monte Carlo integration over posterior
    uncertainty

2.  Posterior Probability: Fast analytical approximation using point
    estimates Perform sample size re-estimation for hybrid trial

Main entry point for SSR calculation. Dispatches to appropriate method
based on theta\$ss_method.

## Usage

``` r
perform_ssr(state, theta, base_args, scenario_params = NULL)
```

## Arguments

- state:

  Current hybrid trial state

- theta:

  Hybrid design parameters

- base_args:

  evolveTrial base configuration

- scenario_params:

  Scenario parameters (for predictive method)

## Value

List with:

- decision: "GO", "NO_GO", or "AMBIGUOUS"

- n_add_recommended: Selected additional N per arm

- pp_at_decision: PP value at decision point

- pp_curve: Full PP curve data frame

- method: Method used ("predictive" or "posterior")

## Details

Implements both predictive probability (recommended) and posterior
probability (fast approximation) methods for determining additional
sample size needs at the conversion decision point.
