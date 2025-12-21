# Update hybrid trial state based on current conditions

This is the main state machine driver. It evaluates the current state
and applies appropriate transitions.

## Usage

``` r
update_hybrid_state(state, theta, base_args, scenario_params)
```

## Arguments

- state:

  Current trial state

- theta:

  Hybrid design parameters

- base_args:

  evolveTrial base configuration

- scenario_params:

  Scenario parameters (true hazards, etc.)

## Value

Updated trial state
