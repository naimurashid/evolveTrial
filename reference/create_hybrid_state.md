# Create initial hybrid trial state

Create initial hybrid trial state

## Usage

``` r
create_hybrid_state(arm_names, reference_arm, theta, base_args)
```

## Arguments

- arm_names:

  Character vector of arm names (e.g., c("Control", "Experimental"))

- reference_arm:

  Name of reference arm for between-arm comparison

- theta:

  Hybrid design parameters (see create_hybrid_theta)

- base_args:

  evolveTrial base configuration

## Value

List containing trial state
