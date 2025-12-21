# Determine if trial concluded with overall success

Success is defined as:

- SA efficacy reached AND conversion not needed (NO_GO due to already
  effective)

- OR BA efficacy reached

## Usage

``` r
is_trial_success(state, theta)
```

## Arguments

- state:

  Final trial state

- theta:

  Hybrid design parameters

## Value

Logical
