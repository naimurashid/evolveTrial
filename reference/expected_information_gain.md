# Expected information gain from additional patients

Computes expected reduction in posterior variance from adding n_add
patients.

## Usage

``` r
expected_information_gain(state, n_add, theta, base_args)
```

## Arguments

- state:

  Current trial state

- n_add:

  Additional patients per arm

- theta:

  Hybrid design parameters

- base_args:

  evolveTrial base configuration

## Value

Expected information gain (relative to current)
