# Optimal N_add finder

Binary search for minimum N_add achieving target PP.

## Usage

``` r
find_optimal_n_add(
  state,
  target_pp,
  theta,
  base_args,
  scenario_params = NULL,
  max_n = 100,
  step = 10
)
```

## Arguments

- state:

  Current trial state

- target_pp:

  Target predictive probability

- theta:

  Hybrid design parameters

- base_args:

  evolveTrial base configuration

- scenario_params:

  Scenario parameters

- max_n:

  Maximum N to consider

- step:

  Step size for search

## Value

Optimal N_add or NA if not achievable
