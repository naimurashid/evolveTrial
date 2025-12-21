# Evaluate between-arm efficacy

Evaluate between-arm efficacy

## Usage

``` r
evaluate_ba_efficacy(
  p_between,
  eff_threshold = 0.975,
  min_events = 30,
  current_events = 0
)
```

## Arguments

- p_between:

  Posterior probability P(HR_AB \< 1 \| data)

- eff_threshold:

  Efficacy threshold (default 0.975)

- min_events:

  Minimum total events required (default 30)

- current_events:

  Current total events across arms

## Value

List with decision and reason
