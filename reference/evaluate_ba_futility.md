# Evaluate between-arm futility

Evaluate between-arm futility

## Usage

``` r
evaluate_ba_futility(
  p_between,
  fut_threshold = 0.05,
  min_events = 30,
  current_events = 0
)
```

## Arguments

- p_between:

  Posterior probability P(HR_AB \< 1 \| data)

- fut_threshold:

  Futility threshold (default 0.05)

- min_events:

  Minimum total events required (default 30)

- current_events:

  Current total events across arms

## Value

List with decision and reason
