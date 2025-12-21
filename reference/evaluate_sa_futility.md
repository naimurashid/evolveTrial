# Evaluate single-arm futility for an arm

Evaluate single-arm futility for an arm

## Usage

``` r
evaluate_sa_futility(
  p_single,
  fut_threshold = 0.1,
  min_events = 15,
  current_events = 0
)
```

## Arguments

- p_single:

  Posterior probability P(HR \< c \| data)

- fut_threshold:

  Futility threshold (default 0.10)

- min_events:

  Minimum events required (default 15)

- current_events:

  Current number of events

## Value

List with decision and reason
