# hybrid_decisions.R Decision rules for hybrid single-arm to between-arm trials

Implements the decision logic for:

1.  Single-arm efficacy/futility decisions

2.  Conversion trigger evaluation

3.  Between-arm efficacy/futility decisions

4.  Overall trial conclusions Evaluate single-arm efficacy for an arm

## Usage

``` r
evaluate_sa_efficacy(
  p_single,
  eff_threshold = 0.9,
  min_events = 15,
  current_events = 0
)
```

## Arguments

- p_single:

  Posterior probability P(HR \< c \| data)

- eff_threshold:

  Efficacy threshold (default 0.90)

- min_events:

  Minimum events required (default 15)

- current_events:

  Current number of events

## Value

List with decision and reason
