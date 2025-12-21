# hybrid_trial.R Core hybrid single-arm to between-arm Bayesian adaptive trial simulator

This module provides the main simulation engine for hybrid trials that:

1.  Begin with single-arm monitoring vs historical control

2.  Evaluate conversion to between-arm comparison based on predictive
    probability

3.  Continue with between-arm monitoring using seamless data Trial
    states for hybrid design

## Usage

``` r
HYBRID_STATES
```

## Format

An object of class `list` of length 4.

## Details

Implements a 4-state machine: STATE_SINGLE -\> STATE_CONSIDER_CONVERSION
-\> STATE_BETWEEN -\> STATE_STOP
