# Binary interim decision check

Evaluates whether to stop for efficacy or futility at an interim look in
a binary endpoint trial.

## Usage

``` r
binary_interim_decision(n_responses, n_total, args, diagnostics = FALSE)
```

## Arguments

- n_responses:

  Number of responders

- n_total:

  Total enrolled

- args:

  Trial arguments containing thresholds

- diagnostics:

  Print diagnostic messages

## Value

List with decision ("continue", "stop_efficacy", "stop_futility") and
probabilities
