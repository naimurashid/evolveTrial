# Make conversion decision based on PP curve

Implements the two-threshold decision rule:

- GO: PP \>= pp_go for some viable N

- NO_GO: max(PP) \< pp_nogo

- AMBIGUOUS: pp_nogo \<= max(PP) \< pp_go (default to NO_GO)

## Usage

``` r
make_conversion_decision(pp_curve, theta)
```

## Arguments

- pp_curve:

  Data frame with n_add and pp columns

- theta:

  Hybrid design parameters

## Value

List with decision, n_add, pp, and reason
