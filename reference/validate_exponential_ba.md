# Validate MC against closed-form for exponential

Validate MC against closed-form for exponential

## Usage

``` r
validate_exponential_ba(a_exp, b_exp, a_ref, b_ref, n_samples = 10000)
```

## Arguments

- a_exp:

  Experimental arm shape

- b_exp:

  Experimental arm rate

- a_ref:

  Reference arm shape

- b_ref:

  Reference arm rate

- n_samples:

  MC sample size

## Value

List with closed_form and mc_estimate
