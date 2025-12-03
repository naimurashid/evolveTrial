# Calculate Simon design operating characteristics analytically

Computes exact operating characteristics for a Simon two-stage design.

## Usage

``` r
simon_oc_exact(n1, r1, n, r, p)
```

## Arguments

- n1:

  Stage 1 sample size

- r1:

  Stage 1 futility boundary (stop if X1 \<= r1)

- n:

  Total sample size (n1 + n2)

- r:

  Total response threshold for efficacy (success if X \> r)

- p:

  True response probability

## Value

List with:

- reject_prob: Probability of rejecting null (power or type I error)

- pet: Probability of early termination at stage 1

- en: Expected sample size
