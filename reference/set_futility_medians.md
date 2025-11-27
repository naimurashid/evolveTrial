# Adjust futility medians for experimental arms

Adjust futility medians for experimental arms

## Usage

``` r
set_futility_medians(
  args,
  null_med,
  alt_med,
  base = c("null+delta", "alt"),
  delta = 0
)
```

## Arguments

- args:

  Argument list whose `futility_median_arms` entry will be updated.

- null_med:

  Numeric scalar representing the null median.

- alt_med:

  Numeric scalar representing the alternative median.

- base:

  Character scalar selecting `"null+delta"` or `"alt"`.

- delta:

  Numeric offset added when `base = "null+delta"`.

## Value

Modified argument list with updated `futility_median_arms`.
