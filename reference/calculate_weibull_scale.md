# Calculate Weibull scale from median survival time

Converts a desired median survival time to the Weibull scale parameter.
R's rweibull(n, shape = k, scale = λ) uses S(t) = exp(-(t/λ)^k). The
population median t~ satisfies 0.5 = exp(-(t~/λ)^k) =\> t~ = λ \* (ln
2)^(1/k). Solving for λ: λ = t~ / (ln 2)^(1/k).

## Usage

``` r
calculate_weibull_scale(desired_median, weibull_shape)
```

## Arguments

- desired_median:

  Target median survival time

- weibull_shape:

  Weibull shape parameter (k)

## Value

Weibull scale parameter (λ)

## Examples

``` r
# For median = 9 months and shape = 1.3:
scale <- calculate_weibull_scale(9, 1.3)  # ~11.93
# Verify: median of rweibull(10000, 1.3, scale) should be ~9
```
