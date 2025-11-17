# Computes the median survival time for a piecewise exponential model. If the 0.5 survival is not reached by the end of the last interval, we *continue past the last cutpoint* with the last interval's hazard (open-ended tail). Only return Inf if the last hazard is exactly zero.

Computes the median survival time for a piecewise exponential model. If
the 0.5 survival is not reached by the end of the last interval, we
*continue past the last cutpoint* with the last interval's hazard
(open-ended tail). Only return Inf if the last hazard is exactly zero.

## Usage

``` r
calculate_median_survival_piecewise(hazard_rates, interval_lengths)
```
