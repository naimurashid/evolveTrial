# evolveTrial

The goal of evolveTrial is to …

## Installation

You can install the development version of evolveTrial from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("naimurashid/evolveTrial")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(evolveTrial)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

![](reference/figures/README-pressure-1.png)

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.

## Diagnosing interim gating choices

When early stopping does not seem to occur during simulations, it is
often because the information gates (minimum events, median follow-up,
or required person-time) postpone the first informative interim look
until late in the trial. The helper
[`estimate_vsref_gate_timing()`](reference/estimate_vsref_gate_timing.md)
gives quick lower-bound heuristics for when those gates can all be
satisfied under deterministic accrual.

``` r
library(evolveTrial)

args <- list(
  arm_names = c("Doublet", "Triplet"),
  reference_arm_name = "Doublet",
  overall_accrual_rate = 3,
  randomization_probs = c(Doublet = 0.5, Triplet = 0.5),
  max_total_patients_per_arm = c(Doublet = 70, Triplet = 70),
  max_follow_up_sim = 24,
  min_events_per_arm = 8,
  min_median_followup_per_arm = 3,
  min_person_time_frac_per_arm = 0.15
)

estimate_vsref_gate_timing(args)
```

For the configuration above the joint lower bound is roughly 18 months,
meaning that with 3-month calendar beats the first eligible interim look
will occur near the sixth look. Adjusting the information gates (for
example by lowering the person-time fraction or the minimum median
follow-up) yields earlier opportunities for efficacy/futility stopping.
