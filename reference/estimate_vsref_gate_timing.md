# Estimate when the vs-reference interim gates can be satisfied

Provides lower-bound heuristics for when the between-arm (vs-reference)
interim gating criteria can all be satisfied under deterministic
accrual.

## Usage

``` r
estimate_vsref_gate_timing(args)
```

## Arguments

- args:

  A named list following the structure passed into
  [`run_simulation_pure()`](run_simulation_pure.md) /
  [`run_scenarios()`](run_scenarios.md). The function uses entries
  relevant to vs-reference gating, namely `arm_names`,
  `reference_arm_name`, `overall_accrual_rate`, `randomization_probs`,
  `max_total_patients_per_arm`, `max_follow_up_sim`,
  `min_events_per_arm`, `min_median_followup_per_arm`, and
  `min_person_time_frac_per_arm`.

## Value

A list with two components:

- `per_arm`: data frame containing the per-arm accrual rate and the
  heuristic lower-bound times (in months) for each gating component.

- `joint_lower_bound`: the maximum of the per-arm lower bounds,
  representing the earliest calendar time at which all gates could
  plausibly be satisfied simultaneously under the heuristics.

## Details

The calculations assume constant accrual with rate
`overall_accrual_rate` split according to `randomization_probs`.
Person-time requirements are approximated using the relationship PT
\\\approx\\ rate \\\times\\ time\\^2\\ / 2 under steady accrual, while
the median follow-up requirement is approximated using the rule-of-thumb
that the median follow-up cannot exceed roughly half of the calendar
time under uniform accrual.

## Examples

``` r
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
#> $per_arm
#>             Arm AccrualRate MinEvents TimeForEvents MinMedianFollowup
#> Doublet Doublet         1.5         8      5.333333                 3
#> Triplet Triplet         1.5         8      5.333333                 3
#>         TimeForMedianFollowup MinPersonTimeMonths TimeForPersonTime
#> Doublet                     6                 252           18.3303
#> Triplet                     6                 252           18.3303
#> 
#> $joint_lower_bound
#> [1] 18.3303
#> 
```
