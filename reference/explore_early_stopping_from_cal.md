# Explore early stopping knobs around a calibrated design

Sweeps futility thresholds, information gates, and interim schedules
around a calibrated design to characterise trade-offs between alpha,
power, and PETs.

## Usage

``` r
explore_early_stopping_from_cal(
  cal,
  base_args,
  null_med,
  alt_med,
  base = c("null+delta", "alt"),
  futility_delta_grid = c(0, 1, 2, 3),
  fut_thr_grid = c(0.6, 0.7, 0.8, 0.9),
  min_events_grid = c(12, 18),
  min_medFU_grid = c(3, 4.5),
  schedule_modes = c("calendar", "persontime"),
  beat_grid = c(3, 6),
  pt_milestones_choices = list(c(0.3, 0.45, 0.6, 0.8, 1)),
  latest_calendar_look_grid = c(Inf),
  min_events_per_arm_grid = c(8, 12),
  min_median_followup_per_arm_grid = c(0, 4.5),
  min_person_time_frac_per_arm_grid = c(0, 0.25),
  sims = 400,
  seed = 123,
  parallel = (.Platform$OS.type == "unix")
)
```

## Arguments

- cal:

  Calibration output from [`grid_calibrate()`](grid_calibrate.md) or
  similar.

- base_args:

  Baseline argument list.

- null_med:

  Control-arm median under null.

- alt_med:

  Experimental median under alternative.

- base:

  Futility baseline mode: "null+delta" or "alt".

- futility_delta_grid:

  Numeric vector of deltas added to null medians.

- fut_thr_grid:

  Numeric vector of futility probability thresholds.

- min_events_grid:

  Global minimum event counts to consider.

- min_medFU_grid:

  Global minimum median follow-up values.

- schedule_modes:

  Character vector: "calendar", "persontime", or both.

- beat_grid:

  Calendar beat schedules to evaluate (months).

- pt_milestones_choices:

  List of person-time milestone vectors.

- latest_calendar_look_grid:

  Backstop calendar times for person-time schedules.

- min_events_per_arm_grid:

  Per-arm minimum events for gating.

- min_median_followup_per_arm_grid:

  Per-arm minimum median follow-up.

- min_person_time_frac_per_arm_grid:

  Per-arm minimum person-time fractions.

- sims:

  Number of simulations per configuration.

- seed:

  RNG seed.

- parallel:

  Logical; use parallel processing.

## Value

data.table of operating characteristics for each configuration.
