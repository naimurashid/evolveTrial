# Evaluate PH-based grid of designs

Helper that sweeps probability thresholds, gate settings, and HR
margins, returning operating characteristics plus expected sample sizes
per arm.

## Usage

``` r
evaluate_ph_grid(
  base_args,
  grid,
  scens,
  sims = 2000,
  seed = 4242,
  parallel = TRUE
)
```

## Arguments

- base_args:

  Baseline argument list passed to
  [`run_scenarios()`](run_scenarios.md).

- grid:

  data.table/data.frame describing the design grid. Must include columns
  `label`, `thr_eff`, `thr_fut`, `margin`, `min_ev`, `min_pt`, and
  optionally `hr_margin`.

- scens:

  Scenario list (e.g., from
  [`scenarios_from_grid()`](scenarios_from_grid.md)).

- sims:

  Number of simulations per grid row.

- seed:

  RNG seed.

- parallel:

  Logical; use parallel execution.

## Value

data.table summarising Type I error, power, PETs, expected N per arm,
control-arm expectations, and total expected N under null/alt. Evaluate
proportional-hazards vs-reference designs over a grid

Runs the supplied grid of thresholds/margins against a scenario list,
collecting operating characteristics for each design.

data.table summarising alpha, power, PETs, and expected sample sizes for
each labelled design.
