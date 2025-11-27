# Render the scenario summary table to a PNG image

Render the scenario summary table to a PNG image

## Usage

``` r
export_scenario_table_to_png(
  results_df,
  file_path = "scenario_summary.png",
  title = "Bayesian Adaptive Design Summary",
  subtitle = NULL,
  highlight_arm = "Triplet",
  snapshot_engine = c("auto", "webshot2", "webshot"),
  vwidth = 1400,
  vheight = 600,
  zoom = 1
)
```

## Arguments

- results_df:

  Data frame returned by [`run_scenarios()`](run_scenarios.md).

- file_path:

  Output path for the PNG image.

- title:

  Main title for the table.

- subtitle:

  Optional subtitle.

- highlight_arm:

  Arm name whose columns should be highlighted.

- snapshot_engine:

  Rendering backend; one of `"auto"`, `"webshot2"`, or `"webshot"`.

- vwidth:

  Viewport width passed to the renderer.

- vheight:

  Viewport height passed to the renderer.

- zoom:

  Zoom factor passed to the renderer.

## Value

Invisibly returns `file_path`. Writes a PNG (and temporary HTML if
`snapshot_engine = "webshot"`).
