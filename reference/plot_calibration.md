# Plot power versus type I error for calibration grids

Visualises the output of [`grid_calibrate()`](grid_calibrate.md) by
plotting power against type I error, optionally highlighting Pareto
frontiers and feasible designs.

## Usage

``` r
plot_calibration(cal, target_alpha = 0.1, label_top_n = 3)
```

## Arguments

- cal:

  List returned by [`grid_calibrate()`](grid_calibrate.md).

- target_alpha:

  Type I error threshold to display.

- label_top_n:

  Number of feasible designs to annotate.

## Value

A ggplot object.
