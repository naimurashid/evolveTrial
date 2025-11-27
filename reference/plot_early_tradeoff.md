# Plot early-stopping trade-offs

Visualises power versus type I error for the early-stopping exploration
grid.

## Usage

``` r
plot_early_tradeoff(
  early_df,
  target_alpha = 0.1,
  fix_min_ev = NULL,
  fix_mfu = NULL,
  fix_beat = NULL
)
```

## Arguments

- early_df:

  data.table/data.frame returned by
  [`explore_early_stopping_from_cal()`](explore_early_stopping_from_cal.md).

- target_alpha:

  Type I error target to show as a reference line.

- fix_min_ev:

  Optional scalar to filter `min_events`.

- fix_mfu:

  Optional scalar to filter `min_medFU`.

- fix_beat:

  Optional scalar to filter calendar beats.

## Value

A ggplot object (or the filtered data when `ggplot2` is unavailable).
