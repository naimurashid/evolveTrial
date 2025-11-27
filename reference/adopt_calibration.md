# Adopt a calibrated design configuration

Picks a selected row from [`grid_calibrate()`](grid_calibrate.md) and
returns updated arguments plus a two-scenario list (null vs alternative)
for subsequent exploration.

## Usage

``` r
adopt_calibration(cal, base_args, null_med, alt_med, which = 1L)
```

## Arguments

- cal:

  Output from [`grid_calibrate()`](grid_calibrate.md).

- base_args:

  Baseline argument list.

- null_med:

  Control-arm median under the null.

- alt_med:

  Experimental median under the alternative.

- which:

  Integer index specifying which row of `cal$top` to adopt.

## Value

A list containing the updated arguments (`args_star`), the selected row
(`pick`), and a two-scenario list (`scens2`).
