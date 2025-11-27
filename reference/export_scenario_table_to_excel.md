# Export a scenario summary table to Excel

Export a scenario summary table to Excel

## Usage

``` r
export_scenario_table_to_excel(pretty_tbl, file_path = "scenario_summary.xlsx")
```

## Arguments

- pretty_tbl:

  Data frame produced by
  [`pretty_scenario_matrix()`](pretty_scenario_matrix.md).

- file_path:

  Output path for the Excel workbook.

## Value

Invisibly returns `file_path`. Writes an `.xlsx` file to disk.
