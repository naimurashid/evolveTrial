# Simulate binary response data for a cohort

Generates binary response outcomes for n patients with true response
probability p.

## Usage

``` r
simulate_binary_response_data(n, p, start_id = 1L)
```

## Arguments

- n:

  Number of patients

- p:

  True response probability (0 to 1)

- start_id:

  Starting patient ID (default 1)

## Value

Data frame with columns: id, response (0/1)
