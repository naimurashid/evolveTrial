# Build scenario overrides from a grid of design choices

Expands a named list of options into a list of scenario override lists
that can be merged with `base_args` prior to simulation.

## Usage

``` r
scenarios_from_grid(choices)
```

## Arguments

- choices:

  Named list where each element is either a vector of scalar values or a
  list whose elements are per-arm vectors.

## Value

A list of scenario override lists. The underlying Cartesian grid is
attached as an attribute named `"grid"`.

## Examples

``` r
scenarios_from_grid(list(
  max_total_patients_per_arm = list(
    c(Doublet = 60, Triplet = 60),
    c(Doublet = 60, Triplet = 70)
  ),
  compare_arms_futility_margin = c(0.3, 0.4)
))
#> [[1]]
#> [[1]]$max_total_patients_per_arm
#> Doublet Triplet 
#>      60      60 
#> 
#> [[1]]$compare_arms_futility_margin
#> [1] 0.3
#> 
#> 
#> [[2]]
#> [[2]]$max_total_patients_per_arm
#> Doublet Triplet 
#>      60      70 
#> 
#> [[2]]$compare_arms_futility_margin
#> [1] 0.3
#> 
#> 
#> [[3]]
#> [[3]]$max_total_patients_per_arm
#> Doublet Triplet 
#>      60      60 
#> 
#> [[3]]$compare_arms_futility_margin
#> [1] 0.4
#> 
#> 
#> [[4]]
#> [[4]]$max_total_patients_per_arm
#> Doublet Triplet 
#>      60      70 
#> 
#> [[4]]$compare_arms_futility_margin
#> [1] 0.4
#> 
#> 
#> attr(,"grid")
#>   max_total_patients_per_arm compare_arms_futility_margin
#> 1                          1                          0.3
#> 2                          2                          0.3
#> 3                          1                          0.4
#> 4                          2                          0.4
```
