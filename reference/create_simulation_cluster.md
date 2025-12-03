# Create a reusable evolveTrial simulation cluster

Creates a parallel cluster optimized for evolveTrial simulations. The
cluster can be reused across multiple calls to
[`run_simulation_pure()`](run_simulation_pure.md) by passing it via the
`cluster` parameter, avoiding the overhead of repeated cluster
creation/destruction.

## Usage

``` r
create_simulation_cluster(
  workers = NULL,
  cluster_type = c("auto", "FORK", "PSOCK")
)
```

## Arguments

- workers:

  Integer; number of worker processes. Defaults to
  `parallel::detectCores() - 1`.

- cluster_type:

  One of `"auto"` (default), `"FORK"`, or `"PSOCK"`. Auto selects FORK
  on Unix systems (faster) and PSOCK on Windows.

## Value

A parallel cluster object that can be passed to
[`run_simulation_pure()`](run_simulation_pure.md) via the `cluster`
parameter.

## Details

This function is particularly useful for Bayesian optimization workflows
where [`run_simulation_pure()`](run_simulation_pure.md) is called
hundreds of times. By creating the cluster once and reusing it, you can
eliminate the 5-10 second cluster spawn/teardown overhead per call.

FORK clusters (default on Unix) are significantly faster because:

- Worker processes inherit the parent environment (no package loading)

- Data is shared via copy-on-write (no serialization overhead)

Remember to call [`release_cluster()`](release_cluster.md) when done to
free resources.

## See also

[`release_cluster()`](release_cluster.md),
[`run_simulation_pure()`](run_simulation_pure.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Create cluster once
cl <- create_simulation_cluster(workers = 8)

# Use in repeated BO evaluations
for (i in 1:100) {
  result <- run_simulation_pure(
    num_simulations = 500,
    ...,
    parallel_replicates = TRUE,
    cluster = cl
  )
}

# Clean up
release_cluster(cl)
} # }
```
