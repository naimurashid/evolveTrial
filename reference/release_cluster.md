# Release an evolveTrial simulation cluster

Stops and releases resources for a cluster created by
[`create_simulation_cluster()`](create_simulation_cluster.md).

## Usage

``` r
release_cluster(cluster)
```

## Arguments

- cluster:

  A cluster object created by
  [`create_simulation_cluster()`](create_simulation_cluster.md) or
  [`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html).

## Value

NULL invisibly.

## See also

[`create_simulation_cluster()`](create_simulation_cluster.md)
