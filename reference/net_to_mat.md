# Create a contingency table from a data.frame

This function generates a contingency table from a two- or three-column
`data.frame`, where each row represents the interaction between two
nodes (e.g., site and species) and an optional third column indicates
the weight of the interaction (if `weight = TRUE`).

## Usage

``` r
net_to_mat(
  net,
  weight = FALSE,
  squared = FALSE,
  symmetrical = FALSE,
  missing_value = 0
)
```

## Arguments

- net:

  A two- or three-column `data.frame` where each row represents the
  interaction between two nodes (e.g., site and species), with an
  optional third column indicating the weight of the interaction.

- weight:

  A `logical` value indicating whether the weight column should be
  considered.

- squared:

  A `logical` value indicating whether the output matrix should be
  square (i.e., containing the same nodes in rows and columns).

- symmetrical:

  A `logical` value indicating whether the resulting matrix should be
  symmetrical. This applies only if `squared = TRUE`. Note that
  different weights associated with opposite pairs already present in
  `net` will be preserved.

- missing_value:

  The value to assign to pairs of nodes not present in `net`. Defaults
  to `0`.

## Value

A `matrix` with the first nodes (from the first column of `net`) as rows
and the second nodes (from the second column of `net`) as columns. If
`squared = TRUE`, the rows and columns will have the same number of
elements, corresponding to the unique union of objects in the first and
second columns of `net`. If `squared = TRUE` and `symmetrical = TRUE`,
the matrix will be forced to be symmetrical based on the upper
triangular part of the matrix.

## See also

For more details illustrated with a practical example, see the vignette:
<https://biorgeo.github.io/bioregion/articles/a2_matrix_and_network_formats.html>.

Associated functions:
[mat_to_net](https://bioRgeo.github.io/bioregion/reference/mat_to_net.md)

## Author

Maxime Lenormand (<maxime.lenormand@inrae.fr>)  
Pierre Denelle (<pierre.denelle@gmail.com>)  
Boris Leroy (<leroy.boris@gmail.com>)

## Examples

``` r
net <- data.frame(
  Site = c(rep("A", 2), rep("B", 3), rep("C", 2)),
  Species = c("a", "b", "a", "c", "d", "b", "d"),
  Weight = c(10, 100, 1, 20, 50, 10, 20)
)

mat <- net_to_mat(net, weight = TRUE)
```
