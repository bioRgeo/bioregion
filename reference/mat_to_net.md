# Create a data.frame from a contingency table

This function generates a two- or three-column `data.frame`, where each
row represents the interaction between two nodes (e.g., site and
species) and an optional third column indicates the weight of the
interaction (if `weight = TRUE`). The input is a contingency table, with
rows representing one set of entities (e.g., site) and columns
representing another set (e.g., species).

## Usage

``` r
mat_to_net(
  mat,
  weight = FALSE,
  remove_zeroes = TRUE,
  include_diag = TRUE,
  include_lower = TRUE
)
```

## Arguments

- mat:

  A contingency table (i.e., a `matrix`).

- weight:

  A `logical` value indicating whether the values in the matrix should
  be interpreted as interaction weights.

- remove_zeroes:

  A `logical` value determining whether interactions with a weight equal
  to 0 should be excluded from the output.

- include_diag:

  A `logical` value indicating whether the diagonal (self-interactions)
  should be included in the output. This applies only to square
  matrices.

- include_lower:

  A `logical` value indicating whether the lower triangular part of the
  `matrix` should be included in the output. This applies only to square
  matrices.

## Value

A `data.frame` where each row represents the interaction between two
nodes. If `weight = TRUE`, the `data.frame` includes a third column
representing the weight of each interaction.

## See also

For more details illustrated with a practical example, see the vignette:
<https://biorgeo.github.io/bioregion/articles/a2_matrix_and_network_formats.html>.

Associated functions:
[net_to_mat](https://bioRgeo.github.io/bioregion/reference/net_to_mat.md)

## Author

Maxime Lenormand (<maxime.lenormand@inrae.fr>)  
Pierre Denelle (<pierre.denelle@gmail.com>)  
Boris Leroy (<leroy.boris@gmail.com>)

## Examples

``` r
mat <- matrix(sample(1000, 50), 5, 10)
rownames(mat) <- paste0("Site", 1:5)
colnames(mat) <- paste0("Species", 1:10)

net <- mat_to_net(mat, weight = TRUE)
```
