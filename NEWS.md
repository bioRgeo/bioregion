# bioregion 1.1.1 (development/github version)

This is a list of changes made in the development/github version of the package
between bioregion 1.1.0 (CRAN release 2024-03-19) and the next CRAN release.

* Updated automated tests [objective code coverage percentage > 60%]
 
# bioregion 1.1.0 

This is a list of changes made in the development/github version of the package
between bioregion 1.0.0 (CRAN release 2023-04-15) and bioregion 1.1.0 
(CRAN release 2024-03-19).

* Added the resolution parameter in the igraph Louvain version

* Added possibility to remove diagonal and lower triangular matrix in 
mat_to_net() for squared matrix with argument include_diag and include_lower

* Added a function to extract a subset of node according to its type (sites or 
species) from a bioregion.clusters object containing both types of nodes (sites 
and species).

* Added a generic function to maintain attributes of bioregion.pairwise.metric
objects + keep track of number of sites and species

* Functions added: nhclu_clara(), nhclu_clarans()  

* The corresponding vignettes are edited to document the new functions.  

* Modification of the way 'bioregion.pairwise.metric' object are controlled.

* Allow to (not) select 'formula metrics' in 
similarity_dissimilarity_conversion() with the new argument "include_formula".

* Allow negative values in similarity() with the Euclidean metric.

# bioregion 1.0.0 

First release on CRAN

