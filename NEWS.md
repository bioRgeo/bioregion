# bioregion 1.0.3 (development/github version)

This is a list of changes made in the development/github version of the package
between bioregion 1.0.0 (CRAN release 2023-04-15) and the next CRAN release.

* Added possibility to remove diagonal and lower triangular matrix in 
mat_to_net() for squared matrix

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

