# bioregion 1.1.3

This is a list of changes made between bioregion 1.1.2 
(CRAN unreleased) and bioregion 1.1.3 (CRAN unreleased).

* Added a new method to construct a consensus tree in hclu_hierarclust. This
method is called Iterative Hierarchical Consensus Tree (IHCT). It conclusively
solves issues related to the order of sites in the distance matrix and 
building a consensus hierarchical tree with a meaningful topology.

* Many changes to functions related to hclu_hierarclust due to this major change

* Updated generic functions to provide plot and print for diana

* Function contribution() renamed as site_species_metrics()

* Function bioregion_metrics() created

* Function subset_node() renamed site_species_subset()

* Indices Cz added to site_species_metrics()

* Updated install_binaries:
   - Archive bin.zip now stored on GitHub and backed up on NextCloud  
   - Added Infomap version 2.8.0
   - Added argument download_only to only execute the download step  
&nbsp;
* Added argument check_install in netclu_infomap, netclu_louvain, and netclu_oslom


# bioregion 1.1.2

This is a list of changes made between bioregion 1.1.1 
(CRAN release 2024-04-19) and bioregion 1.1.2 (CRAN unreleased).

* Affinity propagation algorithm added (function nhclu_affprop)

* Function contribution() added to the package and the workflow.

# bioregion 1.1.1

This is a list of changes made between bioregion 1.1.0 
(CRAN release 2024-03-19) and bioregion 1.1.1 (CRAN release 2024-04-19).

* Added hierarchy for Louvain cpp.

* Added seed argument to stochastic algorithms (except Louvain cpp).

* Added argument cut_weight in netclu_ fonctions.

* Changed value for sites without cluster (0 -> NA).

* Updated automated tests (code coverage > 60%).

* Controls and outputs/inputs standardization.
 
* Fixed a bug in find_optimal_n() in the special case where partition 
metrics did not vary.

# bioregion 1.1.0 

This is a list of changes made between bioregion 1.0.0 
(CRAN release 2023-04-15) and bioregion 1.1.0 (CRAN release 2024-03-19).

* Added the resolution parameter in the igraph Louvain version.

* Added possibility to remove diagonal and lower triangular matrix in 
mat_to_net() for squared matrix with argument include_diag and include_lower.

* Added a function to extract a subset of node according to its type (sites or 
species) from a bioregion.clusters object containing both types of nodes (sites 
and species).

* Added a generic function to maintain attributes of bioregion.pairwise.metric
objects + keep track of number of sites and species.

* Functions added: nhclu_clara(), nhclu_clarans().  

* The corresponding vignettes are edited to document the new functions.  

* Modification of the way 'bioregion.pairwise.metric' object are controlled.

* Allow to (not) select 'formula metrics' in 
similarity_dissimilarity_conversion() with the new argument "include_formula".

* Allow negative values in similarity() with the Euclidean metric.

# bioregion 1.0.0 

First release on CRAN

