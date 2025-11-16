**<span style="color:black"><font size="4">bioregion's TO DO list</span></font>**

#### Issues (in French, some of them may have already been fixed).

bioregion_metrics
    - Problème de contradiction entre la documentation et la fonction pour cluster_object
      (peut-on utiliser un data.frame ou une liste de data.frame ?) ?
    - bug     comat_not_j <- comat[-which(rownames(comat) %in% sites_j), ] quand
      un seul site dans une bioregion j'ai l'impression
#' comat <- matrix(sample(1000, 50), 5, 10)
#' rownames(comat) <- paste0("Site", 1:5)
#' colnames(comat) <- paste0("Species", 1:10)
#'
#' net <- similarity(comat, metric = "Simpson")
#' clust <- netclu_greedy(net)
#' 
#' bioregion_metrics(bioregionalization = clust, 
#'                   comat = comat) 

bioregionalizations_metrics
    - bioregion.partition.metrics ? On utilise les termes cluster_object, clusters, partition,
      partitions, bioregionalizations de manière interchangeable. Faut-il harmoniser les noms ?
    - Ajouter data.frame et dist dans les entrées possibles.
    - Documentation étrange pour dissimilarity (tree ?)

find_optimal_n
    - Changer le nom partitions ?
    - Le data.frame avec les colonnes K et n_clusters trop contraignant ?
    - Test valeurs NA (et d'autres tests) pour le data.frame ?
    - Problème avec les valeurs min et max (voir tests) ?

compare_bioregionalizations
    - Ajouter un test pour vérifier les data.frame et adapter les contrôles pour les listes de
      data.frame ?
    - Ajouter une colonne ID pour chaque data.frame ?
    - Imaginer un méta-objet bioregion.clusters permettant de fusionner les sorties des
      algorithmes ?
      
#### New functions

* Add fuzzy clustering methods

* reassign_cluster to assign and reassign a cluster/bioregion to site and species
based on site_species_metrics() [modify OSLOM] + add inputs$reassign in clusters.

* as_bioregion_cluster() 

* check_cluster() ? to check bioregionalization, comat & (dis)similarity

* bind_cluster()

* between_bioregions_metrics() => 'beta-diversity' between bioregions, output:
matrix/network, produces results like in Lenormand et al. 2019

* Function to visualize network (maybe with visNetwork()) from
between_bioregions_metrics()

#### New features

* Allow to launch OPTICS with a vector of parameter values (as in kmeans).

* Parallelize clu_ functions. 

* https://stats.stackexchange.com/questions/260487/adjusted-rand-index-vs-adjusted-mutual-information

* Code IHCT in Rcpp

* Check hierarchy in igraph ?

#### Miscellaneous

* [Terminology] indices <-> metrics in site_species_metrics()

* [Terminology] clusters <-> bioregionalization in fonctions "Bioregionalization analysis" and
utils and visualization

* [Terminology] define species clusters and species clustering

* Add something abount NA for site/species without bioregion/cluster in vignette
and/or doc

* Check prints in general 

* Check number of clusters in cluster_info and name columns


