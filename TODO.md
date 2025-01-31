**<span style="color:black"><font size="4">bioregion's TO DO list</span></font>**

#### Issues (in French, some of them may have already been fixed).

nhbclu_clara et clarans
    - bien vérifier le problème de seed...

site_species_metrics
    - Ajouter un test pour bipartite_link ?
    - Changer le nom bipartite_link (net dans bioregionalization_metrics) ?
    - Problèmes avec affinity, fidelity et indicator_value :
        - Plantage lorsqu'ils sont utilisés seuls (problème avec names(p_ij) non défini).
        - Seulement rho et cz spécifié dans la documentation.
    - Problème de contradiction entre la documentation et la fonction pour cluster_object
      (peut-on utiliser un data.frame ou une liste de data.frame ?) ?

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

* from_abcABC (to compute metrics from abc and/or ABC).

* reassign_cluster (to assign or reassign cluster based on a partition and a 
(dis)similarity metric) + spatial coherence (reassignment of sites)

* between_bioregions_metrics() => 'beta-diversity' between bioregions, output:
matrix/network, produces results like in Lenormand et al. 2019

* Function to visualize network (maybe with visNetwork()) from
between_bioregions_metrics()

#### New features

* Allow to launch OPTICS with a vector of parameter values (as in kmeans).

* Parallelize clu_ functions. 

* https://stats.stackexchange.com/questions/260487/adjusted-rand-index-vs-adjusted-mutual-information

* Code IHCT in Rcpp

#### Miscellaneous

* data-raw ? 


