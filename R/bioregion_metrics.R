#' Calculate contribution metrics for bioregions
#' 
#' This function calculates the number of sites per bioregion, as well as the
#' the number of species these sites have, the number of endemic species and 
#' the proportion of endemism.
#' 
#' @param cluster_object a `bioregion.clusters` object or a `data.frame` or a 
#' list of `data.frame` containing multiple partitions. At least two partitions
#' are required. If a list of `data.frame` is provided, they should all have
#' the same number of rows (i.e., same items in the clustering for all
#' partitions). 
#' 
#' @param comat a co-occurrence `matrix` with sites as rows and species as
#' columns. 
#' 
#' @details 
#' Endemic species are species only found in the sites belonging to one
#' bioregion.
#'
#' @seealso [site_species_metrics]
#' @return 
#' A `data.frame` with 5 columns.
#'  
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}),
#' Boris Leroy (\email{leroy.boris@gmail.com}) and
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) 
#' 
#' @examples
#' comat_1 <- matrix(sample(0:1000, size = 10*12, replace = TRUE,
#' prob = 1/1:1001), 10, 12)
#' rownames(comat_1) <- paste0("Site", 1:10)
#' colnames(comat_1) <- paste0("Species", 1:12)
#' comat_1 <- cbind(comat_1,
#'                  matrix(0, 10, 8,
#'                         dimnames = list(paste0("Site", 1:10),
#'                                         paste0("Species", 13:20))))
#' 
#' comat_2 <- matrix(sample(0:1000, size = 10*12, replace = TRUE,
#'                          prob = 1/1:1001), 10, 12)
#' rownames(comat_2) <- paste0("Site", 11:20)
#' colnames(comat_2) <- paste0("Species", 9:20)
#' comat_2 <- cbind(matrix(0, 10, 8,
#'                         dimnames = list(paste0("Site", 11:20),
#'                                         paste0("Species", 1:8))),
#'                  comat_2)
#' 
#' comat <- rbind(comat_1, comat_2)
#' 
#' dissim <- dissimilarity(comat, metric = "Simpson")
#' clust1 <- nhclu_kmeans(dissim, n_clust = 3, index = "Simpson")
#' 
#' net <- similarity(comat, metric = "Simpson")
#' com <- netclu_greedy(net)
#' 
#' bioregion_metrics(cluster_object = clust1, comat = comat) 
#' 
#' @export

bioregion_metrics <- function(cluster_object, comat){
  # 1. Controls ---------------------------------------------------------------
  # input can be of format bioregion.clusters
  if (inherits(cluster_object, "bioregion.clusters")) {
    if (inherits(cluster_object$clusters, "data.frame")) {
      has.clusters <- TRUE
      clusters <- cluster_object$clusters
      
      if(ncol(clusters) > 2) {
        stop("This function is designed to be applied on a single partition.",
             "Your cluster_object has multiple partitions (select only one).")
      }
      
    } else {
      if (cluster_object$name == "hierarchical_clustering") {
        stop("No clusters have been generated for your hierarchical tree,
        please extract clusters from the tree before using partition_metrics()
        See ?hclu_hierarclust or ?cut_tree")
      } else {
        stop(
          "cluster_object does not have the expected type of 'clusters' slot")
      }
    }
  } else {
    stop("This function is designed to work on bioregion.clusters objects and
         on a site x species matrix.")
  }
  
  controls(args = NULL, data = comat, type = "input_matrix")
  
  bioregion_df <- NULL
  
  # 2. Function ---------------------------------------------------------------
  # clusters <- clust1$clusters
  # tmp <- unique(clusters[which(clusters[, 2] == 1), "ID"])
  # colnames(comat[tmp, ])
  # dim(comat); dim(comat[tmp, ]); dim(comat[tmp, colSums(comat[tmp, ]) > 0])
  
  # If it is a bipartite object, we can directly count the species
  if(cluster_object$inputs$bipartite == TRUE){
    clusters <- clusters[which(attributes(clusters)$node_type == "site"), ]
  }
  
  bioregion_df <- data.frame()
  
  # Loop over bioregions
  for(j in 1:cluster_object$cluster_info$n_clust){
    focal_j <- unique(clusters[, 2])[j] # bioregion j
    
    # Sites belonging to bioregion j
    sites_j <- clusters[which(clusters[, 2] == focal_j), "ID"]
    
    # Subset site x species matrix with sites belonging to bioregion j
    comat_j <- comat[sites_j, colSums(comat[sites_j, ]) > 0]
    
    # Other sites
    comat_not_j <- comat[-which(rownames(comat) %in% sites_j), ]
    comat_not_j <- comat_not_j[, colSums(comat_not_j) > 0]
    
    # Number of endemics
    endemic_j <- sum(!(colnames(comat_j) %in% colnames(comat_not_j)))
    
    # Output
    bioregion_df <-
      rbind(bioregion_df,
            data.frame(Bioregion = focal_j,
                       Site_number = length(sites_j),
                       Species_number = ncol(comat_j),
                       Endemics = endemic_j,
                       Percentage_Endemic = 100*endemic_j/ncol(comat_j)))
  }
  
  # Test if all bioregions are in the output
  if(length(unique(bioregion_df$Bioregion)) !=
     cluster_object$cluster_info$n_clust){
    warning("Not all bioregions are in the output.")
  }
  
  return(bioregion_df)    
}
