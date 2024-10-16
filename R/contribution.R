#' Calculate contribution metrics of sites and species
#' 
#' Cz metrics
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
#' @param indices a `character` specifying the contribution metric to compute.
#' Available options are `contribution`.
#' 
#' @details 
#' The contribution metric is derived from \insertRef{Lenormand2019}{bioregion}.
#' Its formula is the following:
#' \eqn{(n_ij - ((n_i n_j)/n))/(sqrt(((n - n_j)/(n-1)) (1-(n_j/n)) ((n_i n_j)/n)))}
#' 
#' with n the number of sites, n_i the number of sites in which
#' species i is present, n_j the number of sites belonging to the
#' bioregion j, n_ij the number of occurrences of species i in sites
#' belonging to the bioregion j.
#'
#' @seealso [partition_metrics]
#' @return 
#' A `list` of `data.frames` if multiples indices are selected or a single
#' `data.frame` with three columns if one index is selected. Each `data.frame`
#' has three columns: the species, the bioregion, and the contribution
#' statistics.
#'  
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}),
#' Boris Leroy (\email{leroy.boris@gmail.com}) and
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) 
#' 
#' @examples
#' comat <- matrix(sample(1000, 50), 5, 10)
#' rownames(comat) <- paste0("Site", 1:5)
#' colnames(comat) <- paste0("Species", 1:10)
#' 
#' comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
#'                 20, 25)
#' rownames(comat) <- paste0("Site",1:20)
#' colnames(comat) <- paste0("Species",1:25)
#' 
#' dissim <- dissimilarity(comat, metric = "Simpson")
#' clust1 <- nhclu_kmeans(dissim, n_clust = 3, index = "Simpson")
#' 
#' net <- similarity(comat, metric = "Simpson")
#' com <- netclu_greedy(net)
#' 
#' contribution(cluster_object = clust1, comat = comat,
#' indices = "contribution")
#' 
#' contribution(cluster_object = com, comat = comat,
#' indices = "contribution")
#' 
#' @export

contribution <- function(cluster_object,
                         comat,
                         indices = c("contribution", "Cz")){
  # 1. Controls ---------------------------------------------------------------
  # input can be of format bioregion.clusters
  if (inherits(cluster_object, "bioregion.clusters")) {
    if (inherits(cluster_object$clusters, "data.frame")) {
      has.clusters <- TRUE
      clusters <- cluster_object$clusters
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
  
  if(ncol(clusters) > 2) {
    stop("This function is designed to be applied on a single partition.",
         "Your cluster_object has multiple partitions (select only one).")
  }
  
  controls(args = NULL, data = comat, type = "input_matrix")
  
  controls(args = indices, data = NULL, type = "character")
  if(!(indices %in% c("contribution", "Cz"))){
    stop("Please choose algorithm among the followings values:
    contribution or Cz.", call. = FALSE)
  }
  
  contribution_df <- NULL
  
  # 2. Function ---------------------------------------------------------------
  ## 2.1. Cz ------------------------------------------------------------------
  # only for bipartite cases; not implemented yet
  if("Cz" %in% indices){
    message("Cz indices not yet implemented.")
    Cz_df <- data.frame()
  }
  
  ## 2.2. Contribution --------------------------------------------------------
  if("contribution" %in% indices){
    # Binary site-species matrix
    comat_bin <- comat
    comat_bin[comat_bin > 0] <- 1
    
    # Formula
    n <- nrow(comat) # number of sites
    n_i <- colSums(comat_bin) # number of occurrences per species
    n_j <- table(cluster_object$clusters[, 2]) # number of sites per bioregion
    
    contrib_df <- data.frame()
    
    # Loop over bioregions
    for(j in 1:cluster_object$cluster_info$n_clust){
      focal_j <- unique(cluster_object$clusters[, 2])[j] # bioregion j
      
      # Sites belonging to bioregion j
      focal_sites <- cluster_object$clusters[which(
        cluster_object$clusters[, 2] == focal_j), 1]
      # Number of sites belonging to bioregion j
      n_j <- table(cluster_object$clusters[, 2])[[focal_j]]
      
      # Occurrences per species in each of these sites to get n_ij
      n_ij <- colSums(comat[focal_sites, , drop = FALSE])
      
      # Contribution of species i to bioregion j
      p_ij <- (n_ij - ((n_i*n_j)/n))/(sqrt(((n - n_j)/
                                              (n-1))*(1-(n_j/n))*((n_i*n_j)/n)))
      
      contrib_df <- rbind(contrib_df,
                          data.frame(Bioregion = focal_j,
                                     Species = names(p_ij),
                                     Contribution = as.numeric(p_ij)))
    }
    
    # Controls on the output
    # test if all bioregions are there
    if(length(unique(contrib_df$Bioregion)) !=
       cluster_object$cluster_info$n_clust){
      warning("Not all bioregions are in the output.")
    }
    
    # test if all species are there
    if(length(unique(contrib_df$Species)) != ncol(comat)){
      warning("Not all species are in the output.")
    }
    
    # test if all species are there X times (X = nb of bioregions)
    if(length(unique(table(contrib_df$Species))) != 1 ||
       unique(table(contrib_df$Species)) != cluster_object$cluster_info$n_clust){
      warning("Not all species x bioregions combinations are in the output.")
    }
  }
  
  if(length(indices) == 1){
    if(indices == "Cz"){
      return(Cz_df)
    } else if(indices == "contribution"){
      return(contrib_df)    
    }
  } else{
    if(indices == c("contribution", "Cz") || indices == c("Cz", "contribution")){
      return(list(contribution = contribution_df,
                  cz = Cz_df))
    }
  }
}
