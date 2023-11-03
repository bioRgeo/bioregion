#' Divisive hierarchical clustering based on dissimilarity or beta-diversity
#'
#' This function computes a divisive hierarchical clustering from a
#' dissimilarity (beta-diversity) `data.frame`, calculates the cophenetic correlation
#' coefficient, and can get clusters from the tree if requested by the user.
#' The function implements randomization of the dissimilarity matrix to
#' generate the tree, with a selection method based on the optimal cophenetic
#' correlation coefficient. Typically, the dissimilarity `data.frame` is a
#' `bioregion.pairwise.metric` object obtained by running `similarity`
#' or `similarity` and then `similarity_to_dissimilarity`.
#'
#' @param dissimilarity the output object from [dissimilarity()] or
#'  [similarity_to_dissimilarity()], or a `dist` object. 
#'  If a `data.frame` is used, the first two 
#' columns represent pairs of sites (or any pair of nodes), and the next
#' column(s) are the dissimilarity indices.
#' 
#' @param index name or number of the dissimilarity column to use. By default, 
#' the third column name of `dissimilarity` is used.
#'  
#' @param n_clust an integer or a vector of integers indicating the number of
#' clusters to be obtained from the hierarchical tree, or the output from
#' [partition_metrics]. Should not be used at the same time as
#' `cut_height`.
#' 
#' @param cut_height a numeric vector indicating the height(s) at which the
#' tree should be cut. Should not be used at the same time as `n_clust`.
#' 
#' @param find_h a boolean indicating if the height of cut should be found for
#' the requested `n_clust`.
#' 
#' @param h_max a numeric indicating the maximum possible tree height for
#' the chosen `index`.
#' 
#' @param h_min a numeric indicating the minimum possible height in the tree
#' for the chosen `index`.
#' 
#' @details
#' Chapter 6 of Kaufman and Rousseeuw (1990) fully details the functioning of
#' the diana algorithm.
#'
#' To find an optimal number of clusters, see [partition_metrics()]
#'
#' @return
#' A `list` of class `bioregion.clusters` with five slots:
#' \enumerate{
#' \item{**name**: `character string` containing the name of the algorithm}
#' \item{**args**: `list` of input arguments as provided by the user}
#' \item{**inputs**: `list` of characteristics of the clustering process}
#' \item{**algorithm**: `list` of all objects associated with the
#'  clustering procedure, such as original cluster objects}
#' \item{**clusters**: `data.frame` containing the clustering results}}#'
#'
#' @references
#' Kaufman and Rousseeuw (1990)
#'
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}),
#' Boris Leroy (\email{leroy.boris@gmail.com}) and
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) 
#' 
#' @seealso [cut_tree] 
#' 
#' @examples
#' comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
#' 20, 25)
#' rownames(comat) <- paste0("Site",1:20)
#' colnames(comat) <- paste0("Species",1:25)
#'
#' dissim <- dissimilarity(comat, metric = "all")
#'
#' data("fishmat")
#' fishdissim <- dissimilarity(fishmat)
#' fish_diana <- hclu_diana(fishdissim, index = "Simpson")
#' 
#' # fish_diana <- cluster::diana(fishmat, #dissim,
#' #                              diss = FALSE, #inherits(fishdissim, "dist"),
#' #                              metric = "euclidean",
#' #                              stand = FALSE,
#' #                              stop.at.k = FALSE,
#' #                              # keep.diss = n < 100,
#' #                              # keep.data = !diss,
#' #                              trace.lev = 0)
#' 
#' 
#' @importFrom cluster diana
#' @importFrom stats as.dist cophenetic cor
#' 
#' @export

hclu_diana <- function(dissimilarity,
                       index = names(dissimilarity)[3],
                       n_clust = NULL,
                       cut_height = NULL,
                       find_h = TRUE,
                       h_max = 1,
                       h_min = 0){
  
  # 1. Controls ---------------------------------------------------------------
  if(inherits(dissimilarity, "bioregion.pairwise.metric")){
    if(attr(dissimilarity, "type") == "similarity") {
      stop("dissimilarity seems to be a similarity object.
         hclu_hierarclust() should be applied on dissimilarity, not
         similarities.
         Use similarity_to_dissimilarity() before using hclu_hierarclust()")
    }
    if(is.numeric(index)){
      index <- names(dissimilarity)[index]
    }
    if(!(index %in% colnames(dissimilarity))) {
      stop("Argument index should be one of the column names of dissimilarity")
    }
    
  } else if(!any(inherits(dissimilarity, "bioregion.pairwise.metric"),
                 inherits(dissimilarity, "dist"))){
    if(is.numeric(index)) {
      index <- names(dissimilarity)[index]
    } 
    if(is.null(index) || !(index %in% colnames(dissimilarity))) {
      stop("dissimilarity is not a bioregion.pairwise.metric object, a
           dissimilarity matrix (class dist) or a data.frame with at least 3
           columns (site1, site2, and your dissimilarity index).")
    }
  }
  
  if(!is.null(n_clust)) {
    if(is.numeric(n_clust)) {
      if(any(!(n_clust %% 1 == 0))) {
        stop("n_clust must an integer or a vector of integers determining the
             number of clusters.")
      }
    } else if(inherits(n_clust, "bioregion.partition.metrics")){
      if(!is.null(n_clust$algorithm$optimal_nb_clusters)) {
        n_clust <- n_clust$algorithm$optimal_nb_clusters
      } else {
        stop("n_clust does not have an optimal number of clusters. Did you
        specify partition_optimisation = TRUE in partition_metrics()?")
      }
    } else{
      stop("n_clust must be one of those:
        * an integer determining the number of clusters
        * a vector of integers determining the numbers of clusters for each cut
        * the output from partition_metrics()")
    }
    if(!is.null(cut_height)){
      stop("Please provide either n_clust or cut_height, but not both at the
           same time.")
    }
  }
  
  if(!is.null(cut_height)){
    if(!is.numeric(cut_height) || any(cut_height < 0)){
      stop("cut_height must be a positive integer.")
    }
  }
  
  if(!is.logical(find_h)){
    stop("find_h must be a Boolean.")
  }
  
  if(!is.numeric(h_max) || h_max < 0){
    stop("h_max must be a positive integer.")
  }
  
  if(!is.numeric(h_min) || h_min < 0){
    stop("h_min must be a positive integer.")
  }
  
  if(h_min > h_max){
    stop("h_min must be inferior to h_max.")
  }
  
  # 2. Function ---------------------------------------------------------------
  # Convert dissimilarity into a dist object
  if(!inherits(dissimilarity, "dist")){
    # dist.obj <- .dfToDist(dissimilarity, metric = index)
    dist.obj <- stats::as.dist(
      net_to_mat(dissimilarity[, c(colnames(dissimilarity)[1:2], index)],
                 weight = TRUE, squared = TRUE, symmetrical = TRUE))
    
  } else {
    dist.obj <- dissimilarity
  }
  
  # Output of the function
  outputs <- list(name = "divisive_hierarchical_clustering")
  
  # Adding dynamic_tree_cut = FALSE for compatibility with generic functions
  dynamic_tree_cut <- FALSE
  outputs$args <- list(index = index,
                       # method = method,
                       # randomize = randomize,
                       # n_runs = n_runs,
                       # optimal_tree_method = optimal_tree_method,
                       n_clust = n_clust,
                       cut_height = cut_height,
                       find_h = find_h,
                       h_max = h_max,
                       h_min = h_min,
                       dynamic_tree_cut = dynamic_tree_cut)
  
  outputs$inputs <- list(nb_sites = attr(dist.obj, "Size"))
  
  # DIANA clustering
  diana_clust <- cluster::diana(dist.obj,
                                diss = inherits(dist.obj, "dist"),
                                # metric = "euclidean",
                                # stand = FALSE,
                                # stop.at.k = FALSE,
                                # keep.diss = n < 100,
                                # keep.data = !diss,
                                trace.lev = 0)
  
  outputs$algorithm$final.tree <- diana_clust
  # outputs$diana <- diana_clust
  
  # Cophenetic correlation
  coph <- as.matrix(stats::cophenetic(outputs$algorithm$final.tree))
  coph <- coph[match(attr(dist.obj, "Labels"), rownames(coph)),
               match(attr(dist.obj, "Labels"), colnames(coph))]
  dist.mat <- as.matrix(dist.obj)
  
  outputs$algorithm$final.tree.coph.cor <-
    stats::cor(dist.mat[lower.tri(dist.mat)], coph[lower.tri(coph)],
               method = "spearman")
  
  message(paste0("Output tree has a ",
                 round(outputs$algorithm$final.tree.coph.cor, 2),
                 " cophenetic correlation coefficient with the initial
                   dissimilarity matrix\n"))
  
  class(outputs) <- append("bioregion.clusters", class(outputs))
  
  # Cut tree
  if(any(!is.null(n_clust) | !is.null(cut_height))){
    outputs$clusters <- cut_tree(stats::as.hclust(outputs$algorithm$final.tree), # outputs,
                                 n_clust = n_clust,
                                 cut_height = cut_height,
                                 find_h = find_h,
                                 h_max = h_max,
                                 h_min = h_min,
                                 dynamic_tree_cut = dynamic_tree_cut)
    
    outputs$cluster_info <- data.frame(
      partition_name = names(outputs$clusters)[2:length(outputs$clusters),
                                               drop = FALSE],
      n_clust = length(unique(outputs$clusters[, 2])))
    
    outputs$inputs$hierarchical <- ifelse(ncol(outputs$clusters) > 2,
                                          TRUE,
                                          FALSE)
  } else {
    outputs$clusters <- NA
    outputs$cluster_info <- NA
    outputs$inputs$hierarchical <- FALSE
  }
  
  return(outputs)
  
}
