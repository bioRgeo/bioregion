#' Non hierarchical clustering: CLARANS
#'
#' This function performs non hierarchical clustering on the basis of
#' dissimilarity with partitioning around medoids, using the Clustering Large
#' Applications based on RANdomized Search (CLARANS) algorithm.
#'
#' @param dissimilarity the output object from [dissimilarity()] or
#' [similarity_to_dissimilarity()], or a `dist` object. If a `data.frame` is
#' used, the first two columns represent pairs of sites (or any pair of nodes),
#' and the next column(s) are the dissimilarity indices.
#' 
#' @param index name or number of the dissimilarity column to use. By default, 
#' the third column name of `dissimilarity` is used.
#' 
#' @param n_clust an `integer` or a `vector` of `integers` specifying the
#' requested number(s) of clusters.
#' 
#' @param numlocal an `integer` defining the number of samples to draw.
#' 
#' @param maxneighbor A positive numeric defining the sampling rate.
#' 
#' @param seed an `integer` to define a generator of random numbers.
#' 
#' @details
#' Based on fastkmedoids R package.
#'
#' @return
#' A `list` of class `bioregion.clusters` with five slots:
#' \enumerate{
#' \item{**name**: `character string` containing the name of the algorithm}
#' \item{**args**: `list` of input arguments as provided by the user}
#' \item{**inputs**: `list` of characteristics of the clustering process}
#' \item{**algorithm**: `list` of all objects associated with the
#'  clustering procedure, such as original cluster objects}
#' \item{**clusters**: `data.frame` containing the clustering results}}
#' 
#' @references
#' \insertRef{Schubert2019}{bioregion}
#' 
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}),
#' Boris Leroy (\email{leroy.boris@gmail.com}), and
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) 
#' 
#' @seealso [nhclu_pam] 
#' 
#' @examples
#' comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
#' 20, 25)
#' rownames(comat) <- paste0("Site",1:20)
#' colnames(comat) <- paste0("Species",1:25)
#'
#' dissim <- dissimilarity(comat, metric = "all")
#'
#' clust1 <- nhclu_clarans(dissim, index = "Simpson", n_clust = 5)
#' 
#' partition_metrics(clust1, dissimilarity = dissim,
#' eval_metric = "pc_distance")
#' 
#'    
#' @importFrom stats as.dist
#' @importFrom fastkmedoids fastclarans    
#'                    
#' @export

nhclu_clarans <- function(dissimilarity,
                          index = names(dissimilarity)[3],
                          n_clust = NULL,
                          numlocal = 2L,
                          maxneighbor = 0.025,
                          seed = 123456789L){
  
  # 1. Controls ---------------------------------------------------------------
  if(inherits(dissimilarity, "bioregion.pairwise.metric")){
    if(attr(dissimilarity, "type") == "similarity"){
      stop("dissimilarity seems to be a similarity object.
         nhclu_pam() should be applied on dissimilarity, not similarities.
         Use similarity_to_dissimilarity() before using nhclu_pam()")
    }
    if(is.numeric(index)){
      index <- names(dissimilarity)[index]
    }
    if(!(index %in% colnames(dissimilarity))){
      stop("Argument index should be one of the column names of dissimilarity")
    }
  } else if(!any(inherits(dissimilarity, "bioregion.pairwise.metric"),
                 inherits(dissimilarity, "dist")))
  {
    if(is.numeric(index)) {
      index <- names(dissimilarity)[index]
    }
    if(is.null(index) || !(index %in% colnames(dissimilarity))) {
      stop("dissimilarity is not a bioregion.pairwise.metric object, a
           dissimilarity matrix (class dist) or a data.frame with at least 3
           columns (site1, site2, and your dissimilarity index)")
    }
  }
  
  if(!is.null(n_clust)){
    if(is.numeric(n_clust)) {
      if(any(!(n_clust %% 1 == 0))) {
        stop("n_clust must an integer or a vector of integers determining the
             number of clusters.")
      }
      # Add test to see if n_clust is lower than the number of sites
    } else {
      stop("n_clust must an integer or a vector of integers determining the
           number of clusters.")
    }
  } else{
    stop("n_clust must an integer or a vector of integers determining the
           number of clusters.")
  }
  
  if(!inherits(dissimilarity, "dist")) {
    dist.obj <- stats::as.dist(
      net_to_mat(dissimilarity[, c(colnames(dissimilarity)[1:2], index)],
                 weight = TRUE, squared = TRUE, symmetrical = TRUE))
    
  } else {
    dist.obj <- dissimilarity
  }
  
  controls(args = numlocal, data = NULL, type = "positive_integer")
  controls(args = maxneighbor, data = NULL, type = "positive_numeric")
  controls(args = seed, data = NULL, type = "positive_integer")
  
  # 2. Function ---------------------------------------------------------------
  # Output format
  outputs <- list(name = "clarans")
  
  outputs$args <- list(index = index,
                       n_clust = n_clust,
                       numlocal = numlocal,
                       maxneighbor = maxneighbor,
                       seed = seed)
  
  outputs$inputs <- list(bipartite = FALSE,
                         weight = TRUE,
                         pairwise_metric = TRUE,
                         dissimilarity = TRUE,
                         nb_sites = attr(dist.obj, "Size"),
                         hierarchical = FALSE)
  
  outputs$algorithm <- list()
  
  outputs$clusters <- data.frame(matrix(ncol = 1,
                                        nrow = length(labels(dist.obj)),
                                        dimnames = list(labels(dist.obj),
                                                        "name")))
  
  outputs$clusters$name <- labels(dist.obj)
  
  # CLARANS algorithm
  outputs$algorithm$clarans <-
    lapply(n_clust,
           function(x)
             fastkmedoids::fastclarans(rdist = dist.obj,
                                       n = nrow(dist.obj),
                                       k = x,
                                       numlocal = numlocal,
                                       maxneighbor = maxneighbor,
                                       seed = seed))
  
  names(outputs$algorithm$clarans) <- paste0("K_", n_clust)
  
  outputs$clusters <- data.frame(
    outputs$clusters,
    data.frame(lapply(names(outputs$algorithm$clarans),
                      function(x) outputs$algorithm$clarans[[x]]@assignment)))
  
  outputs$clusters <- knbclu(outputs$clusters, reorder = FALSE)
  
  outputs$cluster_info <- data.frame(
    partition_name = names(outputs$clusters)[2:length(outputs$clusters),
                                             drop = FALSE],
    n_clust = apply(outputs$clusters[, 2:length(outputs$clusters),
                                     drop = FALSE],
                    2, function(x) length(unique(x))))
  
  class(outputs) <- append("bioregion.clusters", class(outputs))
  
  return(outputs)
}
