#' Non hierarchical clustering: partitioning around medoids
#'
#' This function performs non hierarchical clustering on the basis of
#' dissimilarity with partitioning around medoids.
#'
#' @param dissimilarity the output object from [dissimilarity()] or
#' [similarity_to_dissimilarity()], or a `dist` object. If a `data.frame` is
#' used, the first two columns represent pairs of sites (or any pair of nodes),
#' and the next column(s) are the dissimilarity indices.
#' 
#' @param index name or number of the dissimilarity column to use. By default, 
#' the third column name of `dissimilarity` is used.
#' 
#' @param seed for the random number generator (NULL for random by default).
#' 
#' @param n_clust an `integer` vector or a single `integer` value specifying the
#' requested number(s) of clusters.
#' 
#' @param variant a `character` string specifying the variant of pam to use,
#' by default `faster`. Available options are `original`, `o_1`, `o_2`, `f_3`, 
#' `f_4`, `f_5` or `faster`. See [pam][cluster::pam] for more details.
#' 
#' @param nstart an `integer` specifying the number of random start for the
#' pam algorithm. By default, 1 (for the `faster` variant).
#' 
#' @param cluster_only a `boolean` specifying if only the clustering should be
#' returned from the [pam][cluster::pam] function (more efficient).
#' 
#' @param algorithm_in_output a `boolean` indicating if the original output
#' of [pam][cluster::pam] should be returned in the output (`TRUE` by 
#' default, see Value).
#' 
#' @param ... you can add here further arguments to be passed to `pam()`
#' (see [pam][cluster::pam]).
#'
#' @details
#' This method partitions data into the chosen number of cluster on the basis
#' of the input dissimilarity matrix. It is more robust than k-means because it
#' minimizes the sum of dissimilarity between cluster centres and points
#' assigned to the cluster - whereas the k-means approach minimizes the sum of
#' squared euclidean distances (thus k-means cannot be applied directly on the
#' input dissimilarity matrix if the distances are not euclidean).
#'
#' @return
#' A `list` of class `bioregion.clusters` with five slots:
#' \enumerate{
#' \item{**name**: `character` containing the name of the algorithm}
#' \item{**args**: `list` of input arguments as provided by the user}
#' \item{**inputs**: `list` of characteristics of the clustering process}
#' \item{**algorithm**: `list` of all objects associated with the
#'  clustering procedure, such as original cluster objects}
#' \item{**clusters**: `data.frame` containing the clustering results}}
#' 
#' In the `algorithm` slot, if `algorithm_in_output = TRUE`, users can
#' find the output of
#' [pam][cluster::pam].
#' 
#' @references
#' Kaufman L & Rousseeuw PJ (2009) Finding groups in data: An introduction to 
#' cluster analysis. In & Sons. JW (ed.), Finding groups in data: An 
#' introduction to cluster analysis.
#' 
#' @author
#' Boris Leroy (\email{leroy.boris@gmail.com}) \cr
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) \cr
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) 
#' 
#' @seealso [nhclu_kmeans] 
#' 
#' @examples
#' comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
#' 20, 25)
#' rownames(comat) <- paste0("Site",1:20)
#' colnames(comat) <- paste0("Species",1:25)
#'
#' comnet <- mat_to_net(comat)
#' dissim <- dissimilarity(comat, metric = "all")
#' 
#' clust1 <- nhclu_pam(dissim, n_clust = 2:10, index = "Simpson")
#' clust2 <- nhclu_pam(dissim, n_clust = 2:15, index = "Simpson")
#' bioregionalization_metrics(clust2, dissimilarity = dissim,
#' eval_metric = "pc_distance")
#' bioregionalization_metrics(clust2, net = comnet, species_col = "Node2",
#'                    site_col = "Node1", eval_metric = "avg_endemism")
#'    
#' @importFrom stats as.dist
#' @importFrom cluster pam    
#'                    
#' @export

nhclu_pam <- function(dissimilarity,
                      index = names(dissimilarity)[3],
                      seed = NULL,
                      n_clust = c(1,2,3),
                      variant = "faster", 
                      nstart = 1,
                      cluster_only = FALSE, 
                      algorithm_in_output = TRUE,
                      ...){ 
  
  # 1. Controls ---------------------------------------------------------------
  controls(args = NULL, data = dissimilarity, type = "input_nhandhclu")
  if(!inherits(dissimilarity, "dist")){
    controls(args = NULL, data = dissimilarity, type = "input_dissimilarity")
    controls(args = NULL, data = dissimilarity, 
             type = "input_data_frame_nhandhclu")
    controls(args = index, data = dissimilarity, type = "input_net_index")
    net <- dissimilarity
    # Convert tibble into dataframe
    if(inherits(net, "tbl_df")){
      net <- as.data.frame(net)
    }
    net[, 3] <- net[, index]
    net <- net[, 1:3]
    controls(args = NULL, data = net, type = "input_net_index_value")
    dist.obj <- stats::as.dist(
      net_to_mat(net,
                 weight = TRUE, squared = TRUE, symmetrical = TRUE))
  } else {
    controls(args = NULL, data = dissimilarity, type = "input_dist")
    dist.obj <- dissimilarity
    if(is.null(labels(dist.obj))){
      attr(dist.obj, "Labels") <- paste0(1:attr(dist.obj, "Size"))
      message("No labels detected, they have been assigned automatically.")
    }
  }
  
  if(!is.null(seed)){
    controls(args = seed, data = NULL, type = "strict_positive_integer")
  }
  controls(args = n_clust, data = NULL, 
           type = "strict_positive_integer_vector")
  controls(args = variant, data = NULL, type = "character")
  if(!(variant %in% c("original", "o_1", "o_2", "f_3", "f_4", "f_5",
                          "faster"))){
    stop("Please choose variant among the followings values:
original, o_1, o_2, f_3, f_4, f_5 or faster", call. = FALSE)
  }
  controls(args = nstart, data = NULL, type = "strict_positive_integer")
  controls(args = cluster_only, data = NULL, type = "boolean")
  controls(args = algorithm_in_output, data = NULL, type = "boolean")

  # 2. Function ---------------------------------------------------------------
  outputs <- list(name = "nhclu_pam")
  
  outputs$args <- list(index = index,
                       seed = seed,
                       n_clust = n_clust,
                       nstart = nstart,
                       variant = variant,
                       cluster_only = cluster_only,
                       algorithm_in_output = algorithm_in_output,
                       ...)
  
  outputs$inputs <- list(bipartite = FALSE,
                         weight = TRUE,
                         pairwise = TRUE,
                         pairwise_metric = ifelse(!inherits(dissimilarity, 
                                                            "dist"), 
                                                  ifelse(is.numeric(index), 
                                                         names(net)[3], index), 
                                                  NA),
                         dissimilarity = TRUE,
                         nb_sites = attr(dist.obj, "Size"),
                         hierarchical = FALSE)
  
  outputs$algorithm <- list()
  
  outputs$clusters <- data.frame(matrix(ncol = 1,
                                        nrow = length(labels(dist.obj)),
                                        dimnames = list(labels(dist.obj),
                                                        "name")))
  
  outputs$clusters$name <- labels(dist.obj)
  
  if(is.null(seed)){
    outputs$algorithm <- lapply(n_clust,
                                function(x)
                                  cluster::pam(dist.obj,
                                               k = x,
                                               diss = TRUE,
                                               keep.diss = FALSE,
                                               keep.data = FALSE,
                                               nstart = nstart,
                                               variant = variant,
                                               cluster.only = cluster_only,
                                               ...))
  }else{
    set.seed(seed)
    outputs$algorithm <- lapply(n_clust,
                                function(x)
                                  cluster::pam(dist.obj,
                                               k = x,
                                               diss = TRUE,
                                               keep.diss = FALSE,
                                               keep.data = FALSE,
                                               nstart = nstart,
                                               variant = variant,
                                               cluster.only = cluster_only,
                                               ...))
    rm(.Random.seed, envir=globalenv())
  }
  
  names(outputs$algorithm) <- paste0("K_", n_clust)
  
  outputs$clusters <- data.frame(
    outputs$clusters,
    data.frame(lapply(names(outputs$algorithm),
                      function(x) outputs$algorithm[[x]]$clustering)))
  
  outputs$clusters <- knbclu(outputs$clusters, reorder = TRUE)
  
  outputs$cluster_info <- data.frame(
    partition_name = names(outputs$clusters)[2:length(outputs$clusters),
                                             drop = FALSE],
    n_clust = apply(outputs$clusters[, 2:length(outputs$clusters),
                                     drop = FALSE],
                    2, function(x) length(unique(x[!is.na(x)]))))
  
  class(outputs) <- append("bioregion.clusters", class(outputs))
  
  # Set algorithm in output
  if (!algorithm_in_output) {
    outputs$algorithm <- NA
  }
  
  return(outputs)
}
