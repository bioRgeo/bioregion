#' Non hierarchical clustering: k-means analysis
#'
#' This function performs non hierarchical
#' clustering on the basis of dissimilarity with a k-means analysis.
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
#' requested number(s) of clusters
#' 
#' @param iter_max an `integer` specifying the maximum number of
#' iterations for the kmeans method (see [kmeans][stats::kmeans])
#' 
#' @param nstart an `integer` specifying how many random sets of
#' `n_clust` should be selected as starting points for the kmeans analysis
#' (see [kmeans][stats::kmeans])
#' 
#' @param algorithm a `character` specifying the algorithm to use for
#' kmean (see [kmeans][stats::kmeans]). Available options are
#' Hartigan-Wong, Lloyd, Forgy and MacQueen.
#' 
#' @param algorithm_in_output a `boolean` indicating if the original output
#' of [kmeans][stats::kmeans] should be returned in the output (`TRUE` by 
#' default, see Value).
#'
#' @details
#' This method partitions the data into k groups
#' such that that the sum of squares of euclidean distances from points to the
#' assigned cluster centers is minimized. k-means cannot be applied directly
#' on dissimilarity/beta-diversity metrics, because these distances are not
#' euclidean. Therefore, it requires first to transform the dissimilarity
#' matrix with a Principal Coordinate Analysis (using the function
#' [pcoa][ape::pcoa]), and then applying k-means on the coordinates
#' of points in the PCoA. Because this makes an additional transformation of
#' the initial matrix of dissimilarity, the partitioning around medoids method
#' should be preferred ([nhclu_pam])
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
#' [kmeans][stats::kmeans].
#'
#' @author
#' Boris Leroy (\email{leroy.boris@gmail.com}) \cr
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) \cr
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) 
#' 
#' @seealso  [nhclu_pam]
#' @examples
#' comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
#' 20, 25)
#' rownames(comat) <- paste0("Site",1:20)
#' colnames(comat) <- paste0("Species",1:25)
#'
#' comnet <- mat_to_net(comat)
#'
#' dissim <- dissimilarity(comat, metric = "all")
#' 
#' clust1 <- nhclu_kmeans(dissim, n_clust = 2:10, index = "Simpson")
#' clust2 <- nhclu_kmeans(dissim, n_clust = 2:15, index = "Simpson")
#' bioregionalization_metrics(clust2, dissimilarity = dissim,
#'                   eval_metric = "pc_distance")
#' 
#' bioregionalization_metrics(clust2, net = comnet, species_col = "Node2",
#'                   site_col = "Node1", eval_metric = "avg_endemism")
#'
#' @importFrom stats as.dist kmeans
#' @importFrom ape pcoa
#'         
#' @export

nhclu_kmeans <- function(dissimilarity,
                         index = names(dissimilarity)[3],
                         seed = NULL,
                         n_clust = c(1,2,3),
                         iter_max = 10,
                         nstart = 10,
                         algorithm = "Hartigan-Wong",
                         algorithm_in_output = TRUE){
  
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
  controls(args = iter_max, data = NULL, type = "positive_integer")
  controls(args = nstart, data = NULL, type = "positive_integer")
  controls(args = algorithm, data = NULL, type = "character")
  if(!(algorithm %in% c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"))){
    stop("Please choose algorithm among the followings values:
Hartigan-Wong, Lloyd, Forgy or MacQueen", call. = FALSE)
  }
  controls(args = algorithm_in_output, data = NULL, type = "boolean")
  
  # 2. Function ---------------------------------------------------------------
  outputs <- list(name = "nhclu_kmeans")
  
  outputs$args <- list(index = index,
                       seed = seed,
                       n_clust = n_clust,
                       iter_max = iter_max,
                       nstart = nstart,
                       algorithm = algorithm,
                       algorithm_in_output = algorithm_in_output)
  
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
  
  # kmeans only works on Euclidean distances, so the dissimilarity matrix needs
  # to be transformed into a multivariate space with euclidean distances
  # with a Principal Coordinate Analysis
  if(length(unique(dist.obj)) == 1 && unique(dist.obj) == 0){
    stop("All sites are completely dissimilar.")
  } else{
    outputs$clustering_algorithms$pcoa <- ape::pcoa(dist.obj)
  }
  
  # Performing the kmeans on the PCoA with all axes
  if(is.null(seed)){
    outputs$algorithm <- lapply(n_clust,
                                function(x)
                                  stats::kmeans(dist.obj,
                                                centers = x,
                                                iter.max = iter_max,
                                                nstart = nstart,
                                                algorithm = algorithm))
  }else{
    set.seed(seed)
    outputs$algorithm <- lapply(n_clust,
                                function(x)
                                  stats::kmeans(dist.obj,
                                                centers = x,
                                                iter.max = iter_max,
                                                nstart = nstart,
                                                algorithm = algorithm))
    rm(.Random.seed, envir=globalenv())
  }
  names(outputs$algorithm) <- paste0("K_", n_clust)
  
  outputs$clusters <- data.frame(
    outputs$clusters,
    data.frame(lapply(names(outputs$algorithm),
                      function(x) outputs$algorithm[[x]]$cluster)))
  
  outputs$clusters <- knbclu(outputs$clusters, reorder = TRUE)
  
  outputs$cluster_info <- data.frame(
    partition_name = names(outputs$clusters)[2:length(outputs$clusters),
                                             drop = FALSE],
    n_clust = apply(outputs$clusters[, 2:length(outputs$clusters),
                                     drop = FALSE],
                    2, function(x) length(unique(x))))
  
  # Set algorithm in output
  if (!algorithm_in_output) {
    outputs$algorithm <- NA
  }
  
  class(outputs) <-  append("bioregion.clusters", class(outputs))
  
  return(outputs)
}
