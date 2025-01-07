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
#' @param seed for the random number generator (NULL for random by default).
#' 
#' @param n_clust an `integer` vector or a single `integer` value specifying the
#' requested number(s) of clusters.
#' 
#' @param numlocal an `integer` defining the number of samples to draw.
#' 
#' @param maxneighbor a positive `numeric` defining the sampling rate.
#' 
#' @param algorithm_in_output a `boolean` indicating if the original output
#' of [fastclarans][fastkmedoids::fastclarans] should be returned in the output
#'  (`TRUE` by default, see Value).
#' 
#' @details
#' Based on [fastkmedoids](https://cran.r-project.org/package=fastkmedoids)
#' package ([fastclarans][fastkmedoids::fastclarans]).
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
#' [fastclarans][fastkmedoids::fastclarans].
#' 
#' @references 
#' Schubert E & Rousseeuw PJ (2019) Faster k-Medoids Clustering: Improving the 
#' PAM, CLARA, and CLARANS Algorithms. \emph{Similarity Search and Applications}
#' , 11807, 171-187.
#' 
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) \cr
#' Boris Leroy (\email{leroy.boris@gmail.com}) \cr
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
#' bioregionalization_metrics(clust1, dissimilarity = dissim,
#' eval_metric = "pc_distance")
#' 
#'    
#' @importFrom stats as.dist
#' @importFrom fastkmedoids fastclarans    
#'                    
#' @export

nhclu_clarans <- function(dissimilarity,
                          index = names(dissimilarity)[3],
                          seed = NULL,
                          n_clust = c(1,2,3),
                          numlocal = 2,
                          maxneighbor = 0.025,
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
  controls(args = numlocal, data = NULL, type = "positive_integer")
  controls(args = maxneighbor, data = NULL, type = "positive_numeric")
  controls(args = algorithm_in_output, data = NULL, type = "boolean")
  
  # 2. Function ---------------------------------------------------------------
  # Output format
  outputs <- list(name = "nhclu_clarans")
  
  outputs$args <- list(index = index,
                       seed = seed,
                       n_clust = n_clust,
                       numlocal = numlocal,
                       maxneighbor = maxneighbor,
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
  
  # CLARANS algorithm
  if(!is.null(seed)){
    outputs$algorithm <-
      lapply(n_clust,
             function(x)
               fastkmedoids::fastclarans(rdist = dist.obj,
                                         n = nrow(dist.obj),
                                         k = x,
                                         numlocal = numlocal,
                                         maxneighbor = maxneighbor,
                                         seed = seed))
  }else{
    outputs$algorithm <-
      lapply(n_clust,
             function(x)
               fastkmedoids::fastclarans(rdist = dist.obj,
                                         n = nrow(dist.obj),
                                         k = x,
                                         numlocal = numlocal,
                                         maxneighbor = maxneighbor,
                                         seed = seedrng()))
  }  
  names(outputs$algorithm) <- paste0("K_", n_clust)
  
  outputs$clusters <- data.frame(
    outputs$clusters,
    data.frame(lapply(names(outputs$algorithm),
                      function(x) outputs$algorithm[[x]]@assignment)))
  
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
  
  class(outputs) <- append("bioregion.clusters", class(outputs))
  
  return(outputs)
}
