#' Non hierarchical clustering: CLARA
#'
#' This function performs non hierarchical clustering on the basis of
#' dissimilarity with partitioning around medoids, using the Clustering Large
#' Applications (CLARA) algorithm.
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
#' @param n_clust an `integer` or an `integer` vector specifying the
#' requested number(s) of clusters.
#' 
#' @param maxiter an `integer` defining the maximum number of iterations.
#' 
#' @param initializer a `character`, either 'BUILD' (used in classic PAM
#' algorithm) or 'LAB' (linear approximative BUILD).
#' 
#' @param fasttol positive `numeric` defining the tolerance for fast swapping
#' behavior, set to 1 by default.
#' 
#' @param numsamples positive `integer` defining the number of samples to draw.
#' 
#' @param sampling positive `numeric` defining the sampling rate.
#' 
#' @param independent a `boolean` indicating that the previous
#' medoids are not kept in the next sample (FALSE by default).
#' 
#' @param algorithm_in_output a `boolean` indicating if the original output
#' of [fastclara][fastkmedoids::fastclara] should be returned in the output 
#' (`TRUE` by default, see Value).
#' 
#' @details
#' Based on [fastkmedoids](https://cran.r-project.org/package=fastkmedoids)
#' package ([fastclara][fastkmedoids::fastclara]).
#'
#' @return
#' A `list` of class `bioregion.clusters` with five slots:
#' \enumerate{
#' \item{**name**: `character` containing the name of the algorithm}
#' \item{**args**: `list` of input arguments as provided by the user}
#' \item{**inputs**: `list` of characteristics of the clustering process}
#' \item{**algorithm**: `list` of all objects associated with the
#'  clustering procedure, such as original cluster objects (only if
#'  `algorithm_in_output = TRUE`)}
#' \item{**clusters**: `data.frame` containing the clustering results}}
#' 
#' In the `algorithm` slot, if `algorithm_in_output = TRUE`, users can
#' find the output of
#' [fastclara][fastkmedoids::fastclara].
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
#' clust1 <- nhclu_clara(dissim, index = "Simpson", n_clust = 5)
#' 
#' partition_metrics(clust1, dissimilarity = dissim,
#' eval_metric = "pc_distance")
#' 
#'    
#' @importFrom stats as.dist
#' @importFrom fastkmedoids fastclara    
#'                    
#' @export

nhclu_clara <- function(dissimilarity,
                        index = names(dissimilarity)[3],
                        seed = NULL,
                        n_clust = c(1,2,3),
                        maxiter = 0,
                        initializer = "LAB",
                        fasttol = 1,
                        numsamples = 5,
                        sampling = 0.25,
                        independent = FALSE,
                        algorithm_in_output = TRUE){
  
  # 1. Controls ----------------------------------------------------------------
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
  controls(args = maxiter, data = NULL, type = "positive_integer")
  controls(args = initializer, data = NULL, type = "character")
  if(!(initializer %in% c("BUILD", "LAB"))){
    stop("Please choose initializer among the followings values:
BUILD or LAB", call. = FALSE)
  }
  controls(args = fasttol, data = NULL, type = "positive_numeric")
  controls(args = numsamples, data = NULL, type = "positive_integer")
  controls(args = sampling, data = NULL, type = "positive_numeric")
  controls(args = independent, data = NULL, type = "boolean")
  controls(args = algorithm_in_output, data = NULL, type = "boolean")
  
  # 2. Function ---------------------------------------------------------------
  # Output format
  outputs <- list(name = "nhclu_clara")
  
  outputs$args <- list(index = index,
                       seed = seed,
                       n_clust = n_clust,
                       maxiter = maxiter,
                       initializer = initializer,
                       fasttol = fasttol,
                       numsamples = numsamples,
                       sampling = sampling,
                       independent = independent,
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
  
  # CLARA algorithm
  if(is.null(seed)){
    outputs$algorithm <-
      lapply(n_clust,
             function(x)
               fastkmedoids::fastclara(rdist = dist.obj,
                                       n = nrow(dist.obj),
                                       k = x,
                                       maxiter = maxiter,
                                       initializer = initializer,
                                       fasttol = fasttol,
                                       numsamples = numsamples,
                                       sampling = sampling,
                                       independent = independent,
                                       seed = as.numeric(as.POSIXct(
                                         Sys.time())) + sample(-10:10, 1)))
  }else{
    outputs$algorithm <-
      lapply(n_clust,
             function(x)
               fastkmedoids::fastclara(rdist = dist.obj,
                                       n = nrow(dist.obj),
                                       k = x,
                                       maxiter = maxiter,
                                       initializer = initializer,
                                       fasttol = fasttol,
                                       numsamples = numsamples,
                                       sampling = sampling,
                                       independent = independent,
                                       seed = seed))
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
