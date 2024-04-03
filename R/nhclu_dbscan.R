#' dbscan clustering
#'
#' This function performs non hierarchical clustering on the basis of
#' dissimilarity with Density-based Spatial Clustering of Applications with
#' Noise (DBSCAN)
#'
#' @param dissimilarity the output object from [dissimilarity()] or
#' [similarity_to_dissimilarity()], or a `dist` object. If a `data.frame` is
#' used, the first two columns represent pairs of sites (or any pair of nodes),
#' and the next column(s) are the dissimilarity indices.
#' 
#' @param index name or number of the dissimilarity column to use. By default, 
#' the third column name of `dissimilarity` is used.
#' 
#' @param minPts a `numeric` value or a `numeric` vector
#' specifying the minPts argument of [dbscan::dbscan()][dbscan::dbscan]).
#' minPts is the minimum number of points to form a dense region. By default,
#' it is set to the natural logarithm of the number of sites in
#' `dissimilarity`. See details for guidance on choosing this parameter.
#' 
#' @param eps a `numeric` value or a `numeric` vector specifying the
#' eps argument of [dbscan][dbscan::dbscan]). eps specifies how
#' similar points should be to each other to be considered a part of a cluster.
#' See details for guidance on choosing this parameter.
#' 
#' @param plot a `boolean` indicating if the  k-nearest neighbor distance plot
#' should be plotted.
#' 
#' @param algorithm_in_output a `boolean` indicating if the original output
#' of [dbscan][dbscan::dbscan] should be returned in the output (`TRUE` by 
#' default, see Value).
#' 
#' @param ... you can add here further arguments to be passed to `dbscan()`
#' (see [dbscan][dbscan::dbscan]).
#' 
#' @details
#' The dbscan (Density-based spatial clustering of
#' applications with noise) clustering algorithm clusters points on the basis
#' of the density of neighbours around each data points. It necessitates two
#' main arguments, `minPts`, which stands for the minimum number of points to
#' identify a core, and `eps`, which is the radius to find neighbors.
#' `minPts` and `eps` should be defined by the user, which is not
#' straightforward.
#' We recommend reading the help in [dbscan][dbscan::dbscan])
#' to learn how to set these arguments, as well as the paper
#' \insertCite{Hahsler2019}{bioregion}. Note that clusters with a value of 0
#' are points which were deemed as noise by the algorithm.
#'
#' By default the function will select values for `minPts` and `eps`. However,
#' these values can be inadequate and the users is advised to tune these values
#' by running the function multiple times.
#'
#' **Choosing minPts:** how many points should be necessary to make a cluster?
#' i.e., what is the minimum number of sites you expect in a bioregion? Set a
#' value sufficiently large for your dataset and your expectations.
#'
#' **Choosing eps:** how similar should sites be in a cluster?  If `eps` is
#' too small, then a majority of points will be considered too distinct and
#' will not be clustered at all (i.e., considered as noise)? If the value is
#' too high, then clusters will merge together.
#' The value of `eps` depends on the `minPts` argument, and the literature
#' recommends to choose `eps` by identifying a knee in the k-nearest neighbor
#' distance plot. By default
#' the function will try to automatically find a knee in that curve, but the
#' result is uncertain, and so the user should inspect the graph and modify
#' `dbscan_eps` accordingly.To explore eps values, follow the
#' recommendation by the function when you launch it a first time without
#' defining `eps`. Then, adjust depending on your clustering results.
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
#' [dbscan][dbscan::dbscan].
#'
#' @author
#' Boris Leroy (\email{leroy.boris@gmail.com}),
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) and
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) 
#' 
#' @seealso [hclu_optics] 
#'
#' @examples
#' comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
#' 20, 25)
#' rownames(comat) <- paste0("Site",1:20)
#' colnames(comat) <- paste0("Species",1:25)
#' 
#' dissim <- dissimilarity(comat, metric = "all")
#' 
#' clust1 <- nhclu_dbscan(dissim, index = "Simpson")
#' clust2 <- nhclu_dbscan(dissim, index = "Simpson", eps = 0.2)
#' clust3 <- nhclu_dbscan(dissim, index = "Simpson", minPts = c(5, 10, 15, 20),
#'      eps = c(.1, .15, .2, .25, .3))
#'
#' @importFrom stats as.dist
#' @importFrom dbscan kNNdist dbscan
#' @importFrom graphics abline text
#'
#' @export

nhclu_dbscan <- function(dissimilarity,
                         index = names(dissimilarity)[3],
                         minPts = NULL,
                         eps = NULL,
                         plot = TRUE,
                         algorithm_in_output = TRUE,
                         ...){
  
  # 1. Controls ----------------------------------------------------------------
  controls(args = NULL, data = dissimilarity, type = "input_nhandhclu")
  if(!inherits(dissimilarity, "dist")){
    controls(args = NULL, data = dissimilarity, type = "input_dissimilarity")
    controls(args = NULL, data = dissimilarity, 
             type = "input_data_frame_nhandhclu")
    controls(args = index, data = dissimilarity, type = "input_net_index")
    net <- dissimilarity
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
  
  
  if(!is.null(minPts)){
    controls(args = minPts, data = NULL, 
             type = "integer_vector")
    if(sum(minPts <= 2) > 0){
      stop("minPts must be strictly higher than 2.", call. = FALSE)
    }
  }
  if(!is.null(eps)){
    controls(args = eps, data = NULL, type = "positive_numeric_vector")
  }
  controls(args = plot, data = NULL, type = "boolean")
  controls(args = algorithm_in_output, data = NULL, type = "boolean")
  
  # 2. Function ----------------------------------------------------------------
  # Output format
  outputs <- list(name = "nhclu_dbscan")
  
  outputs$args <- list(index = index,
                       minPts = minPts,
                       eps = eps,
                       plot = plot,
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
  
  if(is.null(minPts)){
    # Using a default value of minPts if none provided by the user
    minPts <- log(length(labels(dist.obj)))
  }
  
  if (is.null(eps)) {
    message(
      "Trying to find a knee in the curve to search for an optimal eps value...
       NOTE: this automatic identification of the knee may not work properly
       if the curve has knees and elbows. Please adjust eps manually by
       inspecting the curve, identifying a knee as follows:

                           /
                 curve    /
              ___________/  <- knee
  elbow ->   /
            /
           /")
  }
  
  # DBSCAN algorithm
  cluster_arg_order <- data.frame()
  
  for(minPtsi in minPts){
    knnp <- dbscan::kNNdist(dist.obj, k = minPtsi - 1)
    x_ <- order(knnp)
    
    # Trying to find the knee, and not the elbow
    if (is.null(eps)) {
      knee <- .elbow_finder(x_, max(knnp) - knnp, correct_decrease = TRUE)
      
      eps <- knee[2]
    }
    
    if(plot){
      plot(
        knnp[order(knnp)], type = "l",
        main = "dbscan parameter choice:\nchoose eps where there is a knee in
        the curve",
        xlab = "", ylab = "epsilon")
      graphics::abline(h = eps)
      graphics::text(x = 0,
                     y = eps + 0.025 * (max(knnp) - min(knnp)),
                     labels = paste0("Selected eps value: ",
                                     round(eps, 3)),
                     adj = 0)
    }
    
    for(epsi in eps) {
      outputs$algorithm[[paste0("dbscan_minPts",
                                       minPtsi,
                                       "_eps", epsi)]]<-
        dbscan::dbscan(dist.obj,
                       minPts = minPtsi,
                       eps = epsi)
      cluster_arg_order <- rbind(cluster_arg_order, 
                                 data.frame(minPts = minPtsi,
                                            eps = epsi))
    }
  }
  
  # NOTE: values of 0 mean "noise / no cluster" with dbscan
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
                    2, function(x) length(unique(x))),
    cluster_arg_order)
  
  # Set algorithm in output
  if (!algorithm_in_output) {
    outputs$algorithm <- NA
  }
  
  class(outputs) <-  append("bioregion.clusters", class(outputs))
  
  return(outputs)
}
