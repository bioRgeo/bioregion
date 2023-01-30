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
#' @param minPts a `numeric` value or a vector of `numeric` values
#' specifying the minPts argument of [dbscan::dbscan()][dbscan::dbscan]).
#' minPts is the minimum number of points to form a dense region. By default,
#' it is set to the natural logarithm of the number of sites in
#' `dissimilarity`. See details for guidance on choosing this parameter.
#' 
#' @param eps a `numeric` value or a vector of `numeric` values specifying the
#' eps argument of [dbscan::dbscan()][dbscan::dbscan]). eps specifies how
#' similar points should be to each other to be considered a part of a cluster.
#' See details for guidance on choosing this parameter.
#' 
#' @param plot a `boolean` indicating if the  k-nearest neighbor distance plot
#' should be plotted.
#' 
#' @param ... you can add here further arguments to be passed to `dbscan()`
#' (see [dbscan::dbscan()][dbscan::dbscan])
#'
#' @details
#' The dbscan (Density-based spatial clustering of
#' applications with noise) clustering algorithm clusters points on the basis
#' of the density of neighbours around each data points. It necessitates two
#' main arguments, minPts, which stands for the minimum number of points to
#' identify a core, and eps, which is the radius to find neighbors.
#' minPts and eps should be defined by the user, which is not straightforward.
#' We recommend reading the help in [dbscan::dbscan()][dbscan::dbscan])
#' to learn how to set these arguments, as well as the paper
#' \insertCite{Hahsler2019}{bioRgeo}. Note that clusters with a value of 0
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
#' **Choosing eps:** how similar should sites be in a cluster?  If eps is
#' too small, then a majority of points will be considered to distinct and will
#' not be clustered at all (i.e., considered as noise)? If the value is too
#' high, then clusters will merge together.
#' The value of eps depends on the minPts argument, and the literature
#' recommends to choose eps by identifying a knee in the k-nearest neighbor
#' distance plot. By default
#' the function will try to automatically find a knee in that curve, but the
#' result is uncertain, and so the user should inspect the graph and modify
#' dbscan_eps accordingly.To explore eps values, follow the
#' recommendation by the function when you launch it a first time without
#' defining eps. Then, adjust depending on your clustering results.
#'
#' @return
#' A `list` of class `bioRgeo.clusters` with five slots:
#' \enumerate{
#' \item{**name**: `character string` containing the name of the algorithm}
#' \item{**args**: `list` of input arguments as provided by the user}
#' \item{**inputs**: `list` of characteristics of the input dataset}
#' \item{**algorithm**: `list` of all objects associated with the
#'  clustering procedure, such as original cluster objects}
#' \item{**clusters**: `data.frame` containing the clustering results}}
#'
#' @author
#' Boris Leroy (\email{leroy.boris@gmail.com}),
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) and
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) 
#' 
#' @seealso [hclu_optics] 
#'
#' @examples
#' \dontrun{
#' dissim <- dissimilarity(fishmat, metric = "all")
#' 
#' clust1 <- nhclu_dbscan(dissim, index = "Simpson")
#' clust2 <- nhclu_dbscan(dissim, index = "Simpson", eps = 0.2)
#' clust3 <- nhclu_dbscan(dissim, index = "Simpson", minPts = c(5, 10, 15, 20),
#'      eps = c(.1, .15, .2, .25, .3))
#' }
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
                         ...){
  
  if(inherits(dissimilarity, "bioRgeo.pairwise.metric")){
    if(attr(dissimilarity, "type") == "similarity"){
      stop("dissimilarity seems to be a similarity object.
         nhclu_dbscan() should be applied on dissimilarity, not similarities.
         Use similarity_to_dissimilarity() before using nhclu_dbscan()")
    }
    if(is.numeric(index)){
      index <- names(dissimilarity)[index]
    }
    if(!(index %in% colnames(dissimilarity))) {
      stop("Argument index should be one of the column names of dissimilarity")
    }
  } else if(!any(inherits(dissimilarity, "bioRgeo.pairwise.metric"),
                 inherits(dissimilarity, "dist"))) {
    if(is.numeric(index)){
      index <- names(dissimilarity)[index]
    }
    if(is.null(index) || !(index %in% colnames(dissimilarity))){
      stop("dissimilarity is not a bioRgeo.pairwise.metric object, a
           dissimilarity matrix (class dist) or a data.frame with at least 3
           columns (site1, site2, and your dissimilarity index)")
    }
  }
  
  if(!inherits(dissimilarity, "dist")){
    dist.obj <- stats::as.dist(
      net_to_mat(dissimilarity[, c(colnames(dissimilarity)[1:2], index)],
                 weight = TRUE, squared = TRUE, symmetrical = TRUE))
    
  } else {
    dist.obj <- dissimilarity
  }
  
  if(!is.null(minPts) && !is.numeric(minPts)){
    stop("minPts must be a numeric.")
  }
  
  if(!is.null(eps) && !is.numeric(eps)){
    stop("eps must be a numeric.")
  }
  
  if(!is.logical(plot)){
    stop("plot must be a Boolean.")
  }
  
  outputs <- list(name = "dbscan")
  
  outputs$args <- list(index = index,
                       minPts = minPts,
                       eps = eps,
                       plot = plot,
                       ...)
  
  outputs$inputs <- list(bipartite = FALSE,
                         weight = TRUE,
                         pairwise_metric = TRUE,
                         dissimilarity = TRUE,
                         nb_sites = attr(dist.obj, "Size"))
  
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
      outputs$algorithm$dbscan[[paste0("dbscan_minPts",
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
    data.frame(lapply(names(outputs$algorithm$dbscan),
                      function(x) outputs$algorithm$dbscan[[x]]$cluster)))
  
  outputs$clusters <- knbclu(outputs$clusters, reorder = FALSE)
  
  outputs$cluster_info <- data.frame(
    partition_name = names(outputs$clusters)[2:length(outputs$clusters),
                                             drop = FALSE],
    n_clust = apply(outputs$clusters[, 2:length(outputs$clusters),
                                     drop = FALSE],
                    2, function(x) length(unique(x))),
    cluster_arg_order)
  
  class(outputs) <-  append("bioRgeo.clusters", class(outputs))
  
  return(outputs)
}
