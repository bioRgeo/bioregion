#' Non-hierarchical clustering: DBSCAN
#'
#' This function performs non-hierarchical clustering based on dissimilarity 
#' using the Density-Based Spatial Clustering of Applications with Noise 
#' (DBSCAN) algorithm.
#'
#' @param dissimilarity The output object from [dissimilarity()] or 
#' [similarity_to_dissimilarity()], or a `dist` object. If a `data.frame` is 
#' used, the first two columns should represent pairs of sites (or any pair of 
#' nodes), and the subsequent column(s) should contain the dissimilarity indices.
#'
#' @param index The name or number of the dissimilarity column to use. By 
#' default, the third column name of `dissimilarity` is used.
#'
#' @param minPts A `numeric` vector or a single `numeric` value specifying the 
#' `minPts` argument of [dbscan::dbscan()]. `minPts` is the minimum number of 
#' points to form a dense region. By default, it is set to the natural logarithm 
#' of the number of sites in `dissimilarity`. See Details for guidance on 
#' choosing this parameter.
#'
#' @param eps A `numeric` vector or a single `numeric` value specifying the `eps` 
#' argument of [dbscan::dbscan()]. `eps` specifies how similar points should be 
#' to each other to be considered part of a cluster. See Details for guidance on 
#' choosing this parameter.
#'
#' @param plot A `boolean` indicating whether the k-nearest neighbor distance 
#' plot should be displayed.
#'
#' @param algorithm_in_output A `boolean` indicating whether the original output 
#' of [dbscan::dbscan] should be included in the output. Defaults to `TRUE` (see 
#' Value).
#' 
#' @param verbose A `boolean` indicating whether to 
#' display progress messages. Set to `FALSE` to suppress these messages.
#'
#' @param ... Additional arguments to be passed to `dbscan()` (see 
#' [dbscan::dbscan]).
#'
#' @return
#' A `list` of class `bioregion.clusters` with five components:
#' \enumerate{
#' \item{**name**: A `character` string containing the name of the algorithm.}
#' \item{**args**: A `list` of input arguments as provided by the user.}
#' \item{**inputs**: A `list` of characteristics of the clustering process.}
#' \item{**algorithm**: A `list` of all objects associated with the clustering 
#' procedure, such as original cluster objects (only if 
#' `algorithm_in_output = TRUE`).}
#' \item{**clusters**: A `data.frame` containing the clustering results.}}
#'
#' If `algorithm_in_output = TRUE`, the `algorithm` slot includes the output of 
#' [dbscan::dbscan].
#'
#' @details
#' The DBSCAN (Density-Based Spatial Clustering of Applications with Noise) 
#' algorithm clusters points based on the density of neighbors around each 
#' data point. It requires two main arguments: `minPts`, the minimum number of 
#' points to identify a core, and `eps`, the radius used to find neighbors.
#'
#' **Choosing minPts:** This determines how many points are necessary to form a 
#' cluster. For example, what is the minimum number of sites expected in a 
#' bioregion? Choose a value sufficiently large for your dataset and expectations.
#'
#' **Choosing eps:** This determines how similar sites should be to form a 
#' cluster. If `eps` is too small, most points will be considered too distinct 
#' and marked as noise. If `eps` is too large, clusters may merge. The value of 
#' `eps` depends on `minPts`. It is recommended to choose `eps` by identifying 
#' a knee in the k-nearest neighbor distance plot.
#'
#' By default, the function attempts to find a knee in this curve 
#' automatically, but the result is uncertain. Users should inspect the graph 
#' and modify `eps` accordingly. To explore `eps` values, run the function 
#' initially without defining `eps`, review the recommendations, and adjust 
#' as needed based on clustering results.
#'
#' @author
#' Boris Leroy (\email{leroy.boris@gmail.com}) \cr
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) \cr
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) 
#' 
#' @seealso 
#' For more details illustrated with a practical example, 
#' see the vignette: 
#' \url{https://biorgeo.github.io/bioregion/articles/a4_2_non_hierarchical_clustering.html}.
#' 
#' Associated functions: 
#' [nhclu_clara] [nhclu_clarans] [nhclu_kmeans] [nhclu_pam] [nhclu_affprop] 
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
#' @references
#' Hahsler M, Piekenbrock M & Doran D (2019) Dbscan: Fast density-based 
#' clustering with R. \emph{Journal of Statistical Software}, 91(1), 1--30.
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
                         verbose = TRUE,
                         ...){
  
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
    colnameindex <- index
    if(is.numeric(colnameindex)){
      colnameindex <- colnames(net)[index]
      if(is.null(colnameindex)){
        colnameindex <- NA
      }
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
    colnameindex <- NA
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
  controls(args = verbose, data = NULL, type = "boolean")
  
  # 2. Function ----------------------------------------------------------------
  # Output format
  outputs <- list(name = "nhclu_dbscan")
  
  outputs$args <- list(index = index,
                       minPts = minPts,
                       eps = eps,
                       plot = plot,
                       algorithm_in_output = algorithm_in_output,
                       verbose = verbose,
                       ...)
  
  # Determine pairwise_metric and data_type
  pairwise_metric <- ifelse(!inherits(dissimilarity, "dist"), 
                            colnameindex, 
                            NA)
  data_type <- detect_data_type_from_metric(pairwise_metric)
  
  outputs$inputs <- list(bipartite = FALSE,
                         weight = TRUE,
                         pairwise = TRUE,
                         pairwise_metric = pairwise_metric,
                         dissimilarity = TRUE,
                         nb_sites = attr(dist.obj, "Size"),
                         hierarchical = FALSE,
                         data_type = data_type,
                         node_type = "site")
  
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
  
  if (is.null(eps) & verbose) {
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
  outputs$clusters[,-1][outputs$clusters[,-1]==0]=NA
  outputs$clusters <- knbclu(outputs$clusters, reorder = TRUE)
  
  # Add node_type attribute
  attr(outputs$clusters, "node_type") <- rep("site", dim(outputs$clusters)[1])
  
  outputs$cluster_info <- data.frame(
    partition_name = names(outputs$clusters)[2:length(outputs$clusters),
                                             drop = FALSE],
    n_clust = apply(outputs$clusters[, 2:length(outputs$clusters),
                                     drop = FALSE],
                    2, function(x) length(unique(x[!is.na(x)]))),
    cluster_arg_order)
  
  # Set algorithm in output
  if (!algorithm_in_output) {
    outputs$algorithm <- NA
  }
  
  class(outputs) <-  append("bioregion.clusters", class(outputs))
  
  return(outputs)
}
