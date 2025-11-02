#' OPTICS hierarchical clustering algorithm
#'
#' This function performs semi-hierarchical clustering based on dissimilarity 
#' using the OPTICS algorithm (Ordering Points To Identify the 
#' Clustering Structure).
#'
#' @param dissimilarity The output object from [dissimilarity()] or
#' [similarity_to_dissimilarity()], or a `dist` object. 
#' If a `data.frame` is used, the first two columns represent pairs of
#' sites (or any pair of nodes), and the subsequent column(s) contain the 
#' dissimilarity indices.
#'  
#' @param index The name or number of the dissimilarity column to use. By 
#' default, the third column name of `dissimilarity` is used.
#' 
#' @param minPts A `numeric` value specifying the minPts argument of
#' [dbscan][dbscan::dbscan]. minPts is the minimum number of points required
#' to form a dense region. By default, it is set to the natural logarithm 
#' of the number of sites in `dissimilarity`.
#' 
#' @param eps A `numeric` value specifying the eps argument of
#' [optics][dbscan::optics]. It defines the upper limit of the size
#' of the epsilon neighborhood. Limiting the neighborhood size improves
#' performance and has no or very little impact on the ordering as long as it
#' is not set too low. If not specified (default behavior), the largest
#' minPts-distance in the dataset is used, which gives the same result as
#' infinity.
#' 
#' @param xi A `numeric` value specifying the steepness threshold to
#' identify clusters hierarchically using the Xi method
#' (see [optics][dbscan::optics]).
#' 
#' @param minimum A `boolean` specifying whether the hierarchy should be pruned
#' from the output to only retain clusters at the "minimal" level, i.e.,
#' only leaf / non-overlapping clusters.
#' If `TRUE`, then the argument `show_hierarchy` should be set to `FALSE`.
#' 
#' @param show_hierarchy A `boolean` specifying whether the hierarchy of
#' clusters should be included in the output. By default, the hierarchy is not
#' visible in the clusters obtained from OPTICS; it can only be visualized by
#' plotting the OPTICS object. If `show_hierarchy = TRUE`,
#' the output cluster `data.frame` will contain additional columns
#' showing the hierarchy of clusters.
#' 
#' @param algorithm_in_output A `boolean` indicating whether the original output
#' of [dbscan][dbscan::dbscan] should be returned in the output (`TRUE` by 
#' default, see Value).
#' 
#' @param ... Additional arguments to be passed to `optics()`
#' (see [optics][dbscan::optics]).
#' 
#' @return
#' A `list` of class `bioregion.clusters` with five slots:
#' \enumerate{
#' \item{**name**: A `character` string containing the name of the algorithm.}
#' \item{**args**: A `list` of input arguments as provided by the user.}
#' \item{**inputs**: A `list` describing the characteristics of the clustering process.}
#' \item{**algorithm**: A `list` containing all objects associated with the
#'  clustering procedure, such as the original cluster objects.}
#' \item{**clusters**: A `data.frame` containing the clustering results.}}
#'
#' In the `algorithm` slot, if `algorithm_in_output = TRUE`, users can
#' find the output of [optics][dbscan::optics].
#'
#' @details
#' The OPTICS (Ordering points to identify the clustering structure) is a
#' semi-hierarchical clustering algorithm which orders the points in the
#' dataset such that points which are closest become neighbors, and calculates
#' a reachability distance for each point. Then, clusters can be extracted in a
#' hierarchical manner from this reachability distance, by identifying clusters
#' depending on changes in the relative cluster density. The reachability plot
#' should be explored to understand the clusters and their hierarchical nature,
#' by running plot on the output of the function 
#' if `algorithm_in_output = TRUE`: `plot(object$algorithm)`.
#' We recommend reading (Hahsler et al., 2019) to grasp the
#' algorithm, how it works, and what the clusters mean.
#'
#' To extract the clusters, we use the
#' [extractXi][dbscan::extractXi] function which is based on the
#' steepness of the reachability plot (see
#' [optics][dbscan::optics])
#'
#' @references 
#' Hahsler M, Piekenbrock M & Doran D (2019) Dbscan: Fast density-based 
#' clustering with R. \emph{Journal of Statistical Software} 91, 1--30. 
#' 
#' @seealso 
#' For more details illustrated with a practical example, 
#' see the vignette: 
#' \url{https://biorgeo.github.io/bioregion/articles/a4_1_hierarchical_clustering.html}.
#' 
#' Associated functions: 
#' [nhclu_dbscan] 
#' 
#' @author
#' Boris Leroy (\email{leroy.boris@gmail.com}) \cr
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) \cr
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) 
#' 
#' @examples
#' dissim <- dissimilarity(fishmat, metric = "all")
#'   
#' clust1 <- hclu_optics(dissim, index = "Simpson")
#' clust1
#' 
#' # Visualize the optics plot (the hierarchy of clusters is illustrated at the
#' # bottom)
#' plot(clust1$algorithm)
#'
#' # Extract the hierarchy of clusters
#' clust1 <- hclu_optics(dissim, index = "Simpson", show_hierarchy = TRUE)
#' clust1
#' 
#' @importFrom stats as.dist
#' @importFrom dbscan optics extractXi
#' @importFrom tidyr separate
#' 
#' @export

hclu_optics <- function(dissimilarity,
                        index = names(dissimilarity)[3],
                        minPts = NULL,
                        eps = NULL,
                        xi = 0.05,
                        minimum = FALSE,
                        show_hierarchy = FALSE,
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
    controls(args = minPts, data = NULL, type = "strict_positive_integer")
  }
  if(!is.null(eps)){
    controls(args = eps, data = NULL, type = "strict_positive_integer")
  }
  controls(args = xi, data = NULL, type = "strict_positive_numeric")
  if (xi >= 1) {
    stop("xi must be in the interval (0,1), (see dbscan::optics())", 
         call. = FALSE)
  }
  controls(args = minimum, data = NULL, type = "boolean")
  controls(args = show_hierarchy, data = NULL, type = "boolean")
  controls(args = algorithm_in_output, data = NULL, type = "boolean")
  
  if(minimum & show_hierarchy){
    #warning("When minimum = TRUE, then only the 'minimal'
    #(=leaf/non-overlapping) clusters are returned by optics, hence without any
    #hierarchical structure. In this case, argument show_hierarchy is not
    #relevant - turning it off.")
    show_hierarchy <- FALSE
  }
  
  # 2. Function ----------------------------------------------------------------
  # Output format
  outputs <- list(name = "hclu_optics")
  
  outputs$args <- list(index = index,
                       minPts = minPts,
                       eps = eps,
                       xi = xi,
                       minimum = minimum,
                       show_hierarchy = show_hierarchy,
                       algorithm_in_output = algorithm_in_output,
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
                         hierarchical = show_hierarchy,
                         data_type = data_type,
                         node_type = "site")
  
  outputs$algorithm <- list()
  
  outputs$clusters <- data.frame(matrix(ncol = 1,
                                        nrow = length(labels(dist.obj)),
                                        dimnames = list(labels(dist.obj),
                                                        "name")))
  
  outputs$clusters$name <- labels(dist.obj)
  
  if(is.null(minPts)) {
    # Using a default value of minPts if none provided by the user
    minPts <- log(length(labels(dist.obj)))
  }
  
  outputs$algorithm <- dbscan::optics(x = dist.obj,
                                             minPts = minPts,
                                             eps = eps,
                                             ...)
  outputs$algorithm <-
    dbscan::extractXi(outputs$algorithm,
                      xi = xi,
                      minimum = minimum)
  
  if(!show_hierarchy) {
    outputs$clusters$optics <- outputs$algorithm$cluster
  } else {
    cls_hierarchy <- outputs$algorithm$clusters_xi
    cls_hierarchy$diff <- cls_hierarchy$end - cls_hierarchy$start
    
    cls_order <- cls_hierarchy$cluster_id[order(cls_hierarchy$diff,
                                                decreasing = TRUE)]
    cls_hierarchy$new_cls_id <- NA
    for(cls in cls_order){
      cur_start <- cls_hierarchy$start[cls_hierarchy$cluster_id == cls]
      cur_end <- cls_hierarchy$end[cls_hierarchy$cluster_id == cls]
      
      cur_hier <- cls_hierarchy[- which(cls_hierarchy$cluster_id == cls), ]
      
      cur_hier$cluster_id[which(cur_start >= cur_hier$start &
                                  cur_end <= cur_hier$end)]
      
      if(cls == cls_order[1]) {
        new.id <- cls
      } else if(any(cur_start >= cur_hier$start & cur_end <= cur_hier$end)) {
        sup_lvl <- cur_hier$new_cls_id[which(cur_start >= cur_hier$start &
                                               cur_end <= cur_hier$end)]
        sup_lvl_direct <- sup_lvl[nchar(sup_lvl) == max(nchar(sup_lvl))]
        new.id <- paste0(sup_lvl_direct, ".", cls)
      } else {
        new.id <- cls
      }
      
      cls_hierarchy$new_cls_id[cls_hierarchy$cluster_id == cls] <- new.id
    }
    
    max.col <- max(lengths(
      regmatches(cls_hierarchy$new_cls_id,
                 gregexpr("\\.", cls_hierarchy$new_cls_id)))) + 1
    if("xics" %in% class(cls_hierarchy)) {
      class(cls_hierarchy) <- class(cls_hierarchy)[
        -which(class(cls_hierarchy) == "xics")]
    }
    cls_hierarchy <- 
      as.data.frame(tidyr::separate_wider_delim(data = cls_hierarchy,
                                                cols = "new_cls_id",
                                                delim = ".",
                                                cols_remove = FALSE,
                                                names = paste0("lvl",
                                                               1:max.col),
                                                too_few = "align_start"
      ))
    
    cls_hierarchy[which(is.na(cls_hierarchy), arr.ind = TRUE)] <- 0
    
    for(lvl in grep("lvl", colnames(cls_hierarchy))[2:max.col]){
      cls_hierarchy[, lvl] <- paste(cls_hierarchy[, lvl - 1],
                                    cls_hierarchy[, lvl],
                                    sep = ".")
    }
    
    cls_hierarchy[grep("lvl", colnames(cls_hierarchy))] <-
      lapply(cls_hierarchy[grep("lvl", colnames(cls_hierarchy))],
             function(x) gsub("\\.0", "", x))
    
    outputs$clusters <- data.frame(
      outputs$clusters,
      cls_hierarchy[match(outputs$algorithm$cluster,
                          cls_hierarchy$cluster_id),
                    paste0("lvl", 1:max.col)])
    
  }
  
  outputs$clusters[,-1][outputs$clusters[,-1]==0]=NA
  outputs$clusters <- knbclu(outputs$clusters, reorder = FALSE)
  
  # Add node_type attribute
  attr(outputs$clusters, "node_type") <- rep("site", dim(outputs$clusters)[1])
  
  outputs$cluster_info <- data.frame(
    partition_name = names(outputs$clusters)[2:length(outputs$clusters),
                                             drop = FALSE],
    n_clust = apply(outputs$clusters[, 2:length(outputs$clusters),
                                     drop = FALSE],
                    2, function(x) length(unique(x[!is.na(x)]))))
  
  
  if(show_hierarchy){
    outputs$cluster_info$hierarchical_level <- 1:max.col
  }
  
  # Set algorithm in output
  if (!algorithm_in_output) {
    outputs$algorithm <- NA
  }
  
  class(outputs) <-  append("bioregion.clusters", class(outputs))
  return(outputs)
}
