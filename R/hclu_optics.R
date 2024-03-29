#' OPTICS hierarchical clustering algorithm
#'
#' This function performs semi-hierarchical
#' clustering on the basis of dissimilarity with the OPTICS algorithm (Ordering
#' Points To Identify the Clustering Structure)
#'
#' @param dissimilarity the output object from [dissimilarity()] or
#' [similarity_to_dissimilarity()], or a `dist` object. 
#' If a `data.frame` is used, the first two columns represent pairs of
#' sites (or any pair of nodes), and the next column(s) are the dissimilarity
#' indices.
#'  
#' @param index name or number of the dissimilarity column to use. By default, 
#' the third column name of `dissimilarity` is used.
#' 
#' @param minPts a `numeric` value specifying the minPts argument of
#' [dbscan::dbscan()][dbscan::dbscan]). minPts is the minimum number of
#' points to form a dense region. By default, it is set to the natural
#' logarithm of the number of sites in `dissimilarity`.
#' 
#' @param eps a `numeric` value specifying the eps argument of
#' [dbscan::optics()][dbscan::optics]). It is the upper limit of the size
#' of the epsilon neighborhood. Limiting the neighborhood size improves
#' performance and has no or very little impact on the ordering as long as it
#' is not set too low. If not specified (default behavior), the largest
#' minPts-distance in the data set is used which gives the same result as
#' infinity.
#' 
#' @param xi a `numeric` value specifying the steepness threshold to
#' identify clusters hierarchically using the Xi method
#' (see [dbscan::optics()][dbscan::optics])
#' 
#' @param minimum a `boolean` specifying if the hierarchy should be pruned
#' out from the output to only keep clusters at the "minimal" level, i.e.
#' only leaf / non-overlapping clusters.
#' If `TRUE`, then argument `show_hierarchy` should be `FALSE`
#' 
#' @param show_hierarchy a `boolean` specifying if the hierarchy of
#' clusters should be included in the output. By default, the hierarchy is not
#' visible in the clusters obtained from OPTICS - it can only be visualized by
#' visualising the plot of the OPTICS object. If `show_hierarchy = TRUE`,
#' then the output cluster `data.frame` will contain additional columns
#' showing the hierarchy of clusters.
#' 
#' @param ... you can add here further arguments to be passed to `optics()`
#' (see [dbscan::optics()][dbscan::optics])
#'
#' @details
#' The optics (Ordering points to identify the clustering structure) is a
#' semi-hierarchical clustering algorithm which orders the points in the
#' dataset such that points which are closest become neighbors, and calculates
#' a reachability distance for each point. Then, clusters can be extracted in a
#' hierarchical manner from this reachability distance, by identifying clusters
#' depending on changes in the relative cluster density. The reachability plot
#' should be explored to understand the clusters and their hierarchical nature,
#' by running plot on the output of the function:
#' `plot(object$algorithm$optics)`.
#' We recommend reading \insertCite{Hahsler2019}{bioregion} to grasp the
#' algorithm, how it works, and what the clusters mean.
#'
#' To extract the clusters, we use the
#' [dbscan::extractXi()][dbscan::extractXi] function which is based on the
#' steepness of the reachability plot (see
#' [dbscan::optics()][dbscan::optics])
#'
#' @references 
#' \insertRef{Hahsler2019}{bioregion}
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
#' @author
#' Boris Leroy (\email{leroy.boris@gmail.com}),
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) and
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) 
#' 
#' @seealso [nhclu_dbscan] 
#' @examples
#' dissim <- dissimilarity(fishmat, metric = "all")
#'   
#' clust1 <- hclu_optics(dissim, index = "Simpson")
#' clust1
#' 
#' # Visualize the optics plot (the hierarchy of clusters is illustrated at the
#' # bottom)
#' plot(clust1$algorithm$optics)
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
                        # rename_clusters = TRUE, # to implement?
                        show_hierarchy = FALSE,
                        ...){
  if(inherits(dissimilarity, "bioregion.pairwise.metric")){
    if(attr(dissimilarity, "type") == "similarity"){
      stop("dissimilarity seems to be a similarity object.
         nhclu_dbscan() should be applied on dissimilarity, not similarities.
         Use similarity_to_dissimilarity() before using nhclu_dbscan()")
    }
    if(is.numeric(index)){
      index <- names(dissimilarity)[index]
    }
    if(!(index %in% colnames(dissimilarity))){
      stop("Argument index should be one of the column names of dissimilarity")
    }
    
  } else if(!any(inherits(dissimilarity, "bioregion.pairwise.metric"),
                 inherits(dissimilarity, "dist"))){
    if(is.numeric(index)){
      index <- names(dissimilarity)[index]
    }
    if(!(index %in% colnames(dissimilarity))){
      stop("dissimilarity is not a bioregion.pairwise.metric object, a
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
  
  if(minimum & show_hierarchy){
    warning("When minimum = TRUE, then only the 'minimal'
    (=leaf/non-overlapping) clusters are returned by optics, hence without any
    hierarchical structure. In this case, argument show_hierarchy is not
    relevant - turning it off.")
    show_hierarchy <- FALSE
  }
  
  if(!is.null(minPts)){
    if(minPts %% 1 != 0 || minPts < 1){
      stop("minPts must be a positive integer, indicating the number of points
      to form a dense region (see dbscan::dbscan()).")
    }
  }
  
  if(!is.null(eps)){
    if(eps %% 1 != 0 || eps < 1){
      stop("eps must be a positive integer, indicating the upper limit of the
         size of the epsilon neighborhood (see dbscan::optics()).")
    }
  }
  
  if(!is.null(xi)){
    if(xi <= 0 || xi >= 1){
      stop("xi must be a numeric in the ]0, 1[ interval
           (see dbscan::optics()).")
    }
  }
  
  outputs <- list(name = "optics")
  
  outputs$args <- list(index = index,
                       minPts = minPts,
                       eps = eps,
                       xi = xi,
                       minimum = minimum,
                       show_hierarchy = TRUE,
                       ...)
  
  outputs$inputs <- list(bipartite = FALSE,
                         weight = TRUE,
                         pairwise_metric = TRUE,
                         dissimilarity = TRUE,
                         nb_sites = attr(dist.obj, "Size"))
  outputs$inputs$hierarchical <- ifelse(show_hierarchy,
                                        TRUE,
                                        FALSE)
  
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
  
  outputs$algorithm$optics <- dbscan::optics(x = dist.obj,
                                             minPts = minPts,
                                             eps = eps,
                                             ...)
  outputs$algorithm$optics <-
    dbscan::extractXi(outputs$algorithm$optics,
                      xi = xi,
                      minimum = minimum)
  
  if(!show_hierarchy) {
    outputs$clusters$optics <- outputs$algorithm$optics$cluster
  } else {
    cls_hierarchy <- outputs$algorithm$optics$clusters_xi
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
      cls_hierarchy[match(outputs$algorithm$optics$cluster,
                          cls_hierarchy$cluster_id),
                    paste0("lvl", 1:max.col)])
  }
  
  outputs$clusters <- knbclu(outputs$clusters,
                             method = "length",
                             reorder = FALSE)
  
  outputs$cluster_info <- data.frame(
    partition_name = names(outputs$clusters)[2:length(outputs$clusters),
                                             drop = FALSE],
    n_clust = apply(outputs$clusters[, 2:length(outputs$clusters),
                                     drop = FALSE],
                    2, function(x) length(unique(x))))
  
  if(show_hierarchy){
    outputs$cluster_info$hierarchical_level <- 1:max.col
  }
  
  class(outputs) <-  append("bioregion.clusters", class(outputs))
  return(outputs)
}
