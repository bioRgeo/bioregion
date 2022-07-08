#' Non hierarchical clustering: OPTICS
#'
#' This function performs semi-hierarchical
#' clustering on the basis of distances with the OPTICS algorithm (Ordering
#' Points To Identify the Clustering Structure)
#'
#' @param distances the output object from \code{\link{similarity_to_distance}},
#' a \code{data.frame} with the first columns called "Site1" and "Site2", and
#' the other columns being the distance indices or a \code{dist} object
#' @param index a \code{character} string providing the name of the distance
#' index to use, corresponding to the column
#' name in \code{distances}. By default, the third column name of
#'  \code{distances} is used.
#' @param minPts a \code{numeric} value specifying the minPts argument
#' of \link[dbscan:dbscan]{dbscan::dbscan()}). minPts is the minimum number of
#' points to form a dense region.
#' By default, it is set to the
#' natural logarithm of the number of sites in \code{distances}.
#' @param eps a \code{numeric} value specifying the eps argument
#' of \link[dbscan:optics]{dbscan::optics()}). It is the upper limit of the size
#' of the epsilon neighborhood. Limiting the neighborhood size improves
#' performance and has no or very little impact on the ordering as long as it
#' is not set too low. If not specified (default behaviour), the largest
#' minPts-distance in the
#' data set is used which gives the same result as infinity.
#' @param xi a \code{numeric} value specifying the steepness threshold to
#' identify clusters hierarchically using the Xi method
#' (see \link[dbscan:optics]{dbscan::optics()})
#' @param minimum a \code{boolean} specifying if the hierarchy should be pruned
#' out from the output to only keep clusters at the "minimal" level, i.e.
#' only leaf / non-overlapping clusters.
#' If \code{TRUE}, then argument \code{show_hierarchy} should be \code{FALSE}
#' @param show_hierarchy a \code{boolean} specifying if the hierarchy of
#' clusters should be included in the output. By default, the hierarchy is not
#' visible in the clusters obtained from OPTICS - it can only be visualised by
#' visualising the plot of the OPTICS object. If \code{show_hierarchy = TRUE},
#' then the output cluster \code{data.frame} will contain additional columns
#' showing the hierarchy of clusters.
#' @param ... you can add here further arguments to be passed to \code{optics()}
#' (see \link[dbscan:optics]{dbscan::optics()})
#'
#' @details
#' The optics (Ordering points to identify the clustering
#'  structure) is a semi-hierarchical clustering algorithm which orders the
#'  points in the dataset such that points which are closest become neighbours,
#'  and calculates a reachability distance for each point. Then, clusters
#'  can be extracted in a hierarchical manner from this reachability distance,
#'  by identifying clusters depending on changes in the relative cluster
#'  density. The reachability plot should be explored to understand
#'  the clusters and their hierarchical nature, by running plot on the output
#'  of the function: \code{plot(object$algorithm$optics)}.
#'  We recommend reading \insertCite{Hahsler2019}{bioRgeo} to grasp the
#'  algorithm, how it works, and what the clusters mean.
#'
#'  To extract the clusters, we use the
#'  \link[dbscan:extractXi]{dbscan::extractXi()} function which is based on
#'  the steepness of the reachability plot
#'  (see \link[dbscan:optics]{dbscan::optics()})
#'
#'
#' @return
#' A \code{list} of class \code{bioRgeo.clusters} with five slots:
#' \enumerate{
#' \item{\bold{name}: \code{character string} containing the name of the algorihtm}
#' \item{\bold{args}: \code{list} of input arguments as provided by the user}
#' \item{\bold{inputs}: \code{list} of characteristics of the input dataset}
#' \item{\bold{algorithm}: \code{list} of all objects associated with the
#'  clustering procedure, such as original cluster objects}
#' \item{\bold{clusters}: \code{data.frame} containing the clustering results}}
#'
#' @export
#'
#' @examples
#' simil <- similarity(vegemat, metric = "all")
#' distances <- similarity_to_distance(simil)
#'
#' clust1 <- nhclu_optics(distances,
#'     index = "Simpson")
#' clust1
#' head(clust1$clusters)
#'
#' # Visualise the optics plot (the hierarchy of clusters is illustrated at the bottom)
#' plot(clust1$algorithm$optics)
#'
#' # Extract the hierarchy of clusters
#' clust1 <- nhclu_optics(distances,
#'     index = "Simpson",
#'     show_hierarchy = TRUE)
#' clust1
#' head(clust1$clusters)
nhclu_optics <- function(distances,
                         index = names(distances)[3],
                         minPts = NULL,
                         eps = NULL,
                         xi = 0.05,
                         minimum = FALSE,
                         # rename_clusters = TRUE, # to implement?
                         show_hierarchy = FALSE,
                         ...
)
{
  if(inherits(distances, "bioRgeo.pairwise.metric"))
  {
    if(attr(distances, "type") == "similarity")
    {
      stop("distances seems to be a similarity object.
         nhclu_dbscan() should be applied on distances, not similarities.
         Use similarity_to_distance() before using nhclu_dbscan()")
    }
    if(!(index %in% colnames(distances)))
    {
      stop("Argument index should be one of the column names of distance")
    }

  } else if(!any(inherits(distances, "bioRgeo.pairwise.metric"), inherits(distances, "dist")))
  {
    if(!(index %in% colnames(distances)))
    {
      stop("distances is not a bioRgeo.distance object, a distance matrix (class dist) or a data.frame with at least 3 columns (site1, site2, and your distance index)")
    }
  }

  if(!inherits(distances, "dist"))
  {
    dist.obj <- stats::as.dist(
      net_to_mat(distances[, c(1, 2,
                               which(colnames(distances) == index))],
                 weight = TRUE, squared = TRUE, symmetrical = TRUE))

  }

  if(minimum & show_hierarchy)
  {
    warning("When minimum = TRUE, then only the 'minimal' (=leaf/non-overlapping)
    clusters are returned by optics, hence without any hierarchical structure.
    In this case, argument show_hierarchy is not relevant - turning it off.")
    show_hierarchy <- FALSE
  }


  outputs <- list(name = "optics")

  outputs$args <- list(index = index,
                       minPts = minPts,
                       eps = eps,
                       xi = xi,
                       minimum = minimum,
                       show_hierarchy = TRUE,
                       ...
  )

  outputs$inputs <- list(bipartite = FALSE,
                         weight = TRUE,
                         pairwise_metric = TRUE,
                         distance = TRUE,
                         nb_sites = attr(dist.obj, "Size"))

  outputs$algorithm <- list()

  outputs$clusters <- data.frame(matrix(ncol = 1,
                                        nrow = length(labels(dist.obj)),
                                        dimnames = list(labels(dist.obj),
                                                        "name")))

  outputs$clusters$name <- labels(dist.obj)

  if(is.null(minPts))
  {
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

  # The output cluster numbers do not reflect the hierarchy of clusters
  # >> Find a way to extract the hierarchy of clusters?
  # see outputs$algorithm$optics$clusters_xi
  if(!show_hierarchy) {
    outputs$clusters$optics <-
      outputs$algorithm$optics$cluster
  } else {
    cls_hierarchy <- outputs$algorithm$optics$clusters_xi
    cls_hierarchy$diff <- cls_hierarchy$end - cls_hierarchy$start

    cls_order <- cls_hierarchy$cluster_id[order(cls_hierarchy$diff, decreasing = TRUE)]
    cls_hierarchy$new_cls_id <- NA
    for(cls in cls_order)
    {
      cur_start <- cls_hierarchy$start[cls_hierarchy$cluster_id == cls]
      cur_end <- cls_hierarchy$end[cls_hierarchy$cluster_id == cls]

      cur_hier <- cls_hierarchy[- which(cls_hierarchy$cluster_id == cls), ]

      cur_hier$cluster_id[which(cur_start >= cur_hier$start & cur_end <= cur_hier$end)]

      if(cls == cls_order[1]) {
        new.id <- cls
      } else if(any(cur_start >= cur_hier$start & cur_end <= cur_hier$end)) {
        sup_lvl <- cur_hier$new_cls_id[which(cur_start >= cur_hier$start & cur_end <= cur_hier$end)]
        sup_lvl_direct <- sup_lvl[nchar(sup_lvl) == max(nchar(sup_lvl))]
        new.id <- paste0(sup_lvl_direct,
                         ".",
                         cls)
      } else {
        new.id <- cls
      }

      cls_hierarchy$new_cls_id[cls_hierarchy$cluster_id == cls] <- new.id
    }
    max.col <- max(lengths(regmatches(cls_hierarchy$new_cls_id, gregexpr("\\.", cls_hierarchy$new_cls_id)))) + 1
    cls_hierarchy <- tidyr::separate(data = cls_hierarchy,
                                     col = "new_cls_id",
                                     remove = FALSE,
                                     into = paste0("lvl", 1:max.col),
                                     sep = "\\.",
                                     fill = "right")
    outputs$clusters <- data.frame(outputs$clusters,
                                   cls_hierarchy[match(outputs$algorithm$optics$cluster,
                                                       cls_hierarchy$cluster_id),
                                                 c("cluster_id", "new_cls_id", paste0("lvl", 1:max.col))])
  }

  outputs$clusters <- knbclu(outputs$clusters,
                             method = "length",
                             reorder = FALSE)
  class(outputs) <-  append("bioRgeo.clusters", class(outputs))
  return(outputs)
}
