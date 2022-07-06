#' Non hierarchical clustering: partitioning around medoids
#'
#' This function performs non hierarchical
#' clustering on the basis of distances with partioning around medoids.
#'
#' @param distances the output object from \code{\link{similarity_to_distance}},
#' a \code{data.frame} with the first columns called "Site1" and "Site2", and
#' the other columns being the distance indices or a \code{dist} object
#' @param index a \code{character} string providing the name of the distance
#' index to use, corresponding to the column
#' name in \code{distances}. By default, the third column name of
#'  \code{distances} is used.
#' @param minPts a \code{numeric} value specifying the minPts argument
#' of \link[dbscan:dbscan]{dbscan::dbscan()}). By default, it is set to the
#' natural logarithm of the number of sites in \code{distances}.
#' @param eps a \code{numeric} value specifying the eps argument
#' of \link[dbscan:dbscan]{dbscan::dbscan()}). The value of eps depends on the
#' minPts argument, and should be chosen by identifying a knee in the k-nearest
#' neighbour distance plot. By default the function will try to automatically
#' find a knee, but the result is uncertain, and so the user should inspect the
#' graph and modify \code{dbscan_eps} accordingly.
#' @param plot a \code{boolean} indicating if the  k-nearest
#' neighbour distance plot should be plotted.
#' @param ... you can add here further arguments to be passed to \code{dbscan()}
#' (see \link[dbscan:dbscan]{dbscan::dbscan()})
#'
#' @details
#' The dbscan (Density-based spatial clustering of
#'  applications with noise) clustering algorithm clusters points on the basis
#'  of the density of neighbours around each data points. It necessitates two
#'  main arguments, minPts, which stands for the minimum number of points to
#'  identify a core, and eps, which is the radius to find neighbours.
#'  minPts and eps should be defined by the user, which is not straightforward.
#'  We recommend reading the help in \link[dbscan:dbscan]{dbscan::dbscan()})
#'  to learn how to set these arguments, as well as the paper
#'  \insertCite{Hahsler2019}{bioRgeo}. Note that clusters with a value of 0
#'  are points which were deemed as noise by the algorithm.
#'
#'  By default the function will select values for \code{minPts} and
#'  \code{eps}. However, these values can be inadequate and the users is advised
#'  to tune these values by running the function multiple times.
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
#' clust1 <- nhclu_dbscan(distances,
#'     index = "Simpson")
#' clust2 <- nhclu_dbscan(distances,
#'     n_clust = 2:25,
#'     index = "Simpson")
#' partition_metrics(clust2,
#'                   distances = distances,
#'                   eval_metric = "pc_distance",
#'                   partition_optimisation = TRUE)
#' partition_metrics(clust2,
#'                   sp_site_table = vegemat,
#'                   eval_metric = "avg_endemism",
#'                   partition_optimisation = TRUE)
nhclu_dbscan <- function(distances,
                      index = names(distances)[3],
                      minPts = NULL,
                      eps = NULL,
                      plot = TRUE,
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

  if(!is.null(n_clust)){
    if(is.numeric(n_clust))
    {
      if(any(!(n_clust %% 1 == 0))) # integer testing ain't easy in R
      {
        stop("n_clust must an integer or a vector of integers determining the number of clusters.")
      }
    } else
    {
      stop("n_clust must an integer or a vector of integers determining the number of clusters.")
    }
  }

  if(!inherits(distances, "dist"))
  {
    dist.obj <- stats::as.dist(
      net_to_mat(distances[, c(1, 2,
                               which(colnames(distances) == index))],
                 weight = TRUE, squared = TRUE, symmetrical = TRUE))

  }



  outputs <- list(name = "dbscan")

  outputs$args <- list(index = index,
                       minPts = minPts,
                       eps = eps,
                       plot = plot,
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
                                                        "site")))

  outputs$clusters$site <- labels(dist.obj)

  if(is.null(minPts))
  {
    # Using a default value of minPts if none provided by the user
    minPts <- log(length(labels(dist.obj)))
  }


  knnp <- dbscan::kNNdist(dist.obj,
                          k = minPts - 1)


  x_ <- order(knnp)

  # Trying to find the knee, and not the elbow
  knee <- bioRgeo:::.elbow_finder(x_, max(knnp) - knnp, correct_decrease = TRUE)

  if (is.null(eps)) {
    eps <- knee[2]
  }

  if(plot)
  {
    plot(knnp[order(knnp)], type = "l",
         main = "dbscan parameter choice:\nchoose eps where there is a knee in the curve",
         xlab = "", ylab = "epsilon")
    graphics::abline(h = eps)
    graphics::text(x = 0,
                   y = knee[2] + 0.025 * (max(knnp) - min(knnp)),
                   labels = paste0("Selected eps value: ",
                                   round(eps, 3)),
                   adj = 0)
  }

  outputs$clustering_algorithms$dbscan <-
    dbscan::dbscan(dist.obj,
                   minPts = minPts,
                   eps = eps)

  # NOTE: values of 0 mean "noise / no cluster" with dbscan
  outputs$clusters$dbscan <-
    outputs$clustering_algorithms$dbscan$cluster

  outputs$clusters <- bioRgeo:::knbclu(outputs$clusters)
  class(outputs) <-  append("bioRgeo.clusters", class(outputs))

  return(outputs)
}
