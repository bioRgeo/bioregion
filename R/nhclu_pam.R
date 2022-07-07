#' Non hierarchical clustering: partitioning around medoids
#'
#' This function performs non hierarchical
#' clustering on the basis of distances with partitioning around medoids.
#'
#' @param distances the output object from \code{\link{similarity_to_distance}},
#' a \code{data.frame} with the first columns called "Site1" and "Site2", and
#' the other columns being the distance indices or a \code{dist} object
#' @param index a \code{character} string providing the name of the distance
#' index to use, corresponding to the column
#' name in \code{distances}. By default, the third column name of
#'  \code{distances} is used.
#' @param n_clust an \code{integer} or a \code{vector} of \code{integers}
#' specifying the requested number(s) of clusters
#' @param variant a \code{character} string specifying the variant of pam
#' to use, by default "faster". See \link[cluster:pam]{cluster::pam()} for
#' more details
#' @param nstart an \code{integer} specifying the number of random “starts”
#' for the pam algorithm. By default, 1 (for the \code{"faster"} variant)
#' @param cluster_only a \code{boolean} specifying if only the clustering
#' should be returned from the \link[cluster:pam]{cluster::pam()} function
#' (more efficient)
#' @param ... you can add here further arguments to be passed to \code{pam()}
#' (see \link[cluster:pam]{cluster::pam()})
#'
#' @details
#' This method partitions data into
#'  the chosen number of cluster on the basis of the input distance matrix.
#'  It is more robust than k-means because it minimizes the sum of distances
#'  between cluster centres and points assigned to the cluster -
#'  whereas the k-means approach minimizes the sum of squared euclidean
#'  distances (thus k-means cannot be applied directly on the input distance
#'  matrix if the distances are not euclidean).
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
#' clust1 <- nhclu_pam(distances,
#'     n_clust = 2:10,
#'     index = "Simpson")
#' clust2 <- nhclu_pam(distances,
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
nhclu_pam <- function(distances,
                      index = names(distances)[3],
                      n_clust = NULL,
                      nstart = if(variant == "faster") 1 else NA,
                      variant = "faster", # c("original", "o_1", "o_2", "f_3", "f_4", "f_5", "faster")
                      cluster_only = FALSE,# To reduce computation time & memory, can be provided to cluster functions
                      ... # Further arguments to be passed to cluster::pam
)
{
  if(inherits(distances, "bioRgeo.pairwise.metric"))
  {
    if(attr(distances, "type") == "similarity")
    {
      stop("distances seems to be a similarity object.
         nhclu_pam() should be applied on distances, not similarities.
         Use similarity_to_distance() before using nhclu_pam()")
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



  outputs <- list(name = "pam")

  outputs$args <- list(index = index,
                       n_clust = n_clust,
                       nstart = nstart,
                       variant = variant,
                       cluster_only = cluster_only,
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

  outputs$algorithm$pam <- lapply(n_clust,
                                  function(x)
                                    cluster::pam(dist.obj,
                                                 k = x,
                                                 diss = TRUE,
                                                 keep.diss = FALSE,
                                                 keep.data = FALSE,
                                                 nstart = nstart,
                                                 variant = variant,
                                                 cluster.only = cluster_only,
                                                 ...))

  names(outputs$algorithm$pam) <- paste0("K_", n_clust)

  outputs$clusters <- data.frame(outputs$clusters,
                                 data.frame(lapply(names(outputs$algorithm$pam),
                                                   function(x)
                                                     outputs$algorithm$pam[[x]]$clustering)))
  outputs$clusters <- bioRgeo:::knbclu(outputs$clusters)
  class(outputs) <-  append("bioRgeo.clusters", class(outputs))

  return(outputs)
}
