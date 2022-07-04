#' Non hierarchical clustering based on based on distances or beta-diversity
#'
#' This function includes several algorithms to perform non hierarchical
#' clustering on the basis of distances: partioning around medoids / k-medoids,
#' k-means, dbscan, optics, and model-based clustering based on Gaussian
#'  mixture models.
#'
#' @param distances the output object from \code{\link{similarity_to_distance}},
#' a \code{data.frame} with the first columns called "Site1" and "Site2", and
#' the other columns being the distance indices or a \code{dist} object
#' @param index a \code{character} string providing the name of the distance
#' index to use, corresponding to the column
#' name in \code{distances}. By default, the third column name of
#'  \code{distances} is used.
#' @param method a \code{character} vector specifying the methods to use. Available
#' methods: "pam", "kmeans", "dbscan", "optics", "mclust". Default is "pam".
#' Specifying "all" will run all clustering methods.
#' @param n_clust an \code{integer} specifying the requested number of cluster for pam
#' and kmeans.
#' @param find_nclust_method if n_clust is NULL, specificies the method to be
#' used to find the optimal number of clusters?  ##################################### TO BE IMPLEMENTED
#' @param pam_variant a \code{character} string specifying the variant of pam
#' to use, by default "faster". See \link[cluster:pam]{cluster::pam()} for
#' more details
#' @param pam_nstart an \code{integer} specifying the number of random “starts”
#' for the pam algorithm. By default, 1 (for the \code{"faster"} variant)
#' @param pam_cluster_only a \code{boolean} specifying if only the clustering
#' should be returned from the \link[cluster:pam]{cluster::pam()} function
#' (more efficient)
#' @param kmeans_iter_max an \code{integer} specifying the maximum number of
#' iterations for the kmeans method (see \link[stats:kmeans]{stats::kmeans()})
#' @param kmeans_nstart an \code{integer} specifying how many random sets of
#' \code{n_clust} should be selected as starting points for the kmeans analysis
#' (see \link[stats:kmeans]{stats::kmeans()})
#' @param dbscan_minPts a \code{numeric} value specifying the minPts argument
#' of \link[dbscan:dbscan]{dbscan::dbscan()}). By default, it is set to the
#' natural logarithm of the number of sites in \code{distances}.
#' @param dbscan_eps a \code{numeric} value specifying the eps argument
#' of \link[dbscan:dbscan]{dbscan::dbscan()}). The value of eps depends on the
#' minPts argument, and should be chosen by identifying a knee in the k-nearest
#' neighbour distance plot. By default the function will try to automatically
#' find a knee, but the result is uncertain, and so the user should inspect the
#' graph and modify \code{dbscan_eps} accordingly.
#' @param dbscan_plot a \code{boolean} indicating if the  k-nearest
#' neighbour distance plot should be plotted.
#' @param optics_xi a \code{numeric} value specifying
#'
#' @details
#'
#' \itemize{
#'  \item{\bold{Partitioning around medoids}: This method partitions data into
#'  the chosen number of cluster on the basis of the input distance matrix.
#'  It is more robust than k-means because it minimizes the sum of distances
#'  between cluster centres and points assigned to the cluster -
#'  whereas the k-means approach minimizes the sum of squared euclidean
#'  distances (thus k-means cannot be applied directly on the input distance
#'  matrix if the distances are not euclidean).}
#'  \item{\bold{k-means}: This method partitions the data into k groups
#'  such that that the sum of squares of euclidean distances from points to the
#'  assigned cluster centres is minimized. k-means cannot be applied directly
#'  on dissimilarity/beta-diversity metrics, because these distances are not
#'  euclidean. Therefore, it requires first to transform the distance matrix
#'  with a Principal Coordinate Analysis (using the function
#'  \link[ape:pcoa]{ape::pcoa()}), and then applying k-means on the coordinates
#'  of points in the PCoA.}
#'  \item{\bold{dbscan}: The dbscan (Density-based spatial clustering of
#'  applications with noise) clustering algorithm clusters points on the basis
#'  of the density of neighbours around each data points. It necessitates two
#'  main arguments, minPts, which stands for the minimum number of points to
#'  identify a core, and eps, which is the radius to find neighbours.
#'  minPts and eps should be defined by the user, which is not straightforward.
#'  We recommend reading the helpo in \link[dbscan:dbscan]{dbscan::dbscan()})
#'  to learn how to set these arguments, as well as the paper
#'  \insertCite{Hahsler2019}{bioRgeo}. Note that clusters with a value of 0
#'  are points which were deemed as noise by the algorithm.
#'   }
#'  \item{\bold{optics}: The optics (Ordering points to identify the clustering
#'  structure) is a semi-hierarchical clustering algorithm which orders the
#'  points in the dataset such that points which are closest become neighbours,
#'  and calculates a reachability distance for each point. Then, clusters
#'  can be extracted in a hierarchical manner from this reachability distance,
#'  by identifying clusters depending on changes in the relative cluster
#'  density. The reachability plot should be explored to understand
#'  the clusters and their hierarchical nature, by running plot on the output
#'  of the function: \code{plot(object$outputs$clustering_algorithms$optics)}.
#'  We recommend reading \insertCite{Hahsler2019}{bioRgeo} to grasp the
#'  algorithm, how it works, and what the clusters mean.}
#'  \item{\bold{mclust}: Model-based clustering based on Gaussian mixture
#'  models. This method is normally designed to be applied on the raw dataset,
#'  i.e. not directly on the distance matrix. However, we have found that it
#'  provides insightful results when applied on a distance matrix. Because this
#'  is an unusual application, we recommend users to be careful with the
#'  interpretation. See \insertCite{Scrucca2016}{bioRgeo} for more information
#'  on this method.}
#'
#'  }
#'
#' @return
#' to fill
#'
#' @export
#'
#' @examples
#' simil <- similarity(vegemat, metric = "all")
#' distances <- similarity_to_distance(simil)
#'
#' clust1 <- pam(distances,
#'     n_clust = 5,
#'     index = "Simpson")
#' clust2 <- pam(distances,
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
pam <- function(distances,
                index = names(distances)[3],
                n_clust = NULL,
                nstart = 10,
                variant = "faster", # c("original", "o_1", "o_2", "f_3", "f_4", "f_5", "faster")
                cluster_only = FALSE# , # To reduce computation time & memory, can be provided to cluster functions
                # find_nclust = FALSE,
                # find_nclust_metric = "pc_distance",
                # find_nclust_criterion = "elbow",
                # find_nclust_kmin = 2,
                # find_nclust_kmax = "number of sites",
                # find_nclust_sp_site_table = NULL
)
{
  if(inherits(distances, "bioRgeo.pairwise.metric"))
  {
    if(attr(distances, "type") == "similarity")
    {
      stop("distances seems to be a similarity object.
         clustering_nonhierarchical() should be applied on distances, not similarities.
         Use similarity_to_distance() before using clustering_nonhierarchical()")
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
        stop("n_clust must be an integer determining the number of clusters.")
      }
    } else
    {
      stop("n_clust must be an integer determining the number of clusters.")
    }
  }

  # if(is.null(n_clust)){
  #   # if(!(find_nclust_metric %in% c("pc_distance", "anosim", "avg_endemism",
  #   #                                "tot_endemism"))){
  #   #   stop("find_nclust_method must be one of pc_distance, anosim, avg_endemism or tot_endemism")
  #   # }
  #   # if(!(find_nclust_criterion %in% c("elbow",
  #   #                                   "increasing_step",
  #   #                                   "decreasing_step",
  #   #                                   "mars",
  #   #                                   "min",
  #   #                                   "max",
  #   #                                   "cutoff"))){
  #   #   stop("criterion must be one of elbow, increasing_step, decreasing_step, min, max, cutoff or mars")
  #   # }
  # }

  # Checker les index en input, ou créer une boucle pour toutes les implémenter
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
                       find_nclust_metric = find_nclust_metric,
                       find_nclust_criterion = find_nclust_criterion,
                       find_nclust_kmin = find_nclust_kmin,
                       find_nclust_kmax = find_nclust_kmax
  )

  outputs$inputs <- list(bipartite = FALSE,
                         weight = TRUE,
                         pairwise.metric = TRUE,
                         distance = TRUE,
                         nb_sites = attr(dist.obj, "Size"))

  outputs$algorithm <- list()

  outputs$clusters <- data.frame(matrix(ncol = 1,
                                        nrow = length(labels(dist.obj)),
                                        dimnames = list(labels(dist.obj),
                                                        "site")))

  outputs$clusters$site <- labels(dist.obj)


  # if(find_nclust) {
  #   if(find_nclust_kmax == "number of sites")
  #   {
  #     message("find_nclust_kmax set to nb_sites - 1, such that anosim can be computed")
  #     find_nclust_kmax <- outputs$inputs$nb_sites - 1
  #   } else if(!(find_nclust_kmax %% 1 == 0)) # integer testing ain't easy in R
  #   {
  #     stop("tree_k_max must be an integer determining the number of clusters.")
  #   }
  #
  #   outputs$algorithm$pam <- lapply(find_nclust_kmin:find_nclust_kmax,
  #                                   function(x)
  #                                     cluster::pam(dist.obj,
  #                                                  k = x,
  #                                                  diss = TRUE,
  #                                                  keep.diss = FALSE,
  #                                                  keep.data = FALSE,
  #                                                  nstart = nstart,
  #                                                  variant = variant,
  #                                                  cluster.only = cluster_only))
  #
  #   names(outputs$algorithm$pam) <- paste0("K_", find_nclust_kmin:find_nclust_kmax)
  #
  #   outputs$clusters <- data.frame(outputs$clusters,
  #                                  data.frame(lapply(names(outputs$algorithm$pam),
  #                                                    function(x)
  #                                                      outputs$algorithm$pam[[x]]$clustering)))
  #   outputs$clusters <- bioRgeo:::knbclu(outputs$clusters)
  #   class(outputs) <-  append("bioRgeo.clusters", class(outputs))
  #
  #   a <- partition_metrics(outputs,
  #                          distances = distances,
  #                          distance_index = index,
  #                          sp_site_table = find_nclust_sp_site_table,
  #                          eval_metric = find_nclust_metric,
  #                          partition_optimisation = TRUE,
  #                          criterion = find_nclust_criterion
  #                          )
  # } else {
    outputs$algorithm$pam <- lapply(n_clust,
           function(x)
             cluster::pam(dist.obj,
                          k = x,
                          diss = TRUE,
                          keep.diss = FALSE,
                          keep.data = FALSE,
                          nstart = nstart,
                          variant = variant,
                          cluster.only = cluster_only))

    names(outputs$algorithm$pam) <- paste0("K_", n_clust)

    outputs$clusters <- data.frame(outputs$clusters,
                                   data.frame(lapply(names(outputs$algorithm$pam),
                                                     function(x)
                                                       outputs$algorithm$pam[[x]]$clustering)))
    outputs$clusters <- bioRgeo:::knbclu(outputs$clusters)
    class(outputs) <-  append("bioRgeo.clusters", class(outputs))
  # }

  return(outputs)
}
