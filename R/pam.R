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
#' @param optics_xi a \code{numeric} value specifying
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
#' to fill
#'
#' @export
#'
#' @examples
#' simil <- similarity(vegemat, metric = "all")
#' distances <- similarity_to_distance(simil)
#'
#' clust1 <- pam(distances,
#'     n_clust = 2:10,
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
cl_pam <- function(distances,
                  index = names(distances)[3],
                  n_clust = NULL,
                  nstart = 10,
                  variant = "faster", # c("original", "o_1", "o_2", "f_3", "f_4", "f_5", "faster")
                  cluster_only = FALSE# To reduce computation time & memory, can be provided to cluster functions
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
                       cluster_only = cluster_only
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
