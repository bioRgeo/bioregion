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
#' @param nclust_find_method if n_clust is NULL, specificies the method to be
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
#' comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001), 20, 25)
#' rownames(comat) <- paste0("Site",1:20)
#' colnames(comat) <- paste0("Species",1:25)
#'
#' simil <- similarity(comat, metric = "all")
#' distances <- similarity_to_distance(simil)
#'
#' # User-defined number of clusters
#' # clust1 <- clustering_nonhierarchical(distances,
#' #                                    n_clust = 5,
#' #                                     index = "Simpson")
#'
#' #  clust2 <- clustering_nonhierarchical(distances,
#' #                                     method = c("pam", "kmeans", "dbscan",
#' #                                                 "optics", "mclust"),
#' #                                     n_clust = 5,
#' #                                     index = "Simpson",
#' #                                     dbscan_minPts = 5,
#' #                                     dbscan_eps = 0.06)
clustering_nonhierarchical <- function(distances,
                                       index = names(distances)[3],
                                       method = "pam",
                                       n_clust = NULL,
                                       nclust_find_method = "firstSEmax",
                                       pam_nstart = 10,
                                       pam_variant = "faster", # c("original", "o_1", "o_2", "f_3", "f_4", "f_5", "faster")
                                       pam_cluster_only = FALSE, # To reduce computation time & memory, can be provided to cluster functions
                                       kmeans_iter_max = 10,
                                       kmeans_nstart = 10,
                                       dbscan_minPts = NULL,
                                       dbscan_eps = NULL,
                                       dbscan_plot = TRUE,
                                       optics_xi = 0.05
                                       )
{
  if(!(all(method %in% c("pam", "kmeans", "dbscan",
                         "mclust", "optics")))){
    stop("The chosen clustering method is not available.
     Please chose among the followings:
         'pam', 'kmeans', 'dbscan', 'optic', 'mclust'")
  }

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
      if(!(n_clust %% 1 == 0)) # Integer testing
      {
        stop("n_clust must be an integer determining the number of clusters.")
      }
    } else
    {
      stop("n_clust must be an integer determining the number of clusters.")
    }
  }

  if(is.null(n_clust) & any(
    method %in%
    c("kmeans", "diana", "pam"))){
    if(!(nclust_find_method %in% c("globalmax", "firstmax", "Tibs2001SEmax",
                             "firstSEmax", "globalSEmax."))){
      stop("Chosen gap statistic to determine the optimal number of cluster is
      not available.
       Please chose among the followings:
           globalmax, firstmax, Tibs2001SEmax, firstSEmax or globalSEmax.")
    } else
    {
      message("   - finding an optimal number of clusters")
    }

  }



  if(!(kmeans_nstart %% 1 == 0)){
    stop("nstart must be an integer determining the number of random centroids
           to start k-means analysis.")
  }

  # Unused arguments
  # if(!(kmeans_B %% 1 == 0)){
  #   stop("B must be an integer determining the number of Monte Carlo bootstrap
  #          samples.")
  # }
  # if(!(kmeans_Kmax %% 1 == 0)){
  #   stop("K.max must be a numeric determining the maximum number of clusters
  #          to consider.")
  # }
  #   if (kmeans_Kmax > length(labels(dist.obj))) {
  #     warning("kmeans_Kmax should not be superior to the number of sites,
  # reducing to number of sites.")
  #     kmeans_Kmax <- length(labels(dist.obj))
  #   }


  outputs <- list()

  # Checker les index en input, ou créer une boucle pour toutes les implémenter
  if(!inherits(distances, "dist"))
  {
    # dist.obj <- .dfToDist(distances, metric = index)
    dist.obj <- stats::as.dist(
      net_to_mat(distances[, c(1, 2,
                                      which(colnames(distances) == index))],
                        weight = TRUE, squared = TRUE, symmetrical = TRUE))

  }



  # COMPLETER LES ARGUMENTS A LA FIN
  outputs$args <- list(index = index,
                       method = method,
                       n_clust = n_clust,
                       nclust_find_method = nclust_find_method,
                       pam_nstart = pam_nstart,
                       pam_variant = pam_variant, # c("original", "o_1", "o_2", "f_3", "f_4", "f_5", "faster")
                       pam_cluster_only = pam_cluster_only, # To reduce computation time & memory, can be provided to cluster functions
                       kmeans_iter_max = kmeans_iter_max,
                       kmeans_nstart = kmeans_nstart,
                       dbscan_minPts = dbscan_minPts,
                       dbscan_eps = dbscan_eps,
                       dbscan_plot = dbscan_plot,
                       optics_xi = optics_xi
  )


  outputs$dist.matrix <- dist.obj


  # dist.matrix <- as.matrix(dist.obj)

  outputs$clustering_algorithms <- list()

  outputs$clusters <- data.frame(matrix(ncol = length(method),
                                        nrow = length(labels(dist.obj)),
                                        dimnames = list(labels(dist.obj),
                                                        method)))

  if ("pam" %in% method) {
    outputs$clustering_algorithms$pam <-
      cluster::pam(dist.obj, k = n_clust,
                   diss = TRUE,
                   keep.diss = FALSE,
                   keep.data = FALSE,
                   nstart = pam_nstart,
                   variant = pam_variant,
                   cluster.only = pam_cluster_only)
    if(pam_cluster_only)
    {
      outputs$clusters$pam <- outputs$clustering_algorithms$pam
    } else
    {
      outputs$clusters$pam <- outputs$clustering_algorithms$pam$clustering
    }
  }

  if ("kmeans" %in% method) {
    # kmeans only works on Euclidean distances, so the distance matrix needs to
    # be transformed into a multivariate space with euclidean distances
    # with a Principal Coordinate Analysis
    outputs$clustering_algorithms$pcoa <- ape::pcoa(dist.obj)

    # Performing the kmeans on the PCoA with all axes
    outputs$clustering_algorithms$kmeans <-
      stats::kmeans(
        outputs$clustering_algorithms$pcoa$vectors,
        centers = n_clust,
        iter.max = kmeans_iter_max,
        nstart = kmeans_nstart
      )

    outputs$clusters$kmeans <- outputs$clustering_algorithms$kmeans$cluster
  }


  if("dbscan" %in% method)
  {
    if(is.null(dbscan_minPts))
    {
      # Using a default value of minPts if none provided by the user
      dbscan_minPts <- log(length(labels(dist.obj)))
    }


    knnp <- dbscan::kNNdist(dist.obj,
                            k = dbscan_minPts - 1)


    x_ <- order(knnp)

    # Trying to find the knee, and not the elbow
    knee <- .elbow_finder(x_, max(knnp) - knnp)

    if (is.null(dbscan_eps)) {
      dbscan_eps <- knee[2]
    }

    if(dbscan_plot)
    {
      plot(knnp[order(knnp)], type = "l",
           main = "dbscan parameter choice:\nchoose eps where there is a knee in the curve",
           xlab = "", ylab = "epsilon")
      graphics::abline(h = dbscan_eps)
      graphics::text(x = 0,
           y = knee[2] + 0.025 * (max(knnp) - min(knnp)),
           labels = paste0("Chosen eps value: ",
                           round(dbscan_eps, 3)),
           adj = 0)
    }

    outputs$clustering_algorithms$dbscan <-
      dbscan::dbscan(dist.obj,
                     minPts = dbscan_minPts,
                     eps = dbscan_eps)

    # NOTE: values of 0 mean "noise / no cluster" with dbscan
    outputs$clusters$dbscan <-
      outputs$clustering_algorithms$dbscan$cluster

  }


  if("optics" %in% method)
  {
    if(is.null(dbscan_minPts))
    {
      # Using a default value of minPts if none provided by the user
      dbscan_minPts <- log(length(labels(dist.obj)))
    }

    outputs$clustering_algorithms$optics <- dbscan::optics(x = dist.obj,
                                              minPts = dbscan_minPts)
    outputs$clustering_algorithms$optics <-
      dbscan::extractXi(outputs$clustering_algorithms$optics,
                xi = optics_xi)

    outputs$clusters$optics <-
      outputs$clustering_algorithms$optics$cluster

  }

  if("mclust" %in% method)
  {
    mclustBIC <- mclust::mclustBIC
    outputs$clustering_algorithms$mclust <-
      mclust::Mclust(dist.obj)
    outputs$clusters$mclust <-
      outputs$clustering_algorithms$mclust$classification
  }

  class(outputs) <- append("bioRgeo.nonhierar.cluster", class(outputs))


  return(outputs)
}
