#' Non hierarchical clustering: partitioning around medoids
#'
#' This function performs non hierarchical
#' clustering on the basis of dissimilarity with partitioning around medoids.
#'
#' @param dissimilarity the output object from \code{\link{dissimilarity}} or
#'  \code{\link{similarity_to_dissimilarity}}, or a \code{dist} object. 
#'  If a \code{data.frame} is used, the first two 
#' columns represent pairs of sites (or any pair of nodes), and the next column(s)
#' are the dissimilarity indices. 
#' @param index name or number of the dissimilarity column to use. By default, 
#' the third column name of
#'  \code{dissimilarity} is used.
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
#'  the chosen number of cluster on the basis of the input dissimilarity matrix.
#'  It is more robust than k-means because it minimizes the sum of dissimilarity
#'  between cluster centres and points assigned to the cluster -
#'  whereas the k-means approach minimizes the sum of squared euclidean
#'  distances (thus k-means cannot be applied directly on the input dissimilarity
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
#' @author
#' Boris Leroy (\email{leroy.boris@gmail.com}),
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) and
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) 
#' 
#' @seealso \link{nhclu_kmeans} 
#' @export
#'
#' @examples
#' dissim <- dissimilarity(vegemat, metric = "all")
#' 
#' clust1 <- nhclu_pam(dissim,
#'      n_clust = 2:10,
#'      index = "Simpson")
#' clust2 <- nhclu_pam(dissim,
#'      n_clust = 2:25,
#'      index = "Simpson")
#' partition_metrics(clust2,
#'                    dissimilarity = dissim,
#'                    eval_metric = "pc_distance")
#' partition_metrics(clust2,
#'                    net = vegedf,
#'                    species_col = "Species",
#'                    site_col = "Site",
#'                    eval_metric = "avg_endemism")
nhclu_pam <- function(dissimilarity,
                      index = names(dissimilarity)[3],
                      n_clust = NULL,
                      nstart = if(variant == "faster") 1 else NA,
                      variant = "faster", # c("original", "o_1", "o_2", "f_3", "f_4", "f_5", "faster")
                      cluster_only = FALSE,# To reduce computation time & memory, can be provided to cluster functions
                      ... # Further arguments to be passed to cluster::pam
)
{
  if(inherits(dissimilarity, "bioRgeo.pairwise.metric"))
  {
    if(attr(dissimilarity, "type") == "similarity")
    {
      stop("dissimilarity seems to be a similarity object.
         nhclu_pam() should be applied on dissimilarity, not similarities.
         Use similarity_to_dissimilarity() before using nhclu_pam()")
    }
    if(is.numeric(index))
    {
      index <- names(dissimilarity)[index]
    }
    if(!(index %in% colnames(dissimilarity)))
    {
      stop("Argument index should be one of the column names of dissimilarity")
    }

  } else if(!any(inherits(dissimilarity, "bioRgeo.pairwise.metric"), inherits(dissimilarity, "dist")))
  {
    if(is.numeric(index))
    {
      index <- names(dissimilarity)[index]
    }
    if(!(index %in% colnames(dissimilarity)))
    {
      stop("dissimilarity is not a bioRgeo.pairwise.metric object, a dissimilarity matrix (class dist) or a data.frame with at least 3 columns (site1, site2, and your dissimilarity index)")
    }
  }

  if(!is.null(n_clust)){
    if(is.numeric(n_clust))
    {
      if(any(!(n_clust %% 1 == 0))) # integer testing ain't easy in R
      {
        stop("n_clust must an integer or a vector of integers determining the number of clusters.")
      }
      # Add test to see if n_clust is lower than the number of sites
    } else
    {
      stop("n_clust must an integer or a vector of integers determining the number of clusters.")
    }
  }

  if(!inherits(dissimilarity, "dist"))
  {
    dist.obj <- stats::as.dist(
      net_to_mat(dissimilarity[, c(colnames(dissimilarity)[1:2], index)],
                 weight = TRUE, squared = TRUE, symmetrical = TRUE))

  } else {
    dist.obj <- dissimilarity
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
                         dissimilarity = TRUE,
                         nb_sites = attr(dist.obj, "Size"))

  outputs$algorithm <- list()

  outputs$clusters <- data.frame(matrix(ncol = 1,
                                        nrow = length(labels(dist.obj)),
                                        dimnames = list(labels(dist.obj),
                                                        "name")))

  outputs$clusters$name <- labels(dist.obj)

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
  outputs$clusters <- knbclu(outputs$clusters, reorder = FALSE)
  
  outputs$cluster_info <- data.frame(partition_name = names(outputs$clusters)[2:length(outputs$clusters), drop = FALSE],
                                     n_clust = apply(outputs$clusters[, 2:length(outputs$clusters), drop = FALSE],
                                                     2, function(x) length(unique(x))))
  
  class(outputs) <-  append("bioRgeo.clusters", class(outputs))

  return(outputs)
}
