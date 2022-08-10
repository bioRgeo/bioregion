#' Non hierarchical clustering: k-means analysis
#'
#' This function performs non hierarchical
#' clustering on the basis of dissimilarity with a k-means analysis.
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
#' @param iter_max an \code{integer} specifying the maximum number of
#' iterations for the kmeans method (see \link[stats:kmeans]{stats::kmeans()})
#' @param nstart an \code{integer} specifying how many random sets of
#' \code{n_clust} should be selected as starting points for the kmeans analysis
#' (see \link[stats:kmeans]{stats::kmeans()})
#' @param algorithm a \code{character string} specifying the algorithm to use for
#' kmean (see \link[stats:kmeans]{stats::kmeans()})
#'
#' @details
#' This method partitions the data into k groups
#'  such that that the sum of squares of euclidean distances from points to the
#'  assigned cluster centres is minimized. k-means cannot be applied directly
#'  on dissimilarity/beta-diversity metrics, because these distances are not
#'  euclidean. Therefore, it requires first to transform the dissimilarity matrix
#'  with a Principal Coordinate Analysis (using the function
#'  \link[ape:pcoa]{ape::pcoa()}), and then applying k-means on the coordinates
#'  of points in the PCoA. Because this makes an additional transformation of
#'  the initial matrix of dissimilarity, the partitioning around medoids method
#'  should be prefered (\code{\link{nhclu_pam}})
#'
#' @seealso \link{nhclu_pam} 
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
#' @author
#' Boris Leroy (\email{leroy.boris@gmail.com}),
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) and
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) 
#' 
#' @seealso \link{cut_tree} 
#' @examples
#' dissim <- dissimilarity(vegemat, metric = "all")
#' 
#' clust1 <- nhclu_kmeans(dissim,
#'      n_clust = 2:10,
#'      index = "Simpson")
#' clust2 <- nhclu_kmeans(dissim,
#'      n_clust = 2:25,
#'      index = "Simpson")
#' partition_metrics(clust2,
#'                    dissimilarity = dissim,
#'                    eval_metric = "pc_distance",
#'                    partition_optimisation = TRUE)
#' \donttest{
#' partition_metrics(clust2,
#'                    sp_site_table = vegemat,
#'                    eval_metric = "avg_endemism",
#'                    partition_optimisation = TRUE)
#'                    }
nhclu_kmeans <- function(dissimilarity,
                         index = names(dissimilarity)[3],
                         n_clust = NULL,
                         iter_max = 10,
                         nstart = 10,
                         algorithm = "Hartigan-Wong"
)
{
  if(inherits(dissimilarity, "bioRgeo.pairwise.metric"))
  {
    if(attr(dissimilarity, "type") == "similarity")
    {
      stop("dissimilarity seems to be a similarity object.
         nhclu_kmeans() should be applied on dissimilarity, not similarities.
         Use similarity_to_dissimilarity() before using nhclu_kmeans()")
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



  outputs <- list(name = "kmeans")

  outputs$args <- list(index = index,
                       n_clust = n_clust,
                       iter_max = iter_max,
                       nstart = nstart,
                       algorithm = algorithm
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


  # kmeans only works on Euclidean distances, so the dissimilarity matrix needs to
  # be transformed into a multivariate space with euclidean distances
  # with a Principal Coordinate Analysis
  outputs$clustering_algorithms$pcoa <- ape::pcoa(dist.obj)

  # Performing the kmeans on the PCoA with all axes
  outputs$algorithm$kmeans <- lapply(n_clust,
                                     function(x)
                                       stats::kmeans(dist.obj,
                                                     centers = x,
                                                     iter.max = iter_max,
                                                     nstart = nstart,
                                                     algorithm = algorithm))


  names(outputs$algorithm$kmeans) <- paste0("K_", n_clust)

  outputs$clusters <- data.frame(outputs$clusters,
                                 data.frame(lapply(names(outputs$algorithm$kmeans),
                                                   function(x)
                                                     outputs$algorithm$kmeans[[x]]$cluster)))
  outputs$clusters <- knbclu(outputs$clusters,
                             reorder = FALSE)
  
  outputs$cluster_info <- data.frame(partition_name = names(outputs$clusters)[2:length(outputs$clusters), drop = FALSE],
                                     n_clust = apply(outputs$clusters[, 2:length(outputs$clusters), drop = FALSE],
                                                     2, function(x) length(unique(x))))
  
  class(outputs) <-  append("bioRgeo.clusters", class(outputs))

  return(outputs)
}
