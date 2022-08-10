#' Hierarchical clustering based on dissimilarity or beta-diversity
#'
#' This function generates a hierarchical tree from a dissimilarity (beta-diversity)
#'  \code{data.frame},
#' calculates the cophenetic correlation coefficient,
#' and can get clusters from the tree if requested by the user. The function
#' implements randomization of the dissimilarity matrix to generate the tree, with
#' a selection method based on the optimal cophenetic correlation coefficient.
#' Typically, the dissimilarity \code{data.frame} is a \code{bioRgeo.pairwise.metric} object
#' obtained by running \code{similarity} or \code{similarity} and then 
#' \code{similarity_to_dissimilarity}.
#'
#' @param dissimilarity the output object from \code{\link{dissimilarity}} or
#'  \code{\link{similarity_to_dissimilarity}}, or a \code{dist} object. 
#'  If a \code{data.frame} is used, the first two 
#' columns represent pairs of sites (or any pair of nodes), and the next column(s)
#' are the dissimilarity indices. 
#' @param index name or number of the dissimilarity column to use. By default, 
#' the third column name of
#'  \code{dissimilarity} is used.
#' @param method name of the hierarchical classification method, as in
#' \link[stats:hclust]{stats::hclust()}. Should be one of This should be one of
#' \code{"ward.D"},
#' \code{"ward.D2"}, \code{"single"}, \code{"complete"}, \code{"average"}
#' (= UPGMA), \code{"mcquitty"} (= WPGMA),
#' \code{"median"} (= WPGMC) or \code{"centroid"} (= UPGMC).
#' @param randomize a boolean indicating if the dissimilarity matrix should be
#' randomized, to account for the order of sites in the dissimilarity matrix.
#' @param n_runs number of trials to randomize the dissimilarity matrix.
#' @param keep_trials a boolean indicating if all random trial results.
#' should be stored in the output object (set to FALSE to save space if your
#' \code{dissimilarity} object is large).
#' @param optimal_tree_method a character vector indicating how the final tree
#' should be obtained from all trials. The only option currently is
#' \code{"best"}, which means the tree with the best cophenetic correlation
#' coefficient will be chosen.
#' @param n_clust an integer or a vector of integers indicating the number of
#' clusters to be obtained from the hierarchical tree, or the output from
#' \link{partition_metrics}. Should not be used at the same time as
#' \code{cut_height}.
#' @param cut_height a numeric vector indicating the height(s) at which the tree
#' should be cut. Should not be used at the same time as \code{n_clust}.
#' @param find_h a boolean indicating if the height of cut should be found for
#' the requested \code{n_clust}.
#' @param h_max a numeric indicating the maximum possible tree height for
#' the chosen \code{index}.
#' @param h_min a numeric indicating the minimum possible height in the tree for
#' the chosen \code{index}.
#' @export
#' @details
#' The default method for the hierarchical tree is \code{"average"}, i.e.
#' UPGMA as it has been recommended as the best method to generate a tree
#' from beta diversity dissimilarity \insertCite{Kreft2010}{bioRgeo}
#'
#' Clusters can be obtained by two methods:
#' \itemize{
#' \item{Specifying a desired number of clusters in \code{n_clust}}
#' \item{Specifying one or several heights of cut in \code{cut_height}}}
#'
#' To find an optimal number of clusters, see \code{\link{partition_metrics}}
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
#' In the \code{algorithm} slot, users can find the following elements:
#'
#' \itemize{
#' \item{\code{trials}: a list containing all randomization trials. Each trial
#' containes the dissimilarity matrix, with site order randomized, the associated
#' tree and the cophenetic correlation coefficient (spearman) for that tree}
#' \item{\code{final.tree}: a \code{hclust} object containing the final
#' hierarchical tree to be used}
#' \item{\code{final.tree.coph.cor}: the cophenetic correlation coefficient
#' between the initial dissimilarity matrix and \code{final.tree}}
#' }
#'
#' @references
#' \insertRef{Kreft2010}{bioRgeo}
#' @author
#' Boris Leroy (\email{leroy.boris@gmail.com}),
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) and
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) 
#' 
#' @seealso \link{cut_tree} 
#' @examples
#' comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001), 20, 25)
#' rownames(comat) <- paste0("Site",1:20)
#' colnames(comat) <- paste0("Species",1:25)
#'
#' dissim <- dissimilarity(comat, metric = "all")
#'
#' # User-defined number of clusters
#' tree1 <- hclu_hierarclust(dissim,
#'                                  n_clust = 5)
#' tree1
#' plot(tree1)
#' str(tree1)
#' tree1$clusters
#' 
#' # User-defined height cut
#' # Only one height
#' tree2 <- hclu_hierarclust(dissim,
#'                                  cut_height = .05)
#' tree2
#' tree2$clusters
#' #
#' # Multiple heights
#' tree3 <- hclu_hierarclust(dissim,
#'                                  cut_height = c(.05, .15, .25))
#' tree3
#' tree3$clusters # Mind the order of height cuts: from deep to shallow cuts
#' plot(tree3)
#' 
#' # Recut the tree afterwards
#' tree3.1 <- cut_tree(tree3,
#'                  n = 5)
hclu_hierarclust <- function(dissimilarity,
                             index = names(dissimilarity)[3],
                             method = "average",
                             randomize = TRUE,
                             n_runs = 30,
                             keep_trials = FALSE,
                             optimal_tree_method = "best", # best or consensus
                             n_clust = NULL,
                             cut_height = NULL,
                             find_h = TRUE,
                             h_max = 1,
                             h_min = 0)
{

  if(inherits(dissimilarity, "bioRgeo.pairwise.metric"))
  {
    if(attr(dissimilarity, "type") == "similarity")
    {
      stop("dissimilarity seems to be a similarity object.
         hclu_hierarclust() should be applied on dissimilarity, not similarities.
         Use similarity_to_dissimilarity() before using hclu_hierarclust()")
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
      if(any(!(n_clust %% 1 == 0)))
      {
        stop("n_clust must an integer or a vector of integers determining the number of clusters.")
      }

    } else if(inherits(n_clust, "bioRgeo.partition.metrics"))
    {
      if(!is.null(n_clust$algorithm$optimal_nb_clusters)) {
        n_clust <- n_clust$algorithm$optimal_nb_clusters
      } else {
        stop("n_clust does not have an optimal number of clusters. Did you specify
             partition_optimisation = TRUE in partition_metrics()?")
      }
    } else
    {
      stop("n_clust must be one of those:
        * an integer determining the number of clusters
        * a vector of integers determining the numbers of clusters for each cut
        * the output from partition_metrics()")
    }
    if(!is.null(cut_height))
    {
      stop("Please provide either n_clust or cut_height, but not both at the same time.")
    }
  }


  outputs <- list(name = "hierarchical_clustering")

  if(!inherits(dissimilarity, "dist"))
  {
      # dist.obj <- .dfToDist(dissimilarity, metric = index)
    dist.obj <- stats::as.dist(
      net_to_mat(dissimilarity[, c(colnames(dissimilarity)[1:2], index)],
                        weight = TRUE, squared = TRUE, symmetrical = TRUE))

  } else {
    dist.obj <- dissimilarity
  }

  outputs$args <- list(index = index,
                       method = method,
                       randomize = randomize,
                       n_runs = n_runs,
                       optimal_tree_method = optimal_tree_method,
                       n_clust = n_clust,
                       cut_height = cut_height,
                       find_h = find_h,
                       h_max = h_max,
                       h_min = h_min
  )

  outputs$inputs <- list(bipartite = FALSE,
                         weight = TRUE,
                         pairwise_metric = TRUE,
                         dissimilarity = TRUE,
                         nb_sites = attr(dist.obj, "Size"))

  # outputs$dist.matrix <- dist.obj

  if(randomize)
  {
    message(paste0("Randomizing the dissimilarity matrix with ", n_runs, " trials"))
    for(run in 1:n_runs)
    {
      outputs$algorithm$trials[[run]] <- list()
      outputs$algorithm$trials[[run]]$dist.matrix <- .randomizeDistance(dist.obj)

      # Compute hierarchical tree
      outputs$algorithm$trials[[run]]$hierartree <- fastcluster::hclust(outputs$algorithm$trials[[run]]$dist.matrix,
                                                                     method = method)

      # Calculate cophenetic correlation coefficient
      coph <- as.matrix(stats::cophenetic(outputs$algorithm$trials[[run]]$hierartree))
      # Correct site order to not mess the correlation
      coph <- coph[match(attr(outputs$algorithm$trials[[run]]$dist.matrix, "Labels"),
                         rownames(coph)),
                   match(attr(outputs$algorithm$trials[[run]]$dist.matrix, "Labels"),
                         colnames(coph))]
      dist.mat <- as.matrix(outputs$algorithm$trials[[run]]$dist.matrix)

      outputs$algorithm$trials[[run]]$cophcor <- stats::cor(dist.mat[lower.tri(dist.mat)],
                                                  coph[lower.tri(coph)],
                                                  method = "spearman")
    }

    if(optimal_tree_method == "best")
    {
      coph.coeffs <- sapply(1:n_runs, function(x)
      {
        outputs$algorithm$trials[[x]]$cophcor
      })

      message(paste0(" -- range of cophenetic correlation coefficients among trials: ", round(min(coph.coeffs), 2), " - ", round(max(coph.coeffs), 2)))

      best.run <- which(coph.coeffs == max(coph.coeffs))[1] # There might be multiple trees with the highest cophenetic correlation coefficient, so we arbritrarily take the first one

      outputs$algorithm$final.tree <- outputs$algorithm$trials[[best.run]]$hierartree
      outputs$algorithm$final.tree.coph.cor <- max(coph.coeffs)
    }

    message(paste0("Optimal tree has a ", round(outputs$algorithm$final.tree.coph.cor, 2), " cophenetic correlation coefficient with the initial dissimilarity matrix\n"))

  } else
  {
    outputs$algorithm$final.tree <- fastcluster::hclust(dist.obj, method = method)

    coph <- as.matrix(stats::cophenetic(outputs$algorithm$final.tree))
    coph <- coph[match(attr(dist.obj, "Labels"),
                       rownames(coph)),
                 match(attr(dist.obj, "Labels"),
                       colnames(coph))]
    dist.mat <- as.matrix(dist.obj)

    outputs$algorithm$final.tree.coph.cor <- stats::cor(dist.mat[lower.tri(dist.mat)],
                                              coph[lower.tri(coph)],
                                              method = "spearman")

    message(paste0("Output tree has a ", round(outputs$algorithm$final.tree.coph.cor, 2), " cophenetic correlation coefficient with the initial dissimilarity matrix\n"))
  }


  class(outputs) <- append("bioRgeo.clusters", class(outputs))

  if(any(!is.null(n_clust) | !is.null(cut_height)))
  {
    outputs <- cut_tree(outputs,
                        n_clust = n_clust,
                        cut_height = cut_height,
                        find_h = find_h,
                        h_max = h_max,
                        h_min = h_min)
  } else{
    outputs$clusters <- NA
  }

  if(!keep_trials)
  {
    outputs$algorithm$trials <- "Trials not stored in output"
  }

  return(outputs)
}


.randomizeDistance <- function(distmatrix)
{
  distmatrix <- as.matrix(distmatrix)
  ord <- sample(rownames(distmatrix))
  return(stats::as.dist(distmatrix[ord, ord]))
}
