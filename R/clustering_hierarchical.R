#' Hierarchical clustering from distances or beta-diversity
#'
#' This function generates a hierarchical tree from a distance (beta-diversity)
#'  \code{data.frame},
#' calculates the cophenetic correlation coefficient,
#' and can get clusters from the tree if requested by the user. The function
#' implements randomization of the distance matrix to generate the tree, with
#' a selection method based on the optimal cophenetic correlation coefficient.
#' Typically, the distance \code{data.frame} is a \code{bioRgeo.distance} object
#' obtained by running \code{spproject} and then \code{similarity_to_distance}.
#'
#' @param distances the output object from \code{\link{similarity_to_distance}}
#' or a \code{data.frame} with the first columns called "Site1" and "Site2", and
#' the other columns being the distance indices.
#' @param index name of the distance index to use, corresponding to the column
#' name in \code{distances}. By default, the third column name of
#'  \code{distances} is used.
#' @param method name of the hierarchical classification method, as in
#' \link[stats:hclust]{stats::hclust()}. Should be one of This should be one of
#' \code{"ward.D"},
#' \code{"ward.D2"}, \code{"single"}, \code{"complete"}, \code{"average"}
#' (= UPGMA), \code{"mcquitty"} (= WPGMA),
#' \code{"median"} (= WPGMC) or \code{"centroid"} (= UPGMC).
#' @param randomize a boolean indicating if the distance matrix should be
#' randomized, to account for the order of sites in the distance matrix
#' @param n_runs number of trials to randomize the distance matrix
#' @param optimal_tree_method a character vector indicating how the final tree
#' should be obtained from all trials. The only option currently is
#' \code{"best"}, which means the tree with the best cophenetic correlation
#' coefficient will be chosen.
#' @param n_clust an integer indicating the number of clusters to be obtained
#' from the hierarchical tree. Should not be used at the same time as
#' \code{cut_height} or \code{optim_method}
#' @param cut_height a numeric vector indicating the height(s) at which the tree
#' should be cut. Should not be used at the same time as \code{n_clust} or
#' \code{optim_method}
#' @param find_h a boolean indicating if the height of cut should be found for
#' the requested \code{n_clust}
#' @param h_max a numeric indicating the maximum possible tree height for
#' the chosen \code{index}
#' @param h_min a numeric indicating the minimum possible height in the tree for
#' the chosen \code{index}
#' @param optim_method a character vector indicating the method to find the
#' optimal number of clusters in the tree. Should not be used at the same time
#' as \code{n_clust} or \code{cut_height}
#' @export
#' @details
#' The default method for the hierarchical tree is \code{"average"}, i.e.
#' UPGMA as it has been recommended as the best method to generate a tree
#' from beta diversity distances \insertCite{Kreft2010}{bioRgeo}
#'
#' Clusters can be obtained by three methods:
#' \itemize{
#' \item{Specifying a desired number of clusters in \code{n_clust}}
#' \item{Specifying one or several heights of cut in \code{cut_height}}
#' \item{Specifying an optimization method in \code{optim_method} to find the
#' optimal number of clusters}}
#'
#'
#'
#' @return A \code{list} with additional class \code{bioRgeo.hierar.tree}
#' \itemize{
#' \item{\code{args}: the input arguments}
#' \item{\code{dist.matrix}: the distance/beta diversity matrix as a
#'  \code{dist} object}
#' \item{\code{trials}: a list containing all randomization trials. Each trial
#' containes the distance matrix, with site order randomized, the associated
#' tree and the cophenetic correlation coefficient for that tree}
#' \item{\code{final.tree}: a \code{hclust} object containing the final
#' hierarchical tree to be used}
#' \item{\code{final.tree.coph.cor}: the cophenetic correlation coefficient
#' between the initial distance matrix and \code{final.tree}}
#' \item{\code{clusters}: a \code{data.frame} containing the clusters}
#' }
#'
#' @references
#' \insertRef{Kreft2010}{bioRgeo}
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}),
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) and
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' @seealso \link{distance_to_similarity}
#' @examples
#' comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001), 20, 25)
#' rownames(comat) <- paste0("Site",1:20)
#' colnames(comat) <- paste0("Species",1:25)
#'
#' simil <- spproject(comat, metric = "all")
#' distances <- similarity_to_distance(simil)
#'
#' # User-defined number of clusters
#' tree1 <- clustering_hierarchical(distances,
#'                                  n_clust = 5)
#' tree1
#' plot(tree1)
#' str(tree1)
#' tree1$clusters
#'
#' # User-defined height cut
#' # Only one height
#' tree2 <- clustering_hierarchical(distances,
#'                                  cut_height = .05)
#' tree2
#' tree2$clusters
#'
#' # Multiple heights
#' tree3 <- clustering_hierarchical(distances,
#'                                  cut_height = c(.05, .15, .25))
#' tree3
#' tree3$clusters # Mind the order of height cuts: from deep to shallow cuts
#'
#' # Recut the tree afterwards
#' tree3 <- cut_tree(tree3,
#'                   n = 5)
clustering_hierarchical <- function(distances,
                                    index = names(distances)[3],
                                    method = "average",
                                    randomize = TRUE,
                                    n_runs = 100,
                                    optimal_tree_method = "best", # best or consensus
                                    n_clust = NULL,
                                    cut_height = NULL,
                                    find_h = TRUE,
                                    h_max = 1,
                                    h_min = 0,
                                    optim_method = "firstSEmax")
{

  # Ajouter la possibilité d'avoir une matrice de distance à la place de dist.df

  if(inherits(distances, "bioRgeo.similarity"))
  {
    warning("distances seems to be a bioRgeo.similarity object. clusterHierarch should be applied on distances, not similarities. Consider using similarity_to_distance() before using clusterHierarch")
  }

  if(!is.null(n_clust)){
    if(is.numeric(n_clust))
    {
      if(!(n_clust %% 1 == 0)) # integer testing ain't easy in R
      {
        stop("n_clust must be an integer determining the number of clusters.")
      }
    } else
    {
      stop("n_clust must be an integer determining the number of clusters.")
    }
    if(!is.null(cut_height))
    {
      stop("Please provide either n_clust or cut_height, but not both at the same time.")
    }
  }

  if(!is.null(cut_height))
  {
    if(!is.numeric(cut_height))
    {
      stop("cut.height must be a numeric determing the height(s) of cut on the hierarchical tree.")
    }
  }


  outputs <- list()

  # Checker les index en input, ou créer une boucle pour toutes les implémenter
  dist.obj <- .dfToDist(distances, metric = index)


  outputs$args <- list(index = index,
                       method = method,
                       randomize = randomize,
                       n_runs = n_runs,
                       optimal_tree_method = optimal_tree_method,
                       n_clust = n_clust,
                       cut_height = cut_height,
                       find_h = find_h,
                       h_max = h_max,
                       h_min = h_min,
                       optim_method = optim_method
  )

  outputs$dist.matrix <- dist.obj

  if(randomize)
  {
    message(paste0("Randomizing the distance matrix with ", n_runs, " trials"))
    for(run in 1:n_runs)
    {
      outputs$trials[[run]] <- list()
      outputs$trials[[run]]$dist.matrix <- .randomizeDistance(dist.obj)

      # Compute hierarchical tree
      outputs$trials[[run]]$hierartree <- fastcluster::hclust(outputs$trials[[run]]$dist.matrix,
                                                                     method = method)

      # Calculate cophenetic correlation coefficient
      coph <- as.matrix(cophenetic(outputs$trials[[run]]$hierartree))
      # Correct site order to not mess the correlation
      coph <- coph[match(attr(outputs$trials[[run]]$dist.matrix, "Labels"),
                         rownames(coph)),
                   match(attr(outputs$trials[[run]]$dist.matrix, "Labels"),
                         colnames(coph))]
      dist.mat <- as.matrix(outputs$trials[[run]]$dist.matrix)

      outputs$trials[[run]]$cophcor <- cor(dist.mat[lower.tri(dist.mat)],
                                                  coph[lower.tri(coph)])
    }

    if(optimal_tree_method == "best")
    {
      coph.coeffs <- sapply(1:n_runs, function(x)
      {
        outputs$trials[[x]]$cophcor
      })

      message(paste0(" -- range of cophenetic correlation coefficients among trials: ", round(min(coph.coeffs), 2), " - ", round(max(coph.coeffs), 2)))

      best.run <- which(coph.coeffs == max(coph.coeffs))[1] # There might be multiple trees with the highest cophenetic correlation coefficient, so we arbritrarily take the first one

      outputs$final.tree <- outputs$trials[[best.run]]$hierartree
      outputs$final.tree.coph.cor <- max(coph.coeffs)
    }

    message(paste0("Optimal tree has a ", round(outputs$final.tree.coph.cor, 2), " cophenetic correlation coefficient with the initial distance matrix\n"))

  } else
  {
    outputs$final.tree <- fastcluster::hclust(outputs$dist.matrix, method = method)

    coph <- as.matrix(cophenetic(outputs$final.tree))
    coph <- coph[match(attr(dist.obj, "Labels"),
                       rownames(coph)),
                 match(attr(dist.obj, "Labels"),
                       colnames(coph))]
    dist.mat <- as.matrix(dist.obj)

    outputs$final.tree.coph.cor <- cor(dist.mat[lower.tri(dist.mat)],
                                              coph[lower.tri(coph)])

    message(paste0("Output tree has a ", round(outputs$final.tree.coph.cor, 2), " cophenetic correlation coefficient with the initial distance matrix\n"))
  }

  outputs$clusters <- cut_tree(outputs$final.tree,
                               n_clust = n_clust,
                               cut_height = cut_height,
                               find_h = find_h,
                               h_max = h_max,
                               h_min = h_min)

  class(outputs) <- append("bioRgeo.hierar.tree", class(outputs))
  return(outputs)
}


# Internal functions for clustering_hierarchical
.dfToDist <- function(distancedf, metric)
{
  if(inherits(distancedf, "bioRgeo.distance"))
  {
    other.cols <- colnames(distancedf)[-which(colnames(distancedf) %in% c("Site1", "Site2"))]
  } else
  {
    stop("distancedf should be a distance data.frame of class bioRgeo.distance")
  }

  nodes <- unique(c(distancedf$Site1, distancedf$Site2))

  distancedf <- rbind(data.frame(Site1 = nodes,
                                 Site2 = nodes,
                                 matrix(data = 0,
                                        nr = length(unique(c(distancedf$Site1, distancedf$Site2))),
                                        nc = length(other.cols),
                                        dimnames = list(NULL,
                                                        other.cols))),
                      distancedf)
  distancematrix <- as.dist(xtabs(distancedf[, metric] ~ distancedf$Site2 + distancedf$Site1))

  return(distancematrix)
}

.randomizeDistance <- function(distmatrix)
{
  distmatrix <- as.matrix(distmatrix)
  ord <- sample(rownames(distmatrix))
  return(as.dist(distmatrix[ord, ord]))
}
