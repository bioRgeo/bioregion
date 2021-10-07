#' Cut a hierarchical tree
#'
#' This functions is designed to work on both \code{hclust} objects and
#' \code{bioRgeo.hierar.tree} objects. It cuts the tree for the chosen number(s)
#' of clusters or selected height(s). It also includes a procedure to
#' automatically return the height of cut for the chosen number(s) of clusters.
#'
#' @param tree a \code{bioRgeo.hierar.tree} or a \code{hclust} object
#' @param n_clust an integer indicating the number of clusters to be obtained
#' from the hierarchical tree. Should not be used at the same time as
#' \code{cut_height}
#' @param cut_height a numeric vector indicating the height(s) at which the tree
#' should be cut. Should not be used at the same time as \code{n_clust} or
#' \code{optim_method}
#' @param find_h a boolean indicating if the height of cut should be found for
#' the requested \code{n_clust}
#' @param h_max a numeric indicating the maximum possible tree height for
#' finding the height of cut when \code{find_h = TRUE}
#' @param h_min a numeric indicating the minimum possible height in the tree for
#' finding the height of cut when \code{find_h = TRUE}
#' @param dynamic_tree_cut a boolean indicating if the dynamic tree cut method
#' should be used, in which case \code{n_clust} & \code{cut_height} are ignored
#' @param dynamic_method a character vector indicating the method to be used
#' to dynamically cut the tree: either \code{"tree"} (clusters searched only
#' in the tree) or \code{"hybrid"} (clusters searched on both tree and distance
#' matrix)
#' @param dynamic_minClusterSize an integer indicating the minimum cluster size
#' to use in the dynamic tree cut method (see
#' \link[dynamicTreeCut:cutreeDynamic]{dynamicTreeCut::cutreeDynamic()})
#' @param dist_matrix only useful if \code{tree} is not a
#' \code{bioRgeo.hierar.tree} object and \code{dynamic_method = "hybrid"}.
#' Provide here the distance matrix used to build the \code{tree}
#' @param ... further arguments to be passed to
#' \link[dynamicTreeCut:cutreeDynamic]{dynamicTreeCut::cutreeDynamic()} to
#' customize the dynamic tree cut method.
#'
#' @details
#' The function can cut the tree with two main methods. First, it can cut
#' the entire tree at the same height (either specified by \code{cut_height} or
#' automatically defined for the chosen \code{n_clust}). Second, it can use
#' the dynamic tree cut method \insertRef{Langfelder2008}{bioRgeo}, in which
#' case clusters are detected with an adapative method based on the shape of
#' branches in the tree (thus cuts happen at multiple heights depending on
#' cluster positions in the tree).
#'
#' The dynamic tree cut method has two variants.
#' \itemize{
#' \item{The tree-based only variant
#' (\code{dynamic_method = "tree"}) is a top-down approach which relies only
#' on the tree and follows the order of clustered objects on it}
#' \item{The hybrid variant
#' (\code{dynamic_method = "hybrid"}) is a bottom-up approach which relies on
#' both the tree and the distance matrix to build clusters on the basis of
#' dissimilarity information among sites. This method is useful to detect
#' outlying members in each cluster.}
#' }
#'
#' @note
#' The argument \code{find_h} is ignored if \code{dynamic_tree_cut = TRUE},
#' because heights of cut cannot be estimated in this case.
#'
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}),
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) and
#' Boris Leroy (\email{leroy.boris@gmail.com})
#'
#' @return If \code{tree} is a \code{bioRgeo.hierar.tree}, then the same object
#' is returned with content updated (i.e., \code{args} and \code{clusters}). If
#' \code{tree} is a \code{hclust} object, then a \code{data.frame} containing
#' the clusters is returned.
#' @seealso \link{clustering_hierarchical}
#' @export
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
#' str(tree1)
#'
#' tree2 <- cut_tree(tree1,
#'                   cut_height = .05)
#' tree2
#' tree2$clusters
#'
#' tree3 <- cut_tree(tree1,
#'                   n_clust = c(3, 5, 10))
#' tree3
#' tree3$clusters
#'
#' tree4 <- cut_tree(tree1,
#'                   cut_height = c(.05, .1, .15, .2, .25))
#' tree4
#'
#' tree5 <- cut_tree(tree1,
#'                   n_clust = c(3, 5, 10),
#'                   find_h = FALSE)
#' tree5
#'
#'
#' hclust_tree <- tree2$final.tree
#' clusters_2 <- cut_tree(hclust_tree,
#'                        n_clust = 10)
cut_tree <- function(tree,
                     n_clust = NULL,
                     cut_height = NULL,
                     find_h = TRUE,
                     h_max = 1,
                     h_min = 0,
                     dynamic_tree_cut = FALSE,
                     dynamic_method = "tree",
                     dynamic_minClusterSize = 5,
                     dist_matrix = NULL,
                     ...)
{
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
    if(!is.null(cut_height))
    {
      stop("Please provide either n_clust or cut_height, but not both at the same time.")
    }
  }

  if(dynamic_tree_cut)
  {
    if(!is.null(n_clust))
    {
      message("The dynamic tree cut method was requested, argument n_clust will be ignored")
      n_clust <- NULL
    }
    if(!is.null(cut_height))
    {
      message("The dynamic tree cut method was requested, argument cut_height will be ignored")
      cut_height <- NULL
    }

    if(dynamic_method == "hybrid")
    {
      if(inherits(tree, "bioRgeo.hierar.tree"))
      {
        dist_matrix <- tree$dist.matrix
      } else if(is.null(dist_matrix))
      {
        stop("A distance matrix should be supplied to argument dist_matrix for dynamic_method = 'hybrid'")
      }
    }
  }

  arg_added <- list(...)
  if(inherits(tree, "bioRgeo.hierar.tree"))
  {
    cur.tree <- tree$final.tree
    # Update args
    tree$args[c("n_clust", "cut_height", "find_h", "h_max", "h_min")] <-
      list(n_clust, cut_height, find_h, h_max, h_min)

    if(dynamic_tree_cut)
    {
      tree$args[c("dynamic_tree_cut",
                  "dynamic_method",
                  "dynamic_minClusterSize")] <-
        list(dynamic_tree_cut, dynamic_method,  dynamic_minClusterSize)
      if(length(arg_added))
      {
        tree$args[names(arg_added)] <- arg_added
      }
    }
  } else if (inherits(tree, "hclust"))
  {
    cur.tree <- tree
  } else
  {
    stop("This function is designed to work either on bioRgeo.hierar.tree (output from clustering_hierarchical()) or hclust objects.")
  }


  output_cut_height <- NULL
  output_n_clust <- NULL
  if(dynamic_tree_cut)
  {
    clusters <- data.frame(site = tree$final.tree$labels,
                           cluster = dynamicTreeCut::cutreeDynamic(tree$final.tree,
                                                                   method = dynamic_method,
                                                                   minClusterSize = dynamic_minClusterSize,
                                                                   distM = dist_matrix,
                                                                   ...))
    # Set NAs for unassigned sites
    if(any(clusters$cluster == 0))
    {
      message("Some sites were not assigned to any cluster. They will have a NA in the cluster data.frame.")
      clusters$cluster[which(clusters$cluster == 0)] <- NA
    }
    output_n_clust <- length(unique(stats::na.omit(clusters$cluster)))
  } else if(!is.null(n_clust))
  {
    n_clust <- n_clust[order(n_clust)]
    if(find_h)
    {
      clusters <- data.frame(matrix(nrow = length(cur.tree$labels),
                                    ncol = length(n_clust) + 1,
                                    dimnames = list(cur.tree$labels,
                                                    c("site", paste0("k_", n_clust)))))
      clusters$site <- cur.tree$labels
      for(cur_n in n_clust)
      {
        message("Determining the cut height to reach ", cur_n, " groups...")
        k <- 0
        h1 <- h_max
        h0 <- h_min
        h <- h0 + (h1 - h0) / 2
        # Algorithm to quickly find the height of cut corresponding to the requested number of clusters
        while(k != cur_n & nchar(h) < 50 & h1 != h0)
        {
          h <- h0 + (h1 - h0) / 2
          cls <- dendextend::cutree(cur.tree, h = h)
          k <- max(cls)
          if(k < cur_n)
          {
            h1 <- h
          } else if (k > cur_n)
          {
            h0 <- h
          }
        }
        message(paste0("--> ", h))

        if(k != cur_n)
        {
          warning(paste0("The requested number of cluster could not be found for k = ", cur_n, ". Closest number found: ", k))
        }
        clusters[, paste0("k_", cur_n)] <- as.character(cls)
        output_cut_height <- c(output_cut_height, h)

        output_n_clust <- c(output_n_clust,
                            k)
      }
      names(output_cut_height) <- paste0("k_", n_clust)
    } else
    {
      cls <- dendextend::cutree(cur.tree, k = n_clust)
      clusters <- data.frame(rownames(cls),
                             cluster = cls)
      names(clusters) <- c("site", paste0("k_", n_clust))
      output_cut_height <- "unknown"
      output_n_clust <- sapply(n_clust,
                               function(k, cl){
                                 length(unique(cl[, paste0("k_", k)]))
                               }, cl = clusters)
    }


  } else if(!is.null(cut_height))
  {
    cut_height <- cut_height[order(cut_height, decreasing = TRUE)]
    cls <- dendextend::cutree(cur.tree,
                              h = cut_height)
    if(length(cut_height) == 1)
    {
      clusters <- data.frame(site = names(cls),
                             cluster = as.character(cls))
    } else
    {
      clusters <- data.frame(site = rownames(cls),
                             cluster = cls)
    }
    colnames(clusters) <- c("site", paste0("h_",  cut_height))
    output_n_clust <- sapply(cut_height,
                             function(h, cl){
                               length(unique(cl[, paste0("h_", h)]))
                             }, cl = clusters)
    names(output_n_clust) <- paste0("h_", cut_height)
    output_cut_height <- cut_height
  }

  if(inherits(tree, "bioRgeo.hierar.tree"))
  {
    cur.tree$args$cut_height <- cut_height
    tree$clusters <- clusters
    tree$output_n_clust <- output_n_clust
    tree$output_cut_height <- output_cut_height
    return(tree)
  } else if (inherits(tree, "hclust"))
  {
    return(clusters)
  }
}
