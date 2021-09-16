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
#'
#' @details
#' Specify \code{n_clust} or \code{cut_height}, but not both at the same time.
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
                     h_min = 0)
{
  if(inherits(tree, "bioRgeo.hierar.tree"))
  {
    cur.tree <- tree$final.tree
    # Update args
    tree$args[c("n_clust", "cut_height", "find_h", "h_max", "h_min")] <-
      list(n_clust, cut_height, find_h, h_max, h_min)
  } else if (inherits(tree, "hclust"))
  {
    cur.tree <- tree
  } else
  {
    stop("This function is designed to work either on bioRgeo.hierar.tree (output from clustering_hierarchical()) or hclust objects.")
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
    if(!is.null(cut_height))
    {
      stop("Please provide either n_clust or cut_height, but not both at the same time.")
    }
  }

  output_cut_height <- NULL
  output_n_clust <- NULL
  if(!is.null(n_clust))
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

