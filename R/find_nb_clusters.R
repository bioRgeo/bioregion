#' @note
#' If multiple evaluation metrics are requested, only the first one will
#' be used to determine the optimal number of clusters.
#' @references
#' \insertRef{Holt2013}{bioRgeo}
#' \insertRef{Kreft2010}{bioRgeo}
find_nb_clusters_tree <- function(
  tree,
  dist = NULL,
  k_min = 2,
  k_max = "number of sites",
  eval_metric = "pc_distance",
  criterion = "gap", # gap or cutoff for now
  metric_cutoffs = c(.5, .75, .9, .95, .99, .999),
  gap_quantile = .99,
  plot = TRUE
)
{
  if(inherits(tree, "bioRgeo.hierar.tree"))
  {
    tree_object <- tree$final.tree
    dist_object <- tree$dist.matrix
    # Update args
    # tree$args[c("n_clust", "cut_height", "find_h", "h_max", "h_min")] <-
      # list(n_clust, cut_height, find_h, h_max, h_min)
  } else if (inherits(tree, "hclust"))
  {
    if(is.null(dist))
    {
      stop("If tree is a hclust object, then you must provide the original distance matrix in dist_matrix")
    } else if (!inherits(dist, "dist"))
    {
      stop("dist is not a distance matrix")
    }
    tree_object <- tree
    dist_object <- dist
  } else
  {
    stop("This function is designed to work either on bioRgeo.hierar.tree (output from clustering_hierarchical()) or hclust objects.")
  }

  nb_sites <- attr(dist_object, "Size")

  if(k_max == "number of sites")
  {
    k_max <- nb_sites
  } else if(!(k_max %% 1 == 0)) # integer testing ain't easy in R
  {
    stop("k_max must be an integer determining the number of clusters.")
  }

  # Labels are not in the same order in the tree and in the distance matrix.
  # tree_object$labels == attr(dist_object, "Labels")
  # cbind(tree_object$labels,
  #       attr(dist_object, "Labels"))

  # Create a df with sites in rows and clusters in columns
  clusters <- data.frame(matrix(nrow = nb_sites,
                                ncol = length(k_min:k_max) + 1,
                                dimnames = list(attr(dist_object, "Labels"),
                                                c("site", paste0("k_", 1:length(k_min:k_max))))))

  clusters$site <- attr(dist_object, "Labels")

  dist_mat <- as.matrix(dist_object)
  dist_sum_total <- sum(dist_mat) # Calculation for metric "pc_distance"

  # Prepare evaluation data.frame
  evaluation_df <- data.frame(matrix(nrow = length(k_min:k_max),
                                     ncol = 1 + length(eval_metric),
                                     dimnames = list(k_min:k_max,
                                                     c("n_clusters", eval_metric))))

  evaluation_df$n_clusters <- k_min:k_max

  cur_col <- 2 # Start at the second column and proceed through all columns
  message(paste0("Exploring all clustering results on the tree between ", k_min,
                 " and ", k_max, " groups, this may take a while..."))

  for(k in k_min:k_max)
  {

    # Cut the tree for given k
    cls <- dendextend::cutree(tree_object, k = k)
    # Make sure clusters are assigned to their respective sites
    clusters[, cur_col] <- cls[match(rownames(clusters),
                                     names(cls))]


    if(eval_metric == "pc_distance")
    {
      # Compute total distance for current number of clusters
      # Create a distance matrix with only distances between clusters - not within
      # clusters
      cur_dist_mat <- dist_mat
      for(cur_cls in unique(clusters[, cur_col]))
      {
        # Sites for current cluster
        cur_sites <- clusters[which(clusters[, cur_col] == cur_cls), "site"]
        # Distances between sites of current cluster are set to 0
        cur_dist_mat[cur_sites, cur_sites] <- 0
      }
      # Calculate sum for the current cluster number
      evaluation_df$pc_distance[which(evaluation_df$n_clusters == k)] <- sum(cur_dist_mat) / dist_sum_total
    }

    # Move on to next k
    cur_col <- cur_col + 1
  }

  message(paste0("Clustering explorations finished, finding the optimal number(s) of clusters..."))

  if(criterion == "gap")
  {
    message(" - Gap method")
    # Compute difference between each consecutive nb of clusters
    diffs <- evaluation_df[2:nrow(evaluation_df), eval_metric[1]] -
      evaluation_df[1:(nrow(evaluation_df) - 1), eval_metric[1]]

    qt <- quantile(diffs,
                   gap_quantile)

    evaluation_df$breaks <- FALSE
    evaluation_df$breaks[(which(diffs > qt) + 1)] <- TRUE
  } else if(criterion == "cutoff")
  {
    message(" - Cutoff method")

    evaluation_df$breaks <- FALSE

    hits <- sapply(metric_cutoffs,
           function(x, metric)
             {
             which(metric - x > 0)[1]
           }, metric = evaluation_df[ eval_metric])
    if(all(is.na(hits)))
    {
      stop("No number of cluster satisfied the requested metric_cutoffs")
    } else if(any(is.na(hits)))
    {
      warning("Could not find suitable number of clusters for cutoff(s): ", metric_cutoffs[which(is.na(hits))])
    } else
    evaluation_df$breaks[hits] <- TRUE
  }

  optimal_nb_clusters <- evaluation_df$n_clusters[evaluation_df$breaks]


  if(plot)
  {
    ggplot2::ggplot(evaluation_df, ggplot2::aes_string(x = "n_clusters", y = eval_metric[1])) +
      ggplot2::geom_line(col = "darkgrey") +
      ggplot2::geom_hline(yintercept = evaluation_df[evaluation_df$breaks, eval_metric[1]],
                          linetype = 2) +
      ggplot2::theme_bw()
  }
  outputs <- list(args = list(k_min = k_min,
                              k_max = k_max,
                              eval_metric = eval_metric,
                              criterion = criterion,
                              metric_cutoffs = metric_cutoffs,
                              gap_quantile = gap_quantile),
                  evaluation_df = evaluation_df,
                  optimal_nb_clusters = optimal_nb_clusters)

  class(outputs) <- append("bioRgeo.cluster.optimisation", class(outputs))
  return(outputs)
}

