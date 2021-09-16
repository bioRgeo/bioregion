#' Find an optimal number of clusters in a hierachical tree based on distances
#' or beta-diversity indices.
#'
#' This function aims at finding an optimal number of clusters in a hierarchical
#' tree on the basis of the investigation of how an evaluation metric changes as
#' a function of the number of clusters, and a criterion to find optimal
#' number(s) of clusters upon this relationship between evaluation metric &
#' number of clusters.
#'
#' @param tree tree a \code{bioRgeo.hierar.tree} or a \code{hclust} object
#' @param k_min an integer indicating the minimum number of clusters to be
#' explored
#' @param k_max an integer indicating the maximum number of clusters to be
#' explored, or "number of sites" to use the number of sites as the maximum
#' @param eval_metric metric(s) to be calculated to investigate the effect of
#' different number of clusters. Note that if you request several evaluation
#' metrics, they will all be computed, but only the first will be used to
#' find the optimal number(s) of clusters (order \code{eval_metric} accordingly)
#' @param criterion character string indicating the criterion to be used to
#' identify optimal number(s) of clusters. Possible values include \code{"step"},
#' \code{cutoff}, ...
#' (to be completed)
#' @param step_quantile if \code{criterion = "step"}, specify here the quantile
#' of differences between two consecutive k to be used as the cutoff to identify
#' the most important steps in \code{eval_metric}
#' @param metric_cutoffs if \code{criterion = "cutoff"}, specify here the
#' cutoffs of \code{eval_metric} at which the number of clusters should be
#' extracted
#' @param dist a \code{dist} object or a \code{bioRgeo.distance} object (output
#' from \code{\link{similarity_to_distance}}). Necessary if \code{eval_metric}
#' includes \code{pc_distance}
#' @param dist_index a character string indicating the distance (beta-diversity)
#' index to be used in case \code{dist} is a \code{data.frame} with multiple
#' distance indices
#' @param plot a boolean indicating if a plot of the first \code{eval_metric}
#' should be drawn with the identified optimal numbers of cutoffs
#' @param disable_progress a boolean to enable or disable the progress bar for
#' the exploration of clusters
#'
#' @details
#' This function proceeds in three steps. First, the range of clusterisations
#' between \code{k_min} and \code{k_max} number of clusters are explored on
#'  \code{tree} by cutting the tree for each number of groups \code{k} between
#' \code{k_min} and \code{k_max}. Second, an evaluation metric
#' is calculated for each clusterisation. Third, a criterion is applied on the
#' evaluation metric to identify the optimal number of clusers.
#'
#' Note that multiple evaluation metrics can be requested (e.g., for inspection),
#' but only the first one will be used to apply the criterion for optimal
#' number of clusters.
#'
#' \bold{Evaluation metrics:}
#' \itemize{
#' \item{\code{pc_distance}: this metric is the method used by
#' \insertCite{Holt2013}{bioRgeo}. It is a ratio of the between-cluster sum of distances
#' (beta-diversity)
#' versus the total sum of distances (beta-diversity) for the full distance
#' matrix. In
#' other words, it is calculated on the basis of two elements. First, the total
#' sum of distances is calculated by summing the entire distance matrix
#' (\code{dist}). Second, the between-cluster sum of distances is calculated as
#' follows: for a given number of cluster, the distances are only
#' summed between clusters, not within clusters. To do that efficiently, all
#' pairs of sites within the same clusters have their distance set to zero in
#' the distance matrix, and then the distance matrix is summed. The
#' \code{pc_distance} ratio is obtained by dividing the between-cluster sum
#' of distances by the total sum of distances.}
#' \itemize{\code{pc_endemism}: this metric is the percentage of endemism as
#' recommended by \insertCite{Kreft2010}{bioRgeo}. To be detailed}
#' }
#'
#' \bold{Criteria to find optimal number(s) of clusters}
#' \itemize{
#' \item{\code{step}:
#' This methods consists in identifying clusters at the most important
#' increments, or steps, in the evaluation metric. Therefore, this is relative
#' to the distribution of increments in the evaluation metric over the tested
#' \code{k}. Specify \code{step_quantile} as the quantile cutoff above which
#' increments will be selected as most important (by default, 0.99, i.e. the
#' top 1\% increments will be selected).
#' }
#' \item{\code{cutoffs}:
#' This methods consists in specifying the cutoff value(s) in the evaluation
#' metric from which the number(s) of clusters should be derived.
#' }
#' }
#'
#' @references
#' \insertRef{Holt2013}{bioRgeo}
#'
#' \insertRef{Kreft2010}{bioRgeo}
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}),
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) and
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' @seealso \link{clustering_hierarchical}
find_nclust_tree <- function(
  tree,
  k_min = 2,
  k_max = "number of sites",
  eval_metric = "pc_distance",
  criterion = "step", # step or cutoff for now
  step_quantile = .99,
  metric_cutoffs = c(.5, .75, .9, .95, .99, .999),
  dist = NULL,
  dist_index = names(dist)[3],
  plot = TRUE,
  disable_progress = FALSE
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
    } else if (inherits(dist, "bioRgeo.distance") | dist_index %in% colnames(dist))
    {
      dist_object <- .dfToDist(dist, dist_index)
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

    if(!disable_progress)
    {
      cat(".")
      if((cur_col - 1) == length(k_min:k_max))
      {
        cat("Complete\n")
      } else if((cur_col - 1) %% 50 == 0)
      {
        cat(paste0("k = ", cur_col - 1, " ~ ", round(100 * (cur_col - 1) / k_max, 2), "% complete\n"))
      }
    }
    # Move on to next k
    cur_col <- cur_col + 1
  }

  message(paste0("Clustering explorations finished, finding the optimal number(s) of clusters..."))

  if(criterion == "step")
  {
    message(" - Step method")
    # Compute difference between each consecutive nb of clusters
    diffs <- evaluation_df[2:nrow(evaluation_df), eval_metric[1]] -
      evaluation_df[1:(nrow(evaluation_df) - 1), eval_metric[1]]

    qt <- stats::quantile(diffs,
                          step_quantile)

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
    message("Plotting results...")
    p <- ggplot2::ggplot(evaluation_df, ggplot2::aes_string(x = "n_clusters", y = eval_metric[1])) +
      ggplot2::geom_line(col = "darkgrey") +
      ggplot2::geom_hline(yintercept = evaluation_df[evaluation_df$breaks, eval_metric[1]],
                          linetype = 2) +
      ggplot2::theme_bw()
    print(p)
  }
  outputs <- list(args = list(k_min = k_min,
                              k_max = k_max,
                              eval_metric = eval_metric,
                              criterion = criterion,
                              metric_cutoffs = metric_cutoffs,
                              step_quantile = step_quantile),
                  evaluation_df = evaluation_df,
                  optimal_nb_clusters = optimal_nb_clusters)

  class(outputs) <- append("bioRgeo.cluster.optimisation", class(outputs))
  return(outputs)
}

