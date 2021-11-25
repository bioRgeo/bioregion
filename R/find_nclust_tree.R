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
#' @param eval_metric character string or vector of character strings indicating
#'  metric(s) to be calculated to
#' investigate the effect of different number of clusters. Available options:
#' \code{"pc_distance"}, \code{"anosim"}, \code{"avg_endemism"},
#' \code{"tot_endemism"}
#' Note that if you request several evaluation
#' metrics, they will all be computed, but only the first will be used to
#' find the optimal number(s) of clusters (order \code{eval_metric} accordingly)
#' @param criterion character string indicating the criterion to be used to
#' identify optimal number(s) of clusters. Available methods currently include
#' \code{"step"}, \code{"cutoff"}, \code{"elbow"} or \code{"mars"}.
#' @param step_quantile if \code{criterion = "step"}, specify here the quantile
#' of differences between two consecutive k to be used as the cutoff to identify
#' the most important steps in \code{eval_metric}
#' @param step_levels if \code{criterion = "step"}, specify here the number of
#' steps to keep as cutoffs.
#' @param metric_cutoffs if \code{criterion = "cutoff"}, specify here the
#' cutoffs of \code{eval_metric} at which the number of clusters should be
#' extracted
#' @param dist a \code{dist} object or a \code{bioRgeo.distance} object (output
#' from \code{\link{similarity_to_distance}}). Necessary if \code{eval_metric}
#' includes \code{pc_distance} and \code{tree} is not a
#' \code{bioRgeo.hierar.tree} object
#' @param dist_index a character string indicating the distance (beta-diversity)
#' index to be used in case \code{dist} is a \code{data.frame} with multiple
#' distance indices
#' @param sp_site_table a \code{matrix} with sites in row and species in
#' columns. Should be provided if \code{eval_metric} includes
#' \code{"avg_endemism"} or \code{"tot_endemism"}
#' @param plot a boolean indicating if a plot of the first \code{eval_metric}
#' should be drawn with the identified optimal numbers of cutoffs
#' @param disable_progress a boolean to enable or disable the progress bar for
#' the exploration of clusters
#'
#' @details
#' \loadmathjax
#'
#' This function proceeds in three steps. First, the range of clusterisations
#' between \code{k_min} and \code{k_max} number of clusters are explored on
#' the input object  \code{tree} by cutting the tree for each number of groups
#' \code{k} between
#' \code{k_min} and \code{k_max}. Second, for each clusterisation, the function
#' calculates one or several evaluation metric(s) and stores it. Third, the
#' relation ship evaluation metric ~ number of clusters is explored by the
#' function, and a criterion is applied on the first
#' evaluation metric to return an optimal number of clusters.
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
#' \item{\code{anosim}: This metric is the statistic used in Analysis of
#' Similarities, as suggested in \insertCite{Castro-Insua2018}{bioRgeo} (see
#' \link[vegan:anosim]{vegan::anosim()}). It
#' compares the between-cluster dissimilarities to the within-cluster
#' dissimilarities. It is based based on the difference of mean ranks between
#' groups and within groups with the following formula:
#' \mjeqn{R = (r_B - r_W)/(N (N-1) / 4)}{R = (r_B - r_W)/(N (N-1) / 4)},
#' where \mjeqn{r_B}{r_B} and \mjeqn{r_W}{r_W} are the average ranks
#' between and within clusters respectively, and \mjeqn{N}{N} is the total
#'  number of sites.
#' Note that the function does not estimate the significance here, it only
#' computes the statistic - for significance testing see
#' \link[vegan:anosim]{vegan::anosim()}}.
#' \item{\code{avg_endemism}: this metric is the average percentage of
#' endemism in clusters as
#' recommended by \insertCite{Kreft2010}{bioRgeo}. Calculated as follows:
#' \mjeqn{End_{mean} = \frac{\sum_{i=1}^K E_i / S_i}{K}}{Pc_endemism_mean = sum(Ei / Si) / K}
#'
#'  where \mjeqn{E_i}{Ei} is the number of endemic species in cluster i,
#' \mjeqn{S_i}{Si} is the number of
#' species in cluster i, and K the maximum number of clusters.
#' }
#' \item{\code{tot_endemism}: this metric is the total endemism in the
#' endemism in each cluster as
#' recommended by \insertCite{Kreft2010}{bioRgeo}. Calculated as follows:
#' \mjeqn{End_{tot} = \frac{E}{C}}{Endemism_total = E/C}
#'
#' where \mjeqn{E}{E} is total the number of endemics (i.e., species found in
#' only one cluster) and \mjeqn{C}{C} is number of non-endemics.
#' }
#' }
#'
#' \bold{Please read the note section about the following criteria.}
#'
#' Here we implemented a set of criteria commonly found in the literature or
#' recommended in the bioregionalisation literature. We advocate to move
#' beyond the "Search one optimal number of clusters" paradigm, and consider
#' investigating "multiple optimal numbers of clusters". Using only one cut on
#' a complex tree can overly simplify its information, and ignores the
#' hierarchical nature of the tree. Using multiple cuts likely avoids this
#' oversimplification bias and conveys more information about the hierarchical
#' nature of the tree. See, for example, the reanalysis of Holt et al. (2013)
#' by \insertCite{Ficetola2017}{bioRgeo}, where they used deep, intermediate
#' and shallow cuts. Following this rationale, several of the criteria
#' implemented here can/will return multiple "optimal" cuts.
#'
#' \bold{Criteria to find optimal number(s) of clusters}
#' \itemize{
#' \item{\code{step}:
#' This method consists in identifying clusters at the most important
#' increments, or steps, in the evaluation metric. Therefore, this is relative
#' to the distribution of increments in the evaluation metric over the tested
#' \code{k}. Specify \code{step_quantile} as the quantile cutoff above which
#' increments will be selected as most important (by default, 0.99, i.e. the
#' top 1\% increments will be selected). Basically this method will emphasize the
#' most important changes in the evaluation metric as a first approximation of
#' where important cuts can be chosen.
#' }
#' \item{\code{cutoffs}:
#' This method consists in specifying the cutoff value(s) in the evaluation
#' metric from which the number(s) of clusters should be derived. This is the
#' method used by \insertCite{Holt2013}{bioRgeo}. Note, however, that the
#' cut-offs suggested by Holt et al. (0.9, 0.95, 0.99, 0.999) may be only relevant
#' at very large spatial scales, and lower cut-offs should be considered at
#' finer spatial scales.
#' }
#' \item{\code{elbow}:
#' This method consists in finding one elbow in the evaluation metric curve, as
#' is commonly done in clustering analyses. The idea is to approximate the
#' number of clusters at which the evaluation metric no longer increments.
#' It is based on a fast numerical approximation of the elbow.
#' The code we use here is based on code written by Alan Tseng in the
#' KneeArrower R package \url{https://github.com/agentlans/KneeArrower},
#' and uses code from the
#' \code{signal} R package written in Octave by Paul Kienzle, modified by
#' Pascal Dupuis, and translated to R by Tom Short
#' \url{https://cran.r-project.org/package=signal}}
#' \item{\code{mars}:
#' This method consists in fitting a mars model on the evaluation curve, and
#' using it to identify all cutoffs at which there is no more increase in the
#' evaluation metric. In other words, this method will find cutoffs with the two
#' following conditions: (1) the evaluation metric was increasing before the
#' cutoff and (2) there is no more increase or the increase is slower after the
#' cutoff. This method uses \link[earth:earth]{earth::earth()}.}
#' }
#' @return
#' a \code{list} of class \code{bioRgeo.nclust.tree} with three elements:
#' \itemize{
#' \item{\code{args}: input arguments
#' }
#' \item{\code{evaluation_df}: the data.frame containing \code{eval_metric}
#' for all explored numbers of clusters
#' }
#' \item{\code{optimal_nb_clusters}: a vector containing the optimal number(s)
#' of cluster(s) according to the first computed \code{eval_metric} and the
#' chosen \code{criterion}
#' }}
#' @note Please note that finding the optimal number of clusters is a procedure
#' which normally requires decisions from the users, and as such can hardly be
#' fully automatized. Users are strongly advised to read the references
#' indicated below to look for guidance on how to choose their optimal number(s)
#' of clusters. Consider the "optimal" numbers of clusters returned by this
#' function as first approximation of the best numbers for your problematic.
#' @export
#' @references
#' \insertRef{Castro-Insua2018}{bioRgeo}
#'
#' \insertRef{Ficetola2017}{bioRgeo}
#'
#' \insertRef{Holt2013}{bioRgeo}
#'
#' \insertRef{Kreft2010}{bioRgeo}
#'
#' \insertRef{Langfelder2008}{bioRgeo}
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}),
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) and
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' @seealso \link{clustering_hierarchical}
#' @examples
#' simil <- spproject(vegemat, metric = "all")
#' distances <- similarity_to_distance(simil)
#'
#' # User-defined number of clusters
#' tree1 <- clustering_hierarchical(distances,
#'                                  n_clust = 5,
#'                                  index = "Simpson")
#' tree1
#'
#' find_nclust_tree(tree1)
#'
#' find_nclust_tree(tree1, step_levels = 5)
#'
#' find_nclust_tree(tree1, eval_metric = c("tot_endemism",
#'                                         "avg_endemism",
#'                                         "pc_distance",
#'                                         "anosim"),
#'                  criterion = "step",
#'                  sp_site_table = vegemat,
#'                  k_max = 50) -> b
#' find_nclust_tree(tree1, eval_metric = "pc_distance", criterion = "elbow")
#' find_nclust_tree(tree1, eval_metric = "pc_distance", criterion = "mars")
find_nclust_tree <- function(
  tree,
  k_min = 2,
  k_max = "number of sites",
  eval_metric = "pc_distance",
  criterion = "step", # step, elbow or cutoff for now
  step_quantile = .99,
  step_levels = NULL,
  metric_cutoffs = c(.5, .75, .9, .95, .99, .999),
  dist = NULL,
  dist_index = names(dist)[3],
  sp_site_table = NULL,
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
    if("anosim" %in% eval_metric)
    {
      message("k_max set to nb_sites - 1, such that anosim can be computed")
      k_max <- nb_sites - 1
    } else
    {
      k_max <- nb_sites
    }
  } else if(!(k_max %% 1 == 0)) # integer testing ain't easy in R
  {
    stop("k_max must be an integer determining the number of clusters.")
  }

  if(any(c("avg_endemism", "tot_endemism") %in% eval_metric))
  {
    if(is.null(sp_site_table))
    {
      stop("To calculate endemism metrics, you have to provide the species-site table in the sp_site_table argument.")
    }
    if(any(!rownames(sp_site_table) %in% tree_object$labels))
    {
      stop("sp_site_table should be a matrix with sites in rows, species in columns.\nRow names should be site names, and should be identical to names in the tree")
    }
  }

  # Labels are not in the same order in the tree and in the distance matrix.
  # tree_object$labels == attr(dist_object, "Labels")
  # cbind(tree_object$labels,
  #       attr(dist_object, "Labels"))

  # Create a df with sites in rows and clusters in columns
  clusters <- data.frame(matrix(nrow = nb_sites,
                                ncol = length(k_min:k_max) + 1,
                                dimnames = list(attr(dist_object, "Labels"),
                                                c("site", paste0("k_", (k_min:k_max))))))

  clusters$site <- attr(dist_object, "Labels")

  dist_mat <- as.matrix(dist_object)
  dist_sum_total <- sum(dist_mat) # Calculation for metric "pc_distance"

  # Prepare evaluation data.frame
  evaluation_df <- data.frame(matrix(nrow = length(k_min:k_max),
                                     ncol = 1 + length(eval_metric),
                                     dimnames = list(k_min:k_max,
                                                     c("n_clusters", eval_metric))))

  evaluation_df$n_clusters <- k_min:k_max


  # Create a vector of positions in the distance matrix indicating to what
  # row/column each element corresponds
  # Will be used to distinguish within vs. between clusters
  rownames_dist <- attr(dist_object,
                        "Labels")[as.vector(
                          stats::as.dist(row(matrix(nrow = nb_sites,
                                                    ncol = nb_sites))))]
  colnames_dist <- attr(dist_object,
                        "Labels")[as.vector(
                          stats::as.dist(col(matrix(nrow = nb_sites,
                                                    ncol = nb_sites))))]

  cur_col <- 2 # Start at the second column and proceed through all columns
  message(paste0("Exploring all clustering results on the tree between ", k_min,
                 " and ", k_max, " groups, this may take a while..."))

  for(k in k_min:k_max)
  {
    # Cut the tree for given k
    cls <- dendextend::cutree(tree_object, k = k)
    # Make sure clusters are assigned to their respective sites
    clusters[, paste0("k_", k)] <- cls[match(rownames(clusters),
                                             names(cls))]

    # cur_clusters <- clusters[, cur_col]

    # The next line will create, for each element of the distance matrix, a
    # vector indicating whether if each distance is within or between clusters
    within_clusters <- clusters[rownames_dist,
                                cur_col] == clusters[colnames_dist,
                                                     cur_col]

    if("pc_distance" %in% eval_metric)
    {
      # Compute total distance for current number of clusters
      # Create a distance matrix with only distances between clusters - not within
      # clusters
      cur_dist_object <- dist_object
      # Distances between sites of current cluster are set to 0
      cur_dist_object[within_clusters] <- 0

      # Calculate sum for the current cluster number
      evaluation_df$pc_distance[which(evaluation_df$n_clusters == k)] <- sum(cur_dist_object) / sum(dist_object)
    }

    if("anosim" %in% eval_metric)
    {

      dist_ranks <- rank(dist_object)

      denom <- nb_sites * (nb_sites - 1) / 4

      wb_average_rank <- tapply(dist_ranks, within_clusters, mean)

      # Analysis of similarities, formula from the vegan package
      evaluation_df$anosim[which(evaluation_df$n_clusters == k)] <- -diff(wb_average_rank) / denom
    }

    if(any(c("avg_endemism", "tot_endemism") %in% eval_metric))
    {
      cur_contin <- as.data.frame(sp_site_table)
      cur_contin$cluster <- clusters[match(rownames(cur_contin), clusters$site), paste0("k_", k)]
      cluster_contin <- stats::aggregate(. ~ cluster, cur_contin, sum)
      rownames(cluster_contin) <- cluster_contin[, 1]
      cluster_contin <- cluster_contin[, -1]
      cluster_contin[cluster_contin > 0] <- 1
      occ <- colSums(cluster_contin)
      cluster_contin_end <- cluster_contin[, which(occ == 1)]

      richness <- rowSums(cluster_contin)
      richness_end <- rowSums(cluster_contin_end)



      pc_endemism_per_cluster <- richness_end / richness

      # Average endemism per cluster
      if("avg_endemism" %in% eval_metric)
      {
        evaluation_df$avg_endemism[which(evaluation_df$n_clusters == k)] <- mean(pc_endemism_per_cluster)
      }
      # Total endemism
      if("tot_endemism" %in% eval_metric)
      {
        evaluation_df$tot_endemism[which(evaluation_df$n_clusters == k)] <- length(which(occ == 1)) / length(which(occ != 1))
      }
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

  evaluation_df$optimal_nclust <- FALSE
  if(criterion == "mars")
  {
    mars_model <- earth::earth(evaluation_df[, eval_metric[1]] ~ evaluation_df$n_clusters)
    mars_sum <- summary(mars_model)
    if(nrow(mars_sum$coefficients) > 1)
    {
      funs <- rownames(mars_sum$coefficients)[2:nrow(mars_sum$coefficients)]
      funs <- gsub(")", "", gsub("h(", "", funs, fixed = TRUE), fixed = TRUE)
      funs <- gsub(paste0("evaluation_df$n_clusters"), "", funs, fixed = TRUE)

      mars_hinges <- as.numeric(gsub("-", "", funs))

      mars_cutoffs <- data.frame(cutoff = unique(mars_hinges),
                                 before = NA, after = NA)
      mars_cutoffs <- mars_cutoffs[order(mars_cutoffs$cutoff), ]

      mars_preds <- data.frame(evaluation_df,
                               stats::predict(mars_model,
                                       evaluation_df$n_clusters))
      def_cut <- min(mars_preds$n_clusters)

      for (cur_cut in mars_cutoffs$cutoff)
      {

        mars_cutoffs$before[mars_cutoffs$cutoff == cur_cut] <-
          (mars_preds[which(mars_preds$n_clusters == cur_cut), ncol(mars_preds)] -
             mars_preds[which(mars_preds$n_clusters == def_cut), ncol(mars_preds)]) /
          (cur_cut - def_cut)

        def_cut <- cur_cut
      }
      if(nrow(mars_cutoffs) > 1)
      {
        mars_cutoffs$after[1:(nrow(mars_cutoffs) - 1)] <- mars_cutoffs$before[2:nrow(mars_cutoffs)]
      }
      mars_cutoffs$after[nrow(mars_cutoffs)] <-
        (mars_preds[which(mars_preds$n_clusters == max(mars_preds$n_clusters)), ncol(mars_preds)] -
           mars_preds[which(mars_preds$n_clusters == def_cut), ncol(mars_preds)]) /
        (max(mars_preds$n_clusters) - def_cut)
      mars_cutoffs$change_slope <- mars_cutoffs$after - mars_cutoffs$before

      mars_cutoffs$breaks <- ifelse(mars_cutoffs$before > 0 & mars_cutoffs$change_slope < 0, TRUE, FALSE)

      evaluation_df$optimal_nclust[which(evaluation_df$n_clusters %in% mars_cutoffs$cutoff[mars_cutoffs$breaks])] <- TRUE
    } else
    {
      stop("No cutoff point was found with the MARS method")
    }
  }

  if(criterion == "elbow")
  {
    message(" - Elbow method")

    elbow <- .findCutoff(evaluation_df$n_clusters,
                            evaluation_df[, eval_metric[1]])
    message(paste0("   * elbow found at ",
                   elbow$x, " clusters, rounding to ", round(elbow$x)))


    evaluation_df$optimal_nclust[which(evaluation_df$n_clusters == round(elbow$x))] <- TRUE
  } else if(criterion == "step")
  {
    message(" - Step method")
    # Compute difference between each consecutive nb of clusters
    diffs <- evaluation_df[2:nrow(evaluation_df), eval_metric[1]] -
      evaluation_df[1:(nrow(evaluation_df) - 1), eval_metric[1]]
    # First difference is considered to be an increment from 0
    diffs <- c(evaluation_df[1, eval_metric[1]],
               diffs)


    if(!is.null(step_levels))
    {
      message(paste0("   * step_levels provided, selecting the ",
                     step_levels, " highest step(s) as cutoff(s)"))
      level_diffs <- diffs[order(diffs, decreasing = TRUE)][1:step_levels]
      evaluation_df$optimal_nclust[which(diffs %in% level_diffs)] <- TRUE
    } else if(!is.null(step_quantile))
    {
      message(paste0("   * selecting the highest step(s) as cutoff(s) based on
                     the quantile ", step_quantile))
      qt <- stats::quantile(diffs,
                            step_quantile)
      evaluation_df$optimal_nclust[which(diffs > qt)] <- TRUE
    }

  } else if(criterion == "cutoff")
  {
    message(" - Cutoff method")


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
    evaluation_df$optimal_nclust[hits] <- TRUE
  }

  optimal_nb_clusters <- evaluation_df$n_clusters[evaluation_df$optimal_nclust]


  if(plot)
  {
    message("Plotting results...")
    p <- ggplot2::ggplot(evaluation_df, ggplot2::aes_string(x = "n_clusters", y = eval_metric[1])) +
      ggplot2::geom_line(col = "darkgrey") +
      # ggplot2::geom_hline(yintercept = evaluation_df[evaluation_df$optimal_nclust, eval_metric[1]],
      #                     linetype = 2) +
      ggplot2::geom_vline(xintercept = evaluation_df[evaluation_df$optimal_nclust, "n_clusters"],
                          linetype = 2) +
      ggplot2::theme_bw()
    if(criterion == "mars")
    {
      message("   (the red line is the prediction from MARS models)")
      p <- p + ggplot2::geom_line(data = mars_preds,
                                  ggplot2::aes_string(x = "n_clusters", y = mars_preds[, ncol(mars_preds)]),
                                  col = "red")
    }
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

  class(outputs) <- append("bioRgeo.nclust.tree", class(outputs))
  return(outputs)
}

