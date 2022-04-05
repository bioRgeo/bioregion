#' Find an optimal number of clusters in a hierachical tree based on distances
#' or beta-diversity indices.
#'
#' This function aims at finding an optimal number of clusters in a hierarchical
#' tree on the basis of the investigation of how an evaluation metric changes as
#' a function of the number of clusters, and a criterion to find optimal
#' number(s) of clusters upon this relationship between evaluation metric &
#' number of clusters.
#'
#' @param cluster_object tree a \code{bioRgeo.hierar.tree} or a \code{hclust} object
#' @param eval_metric character string or vector of character strings indicating
#'  metric(s) to be calculated to
#' investigate the effect of different number of clusters. Available options:
#' \code{"pc_distance"}, \code{"anosim"}, \code{"avg_endemism"},
#' \code{"tot_endemism"}
#' @param distances a \code{dist} object or a \code{bioRgeo.distance} object (output
#' from \code{\link{similarity_to_distance}}). Necessary if \code{eval_metric}
#' includes \code{pc_distance} and \code{tree} is not a
#' \code{bioRgeo.hierar.tree} object
#' @param distance_index a character string indicating the distance (beta-diversity)
#' index to be used in case \code{dist} is a \code{data.frame} with multiple
#' distance indices
#' @param tree_k_min an integer indicating the minimum number of clusters to be
#' explored. Only useful if
#' \code{cluster_object} is based on \code{\link{clustering_hierarchical}}.
#' @param tree_k_max an integer indicating the maximum number of clusters to be
#' explored, or "number of sites" to use the number of sites as the maximum.
#' Only useful if
#' \code{cluster_object} is based on \code{\link{clustering_hierarchical}}.
#' @param tree_force_repartitioning a boolean indicating if the function should
#' re-partition clusters based on \code{tree_k_min} and \code{tree_k_max}. Only useful if
#' \code{cluster_object} is based on \code{\link{clustering_hierarchical}}.
#' @param partition_optimisation a boolean specifying if the function should
#' find an optimal number of clusters, based on the chosen \code{criterion}
#' @param criterion character string indicating the criterion to be used to
#' identify optimal number(s) of clusters. Available methods currently include
#' \code{"step"}, \code{"cutoff"}, \code{"elbow"}, \code{"mars"}, \code{"min"} or
#' \code{"max"}.
#' @param step_quantile if \code{criterion = "step"}, specify here the quantile
#' of differences between two consecutive k to be used as the cutoff to identify
#' the most important steps in \code{eval_metric}
#' @param step_levels if \code{criterion = "step"}, specify here the number of
#' steps to keep as cutoffs.
#' @param metric_cutoffs if \code{criterion = "cutoff"}, specify here the
#' cutoffs of \code{eval_metric} at which the number of clusters should be
#' extracted
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
#' between \code{tree_k_min} and \code{tree_k_max} number of clusters are explored on
#' the input object  \code{tree} by cutting the tree for each number of groups
#' \code{k} between
#' \code{tree_k_min} and \code{tree_k_max}. Second, for each clusterisation, the function
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
#' \item{\code{increasing_step} or \code{decreasing_step}:
#' This method consists in identifying clusters at the most important
#' changes, or steps, in the evaluation metric. The objective can be to either
#' look for largest increases (\code{increasing_step}) or largest decreases
#' \code{decreasing_step}. Steps are calculated based on the pairwise
#' differences between partitions.
#' Therefore, this is relative
#' to the distribution of changes in the evaluation metric over the tested
#' \code{k}. Specify \code{step_quantile} as the quantile cutoff above which
#' steps will be selected as most important (by default, 0.99, i.e. the
#' largest 1\% steps will be selected).Alternatively, you can also choose to
#' specify the number of top steps to keep, e.g. to keep the largest three
#' steps, specify \code{step_level = 3}
#' Basically this method will emphasize the
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
#' It is based on a fast method finding the maximum distance between the curve
#' and a straight line linking the minimum and maximum number of points.
#' The code we use here is based on code written by Esben Eickhardt available
#' here \url{https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve/42810075#42810075}.
#' The code has been modified to work on both increasing and decreasing
#' evaluation metrics.}
#' \item{\code{mars}:
#' This method consists in fitting a mars model on the evaluation curve, and
#' using it to identify all cutoffs at which there is no more increase in the
#' evaluation metric. In other words, this method will find cutoffs with the two
#' following conditions: (1) the evaluation metric was increasing before the
#' cutoff and (2) there is no more increase or the increase is slower after the
#' cutoff. This method uses \link[earth:earth]{earth::earth()}.}
#'
#' \item{\code{min} & \code{max}:
#' Picks the optimal partition(s) respectively at the minimum or maximum value
#' of the evaluation metric.}
#' }
#' @return
#' a \code{list} of class \code{bioRgeo.partition.metrics} with three elements:
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
#' simil <- similarity(vegemat, metric = "all")
#' distances <- similarity_to_distance(simil)
#'
#' # User-defined number of clusters
#' tree1 <- clustering_hierarchical(distances,
#'                                  n_clust = 5,
#'                                  index = "Simpson")
#' tree1
#'
#' a <- partition_metrics(tree1,
#'                   distances = distances,
#'                   eval_metric = c("tot_endemism",
#'                                   "avg_endemism",
#'                                   "pc_distance",
#'                                   "anosim"))
#'
#' a <- partition_metrics(tree1, tree_k_max = 50,
#'                   tree_force_repartitioning = TRUE,
#'                   distances = distances,
#'                   step_levels = 5)
#'  a <- partition_metrics(tree1, tree_k_max = 50,
#'                   eval_metric = c("tot_endemism",
#'                                   "avg_endemism",
#'                                   "pc_distance",
#'                                   "anosim"),
#'                   tree_force_repartitioning = TRUE,
#'                   distances = distances,
#'                   sp_site_table = vegemat)
#'
#' a <- partition_metrics(tree1,
#'                   eval_metric = c("tot_endemism",
#'                                   "avg_endemism",
#'                                   "pc_distance",
#'                                   "anosim"),
#'                   tree_k_max = 25,
#'                   tree_force_repartitioning = TRUE,
#'                   partition_optimisation = TRUE,
#'                   distances = distances,
#'                   sp_site_table = vegemat,
#'                   criterion = "decreasing_step",
#'                   step_levels = 5)
#'
#' partition_metrics(tree1,
#'                  distances = distances,
#'                  eval_metric = "pc_distance",
#'                  tree_k_max = 50,
#'                  tree_force_repartitioning = TRUE,
#'                  partition_optimisation = TRUE,
#'                  criterion = "elbow")
#'
#' partition_metrics(tree1,
#'                  distances = distances,
#'                  eval_metric = "pc_distance",
#'                  tree_k_max = 50,
#'                  tree_force_repartitioning = TRUE,
#'                  partition_optimisation = TRUE,
#'                  criterion = "increasing_step")
#'
#' partition_metrics(tree1,
#'                  distances = distances,
#'                  eval_metric = "pc_distance",
#'                  tree_k_max = 50,
#'                  tree_force_repartitioning = TRUE,
#'                  partition_optimisation = TRUE,
#'                  criterion = "mars")
#'
#' partition_metrics(tree1,
#'                  sp_site_table = vegemat,
#'                  eval_metric = "tot_endemism",
#'                  tree_k_max = 50,
#'                  tree_force_repartitioning = TRUE,
#'                  partition_optimisation = TRUE,
#'                  criterion = "max")
partition_metrics <- function(
  cluster_object,
  distances = NULL,
  distance_index = names(distances)[3],
  eval_metric = "pc_distance",
  tree_k_min = 2,
  tree_k_max = "number of sites",
  tree_force_repartitioning = FALSE,
  partition_optimisation = FALSE,
  criterion = "elbow", # step, elbow, cutoff, min, max and mars for now
  step_quantile = .99,
  step_levels = NULL,
  metric_cutoffs = c(.5, .75, .9, .95, .99, .999),
  sp_site_table = NULL,
  plot = TRUE,
  disable_progress = FALSE
)
{
  distance_based_metrics <- c("pc_distance",
                              "anosim")
  compo_based_metrics <- c("avg_endemism",
                           "tot_endemism")

  # 1. Does input object contains partitions already? ---------------
  if (inherits(cluster_object, "bioRgeo.clusters")) {
    if (inherits(cluster_object$clusters, "data.frame") & # Partitions already exist?
        !all(tree_force_repartitioning & # Or should we force repartition, if object comes from hierarchical clustering?
             cluster_object$name == "hierarchical_clustering")) {
      has.clusters <- TRUE
    } else {
      if (cluster_object$name == "hierarchical_clustering") {
        has.clusters <- FALSE
        tree_object <- cluster_object$algorithm$final.tree
      } else {
        stop("cluster_object does not have the expected type of 'clusters' slot")
      }
    }
  } else if (!inherits(cluster_object, "hclust")) {
    stop("This function is designed to work either on bioRgeo.clusters objects (outputs from clustering functions) or hclust objects.")
  } else {
    tree_object <- cluster_object
  }


  if (is.null(distances)) {
    has.distances <- FALSE
    if(any(eval_metric %in% distance_based_metrics))
    {
      warning(paste0("No distance matrix provided, so metrics ",
                     paste(eval_metric[which(eval_metric %in% distance_based_metrics)],
                           collapse = ", "),
                     " will not be computed"))
      eval_metric <- eval_metric[-which(eval_metric %in% distance_based_metrics)]
    }
  } else if (inherits(distances, "bioRgeo.pairwise.metric")) {
    if (attr(distances, "type") == "distance") {
      dist_object <- stats::as.dist(
        net_to_mat(distances[, c(
          1, 2,
          which(colnames(distances) == distance_index)
        )],
        weight = TRUE, squared = TRUE, symmetrical = TRUE
        )
      )
      has.distances <- TRUE
    } else
    {
      stop("distances must be an object containing distances from distance() or similarity_to_distance(), or an object of class dist")
    }
  } else if (!inherits(distances, "dist")) {
    stop("distances must be an object containing distances from distance() or similarity_to_distance(), or an object of class dist")
  }

  if (is.null(sp_site_table)) {
    has.contin <- FALSE
    if(any(eval_metric %in% compo_based_metrics))
    {
      warning(paste0("No species-site table provided, so metrics ",
                     paste(eval_metric[which(eval_metric %in% compo_based_metrics)],
                           collapse = ", "),
                     " will not be computed"))
      eval_metric <- eval_metric[-which(eval_metric %in% compo_based_metrics)]
    }
  } else {
    if(has.clusters)
    {
      if(any(!rownames(sp_site_table) %in% cluster_object$clusters$site))
      {
        stop("sp_site_table should be a matrix with sites in rows, species in columns.\nRow names should be site names, and should be identical to site names in the cluster object")
      }
    } else
    {
      if(any(!rownames(sp_site_table) %in% tree_object$labels))
      {
        stop("sp_site_table should be a matrix with sites in rows, species in columns.\nRow names should be site names, and should be identical to names in the tree")
      }
    }
    has.contin <- TRUE
  }

  if(!length(eval_metric))
  {
    stop("No evaluation metric can be computed because of missing arguments. Check arguments distances and sp_site_table")
  }

  if(partition_optimisation & length(criterion) > 2)
  {
    stop("Please provide only one method to select the optimal number of clusters, in argument 'criterion'")
  }


  # 2. Create partitions on hclust objects if need be --------------
  nb_sites <- cluster_object$inputs$nb_sites

  if(!has.clusters) {

    if(tree_k_max == "number of sites")
    {
      if("anosim" %in% eval_metric)
      {
        message("tree_k_max set to nb_sites - 1, such that anosim can be computed")
        tree_k_max <- nb_sites - 1
      } else
      {
        tree_k_max <- nb_sites
      }
    } else if(!(tree_k_max %% 1 == 0)) # integer testing ain't easy in R
    {
      stop("tree_k_max must be an integer determining the number of clusters.")
    }

    cluster_object <- suppressMessages(cut_tree(cluster_object,
                                                n_clust = tree_k_min:tree_k_max))

  }

  clusters <- cluster_object$clusters



  # 3. Calculate metrics ---------------------
  # Labels are not in the same order in the tree and in the distance matrix.
  # tree_object$labels == attr(dist_object, "Labels")
  # cbind(tree_object$labels,
  #       attr(dist_object, "Labels"))

  if(has.distances) {
    dist_mat <- as.matrix(dist_object)
    dist_sum_total <- sum(dist_mat) # Calculation for metric "pc_distance"

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
  }


  # Prepare evaluation data.frame
  evaluation_df <- data.frame(matrix(nrow = ncol(cluster_object$clusters) - 1,
                                     ncol = 2 + length(eval_metric),
                                     dimnames = list(colnames(cluster_object$clusters)[2:ncol(cluster_object$clusters)],
                                                     c("K", "n_clusters", eval_metric))))
  evaluation_df$K <- colnames(cluster_object$clusters)[2:ncol(cluster_object$clusters)]
  evaluation_df$n_clusters <- apply(cluster_object$clusters[, 2:(ncol(cluster_object$clusters)), drop = FALSE],
                                    2,
                                    function (x) length(unique(x)))

  evaluation_df <- evaluation_df[order(evaluation_df$n_clusters), ]



  cur_col <- 2 # Start at the second column and proceed through all columns
  message("Exploring metrics for all partitions, this may take a while...")
  for(cls_col in 2:ncol(clusters))
  {

    # cur_clusters <- clusters[, cur_col]

    if(has.distances){
      # The next line will create, for each element of the distance matrix, a
      # vector indicating whether if each distance is within or between clusters
      within_clusters <- clusters[rownames_dist,
                                  cls_col] == clusters[colnames_dist,
                                                       cls_col]

      if("pc_distance" %in% eval_metric)
      {
        # Compute total distance for current number of clusters
        # Create a distance matrix with only distances between clusters - not within
        # clusters
        cur_dist_object <- dist_object
        # Distances between sites of current cluster are set to 0
        cur_dist_object[within_clusters] <- 0

        # Calculate sum for the current cluster number
        evaluation_df$pc_distance[cls_col - 1] <- sum(cur_dist_object) / sum(dist_object)
      }

      if("anosim" %in% eval_metric)
      {
        dist_ranks <- rank(dist_object)

        denom <- nb_sites * (nb_sites - 1) / 4

        wb_average_rank <- tapply(dist_ranks, within_clusters, mean)

        # Analysis of similarities, formula from the vegan package
        evaluation_df$anosim[cls_col - 1] <- -diff(wb_average_rank) / denom
      }
    }


    if(has.contin)
    {
      cur_contin <- as.data.frame(sp_site_table)
      cur_contin$cluster <- clusters[match(rownames(cur_contin), clusters$site), cls_col]
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
        evaluation_df$avg_endemism[cls_col - 1] <- mean(pc_endemism_per_cluster)
      }
      # Total endemism
      if("tot_endemism" %in% eval_metric)
      {
        evaluation_df$tot_endemism[cls_col - 1] <- length(which(occ == 1)) / length(which(occ != 1))
      }
    }

    if(!disable_progress)
    {
      cat(".")
      if(cls_col == ncol(clusters))
      {
        cat("Complete\n")
      } else if((cls_col - 1) %% 50 == 0)
      {
        cat(paste0("k = ", cls_col - 1, " ~ ", round(100 * (cls_col - 1) / tree_k_max, 2), "% complete\n"))
      }
    }
  }


  # 4 find optimal number of clusters? ------------------
  if(partition_optimisation)
  {
    if(ncol(cluster_object$clusters) == 2)
    {
      message("Only one partition in the cluster table, skipping optimal number of cluster analysis")
      optimal_nb_clusters = NULL
    } else
    {
      message(paste0("Number of partitions: ", ncol(cluster_object$clusters) - 1), "\n")

      if(criterion %in% c("elbow",
                          "increasing_step",
                          "decreasing_step",
                          "mars")) {
        if(ncol(cluster_object$clusters) <= 5) {
          message(paste0("...Caveat: be cautious with the interpretation of metric analyses with such a low number of partitions"))
        }
      } else if(!(criterion %in% c("min", "max", "cutoff"))) {
        stop("criterion must be one of elbow, increasing_step, decreasing_step, min, max, cutoff or mars")
      }

      message(paste0("Searching for potential optimal number(s) of clusters based on the ",
                     criterion, " method"))

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

        elbow <- .elbow_finder(evaluation_df$n_clusters,
                               evaluation_df[, eval_metric[1]],
                               correct_decrease = TRUE)
        message(paste0("   * elbow found at ",
                       elbow[1], " clusters, rounding to ", round(elbow[1])))


        evaluation_df$optimal_nclust[which(evaluation_df$n_clusters == round(elbow[1]))] <- TRUE
      }
      if(criterion %in% c("increasing_step",
                          "decreasing_step"))
      {
        message(" - Step method")
        # Compute difference between each consecutive nb of clusters
        diffs <- evaluation_df[2:nrow(evaluation_df), eval_metric[1]] -
          evaluation_df[1:(nrow(evaluation_df) - 1), eval_metric[1]]
        # Should the first difference is considered to be an increment from 0?
        # diffs <- c(evaluation_df[1, eval_metric[1]],
        #                diffs)
        if(criterion == "decreasing_step"){
          diffs <- - diffs
        }

        if(!is.null(step_levels))
        {
          message(paste0("   * step_levels provided, selecting the ",
                         step_levels, " ",
                         "largest step(s) as cutoff(s)"))
          level_diffs <- diffs[order(diffs, decreasing = TRUE)][1:step_levels]
          evaluation_df$optimal_nclust[which(diffs %in% level_diffs) + 1] <- TRUE
        } else if(!is.null(step_quantile))
        {
          message(paste0("   * selecting the largest step(s) as cutoff(s) based on
                     the quantile ", step_quantile))
          qt <- stats::quantile(diffs,
                                step_quantile)
          evaluation_df$optimal_nclust[which(diffs > qt) + 1] <- TRUE
        }

      }
      if(criterion == "cutoff")
      {
        message(" - Cutoff method")


        hits <- sapply(metric_cutoffs,
                       function(x, metric)
                       {
                         which(metric - x > 0)[1]
                       }, metric = evaluation_df[, eval_metric])
        if(all(is.na(hits)))
        {
          stop("No number of cluster satisfied the requested metric_cutoffs")
        } else if(any(is.na(hits)))
        {
          warning("Could not find suitable number of clusters for cutoff(s): ", metric_cutoffs[which(is.na(hits))])
        } else
          evaluation_df$optimal_nclust[hits] <- TRUE
      }

      if(criterion == "max")
      {
        message(" - Max value method")

        evaluation_df$optimal_nclust[
          which(evaluation_df[, eval_metric] ==
                  max(evaluation_df[, eval_metric]))] <- TRUE
      }

      if(criterion == "min")
      {
        message(" - Min value method")

        evaluation_df$optimal_nclust[
          which(evaluation_df[, eval_metric] ==
                  max(evaluation_df[, eval_metric]))] <- TRUE
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
    }
  } else
  {
    optimal_nb_clusters = NULL
  }

  outputs <- list(args = list(eval_metric = eval_metric,
                              tree_k_min = ifelse(has.clusters, NA,
                                                  tree_k_min),
                              tree_k_max = ifelse(has.clusters, NA,
                                                  tree_k_max),
                              tree_force_repartitioning = tree_force_repartitioning,
                              partition_optimisation = partition_optimisation,
                              criterion = criterion,
                              metric_cutoffs = metric_cutoffs,
                              step_quantile = step_quantile,
                              step_levels = step_levels,
                              has_clusters = has.clusters),
                  evaluation_df = evaluation_df,
                  optimal_nb_clusters = optimal_nb_clusters)

  class(outputs) <- append("bioRgeo.partition.metrics", class(outputs))
  return(outputs)
}

