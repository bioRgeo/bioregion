#' Calculate metrics for one or several partitions
#' 
#' This function aims at calculating metrics for one or several partitions, 
#' usually on outputs from \code{netclu_}, \code{hclu_} or \code{nhclu_} 
#' functions. Metrics may require the users to provide either a similarity or
#' dissimilarity matrix, or to provide the initial species-site table.
#'
#' @param cluster_object tree a \code{bioRgeo.hierar.tree} or a \code{hclust} object
#' @param eval_metric character string or vector of character strings indicating
#'  metric(s) to be calculated to
#' investigate the effect of different number of clusters. Available options:
#' \code{"pc_distance"}, \code{"anosim"}, \code{"avg_endemism"},
#' \code{"tot_endemism"}
#' @param dissimilarity a \code{dist} object or a \code{bioRgeo.pairwise.metric} object (output
#' from \code{\link{similarity_to_dissimilarity}}). Necessary if \code{eval_metric}
#' includes \code{pc_distance} and \code{tree} is not a
#' \code{bioRgeo.hierar.tree} object
#' @param dissimilarity_index a character string indicating the dissimilarity (beta-diversity)
#' index to be used in case \code{dist} is a \code{data.frame} with multiple
#' dissimilarity indices
#' @param sp_site_table a \code{matrix} with sites in row and species in
#' columns. Should be provided if \code{eval_metric} includes
#' \code{"avg_endemism"} or \code{"tot_endemism"}
#' @param disable_progress a boolean to enable or disable the progress bar for
#' the exploration of clusters
#'
#' @details
#' \loadmathjax
#'
#' \bold{Evaluation metrics:}
#' \itemize{
#' \item{\code{pc_distance}: this metric is the method used by
#' \insertCite{Holt2013}{bioRgeo}. It is a ratio of the between-cluster sum of dissimilarity
#' (beta-diversity)
#' versus the total sum of dissimilarity (beta-diversity) for the full dissimilarity
#' matrix. In
#' other words, it is calculated on the basis of two elements. First, the total
#' sum of dissimilarity is calculated by summing the entire dissimilarity matrix
#' (\code{dist}). Second, the between-cluster sum of dissimilarity is calculated as
#' follows: for a given number of cluster, the dissimilarity is only
#' summed between clusters, not within clusters. To do that efficiently, all
#' pairs of sites within the same clusters have their dissimilarity set to zero in
#' the dissimilarity matrix, and then the dissimilarity matrix is summed. The
#' \code{pc_distance} ratio is obtained by dividing the between-cluster sum
#' of dissimilarity by the total sum of dissimilarity.}
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
#' @return
#' a \code{list} of class \code{bioRgeo.partition.metrics} with three elements:
#' \itemize{
#' \item{\code{args}: input arguments
#' }
#' \item{\code{evaluation_df}: the data.frame containing \code{eval_metric}
#' for all explored numbers of clusters
#' }
#' }
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
#' @seealso \link{hclu_hierarclust}
#' @examples
#' dissim <- dissimilarity(vegemat, metric = "all")
#'
#' # User-defined number of clusters
#' tree1 <- hclu_hierarclust(dissim,
#'                           n_clust = 2:20,
#'                           index = "Simpson")
#' tree1
#'
#' a <- partition_metrics(tree1,
#'                   dissimilarity = dissim,
#'                   sp_site_table = vegemat,
#'                   eval_metric = c("tot_endemism",
#'                                   "avg_endemism",
#'                                   "pc_distance",
#'                                   "anosim"))
#'

partition_metrics <- function(
  cluster_object,
  dissimilarity = NULL,
  dissimilarity_index = names(dissimilarity)[3],
  sp_site_table = NULL,
  eval_metric = "pc_distance",
  disable_progress = FALSE
)
{
  dissimilarity_based_metrics <- c("pc_distance",
                                   "anosim")
  compo_based_metrics <- c("avg_endemism",
                           "tot_endemism")

  # 1. Does input object contains partitions? ---------------
  if (inherits(cluster_object, "bioRgeo.clusters")) {
    if (inherits(cluster_object$clusters, "data.frame")) {
      has.clusters <- TRUE
    } else {
      if (cluster_object$name == "hierarchical_clustering") {
        stop("No clusters have been generated for your hierarchical tree, please extract clusters from the tree before using partition_metrics()
        See ?hclu_hierarclust or ?cut_tree")
      } else {
        stop("cluster_object does not have the expected type of 'clusters' slot")
      }
    }
  } else if (!inherits(cluster_object, "hclust")) {
    stop("This function is designed to work either on bioRgeo.clusters objects (outputs from clustering functions)")
  } 


  if (is.null(dissimilarity)) {
    has.dissimilarity <- FALSE
    if(any(eval_metric %in% dissimilarity_based_metrics))
    {
      warning(paste0("No dissimilarity matrix provided, so metrics ",
                     paste(eval_metric[which(eval_metric %in% dissimilarity_based_metrics)],
                           collapse = ", "),
                     " will not be computed"))
      eval_metric <- eval_metric[-which(eval_metric %in% dissimilarity_based_metrics)]
    }
  } else if (inherits(dissimilarity, "bioRgeo.pairwise.metric")) {
    if (attr(dissimilarity, "type") == "dissimilarity") {
      dist_object <- stats::as.dist(
        net_to_mat(dissimilarity[, c(colnames(dissimilarity)[1:2], dissimilarity_index)],
                   weight = TRUE, squared = TRUE, symmetrical = TRUE))
      has.dissimilarity <- TRUE
    } else
    {
      stop("dissimilarity must be an object containing dissimilarity indices from dissimilarity() or similarity_to_dissimilarity(), or an object of class dist")
    }
  } else if (!inherits(dissimilarity, "dist")) {
    stop("dissimilarity must be an object containing dissimilarity indices from dissimilarity() or similarity_to_dissimilarity(), or an object of class dist")
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
      if(any(!rownames(sp_site_table) %in% cluster_object$clusters$ID))
      {
        stop("sp_site_table should be a matrix with sites in rows, species in columns.\nRow names should be site names, and should be identical to site names in the cluster object")
      }
    } 
    has.contin <- TRUE
  }

  if(!length(eval_metric))
  {
    stop("No evaluation metric can be computed because of missing arguments. Check arguments dissimilarity and sp_site_table")
  }

  nb_sites <- cluster_object$inputs$nb_sites

  # 2. Calculate metrics ---------------------

  if(has.dissimilarity) {
    dist_mat <- as.matrix(dist_object)
    dist_sum_total <- sum(dist_mat) # Calculation for metric "pc_distance"

    # Create a vector of positions in the dissimilarity matrix indicating to what
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
  for(cls_col in 2:ncol(cluster_object$clusters))
  {

    # cur_clusters <- clusters[, cur_col]

    if(has.dissimilarity){
      # The next line will create, for each element of the dissimilarity matrix, a
      # vector indicating whether if each dissimilarity is within or between clusters
      within_clusters <- cluster_object$clusters[rownames_dist,
                                  cls_col] == cluster_object$clusters[colnames_dist,
                                                       cls_col]

      if("pc_distance" %in% eval_metric)
      {
        # Compute total dissimilarity for current number of clusters
        # Create a dissimilarity matrix with only dissimilarity between clusters - not within
        # clusters
        cur_dist_object <- dist_object
        # dissimilarity between sites of current cluster are set to 0
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
      cur_contin$cluster <- cluster_object$clusters[match(rownames(cur_contin), cluster_object$clusters$ID), cls_col]
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
      if(cls_col == ncol(cluster_object$clusters))
      {
        cat("Complete\n")
      } else if((cls_col - 1) %% 50 == 0)
      {
        cat(paste0("k = ", cls_col - 1, " ~ ", round(100 * (cls_col - 1) / ncol(cluster_object$clusters), 2), "% complete\n"))
      }
    }
  }


  outputs <- list(args = list(eval_metric = eval_metric),
                  evaluation_df = evaluation_df)

  class(outputs) <- append("bioRgeo.partition.metrics", class(outputs))
  return(outputs)
}

