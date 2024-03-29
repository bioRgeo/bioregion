#' Calculate metrics for one or several partitions
#' 
#' This function aims at calculating metrics for one or several partitions, 
#' usually on outputs from `netclu_`, `hclu_` or `nhclu_` functions. Metrics
#' may require the users to provide either a similarity or dissimilarity
#' matrix, or to provide the initial species-site table.
#'
#' @param cluster_object a `bioregion.clusters` object
#' 
#' @param eval_metric character string or vector of character strings indicating
#'  metric(s) to be calculated to investigate the effect of different number
#'  of clusters. Available options: `"pc_distance"`, `"anosim"`,
#'  `"avg_endemism"` and `"tot_endemism"`
#'  
#' @param dissimilarity a `dist` object or a `bioregion.pairwise.metric` object (output
#' from [similarity_to_dissimilarity()]). Necessary if `eval_metric`
#' includes `pc_distance` and `tree` is not a
#' `bioregion.hierar.tree` object
#' 
#' @param dissimilarity_index a character string indicating the dissimilarity
#' (beta-diversity) index to be used in case `dist` is a `data.frame` with
#' multiple dissimilarity indices
#' 
#' @param net the species-site network (i.e., bipartite network). Should be
#' provided if `eval_metric` includes `"avg_endemism"` or `"tot_endemism"`
#' 
#' @param site_col name or number for the column of site nodes (i.e. primary
#' nodes). Should be provided if `eval_metric` includes `"avg_endemism"` or
#' `"tot_endemism"`
#' 
#' @param species_col name or number for the column of species nodes (i.e.
#' feature nodes). Should be provided if `eval_metric` includes
#' `"avg_endemism"` or `"tot_endemism"`
#'
#' @details
#' \loadmathjax
#'
#' **Evaluation metrics:**
#' \itemize{
#' \item{`pc_distance`: this metric is the method used by
#' \insertCite{Holt2013}{bioregion}. It is a ratio of the between-cluster sum of
#' dissimilarity (beta-diversity) versus the total sum of dissimilarity
#' (beta-diversity) for the full dissimilarity matrix. In other words, it is
#' calculated on the basis of two elements. First, the total sum of
#' dissimilarity is calculated by summing the entire dissimilarity matrix
#' (`dist`). Second, the between-cluster sum of dissimilarity is calculated as
#' follows: for a given number of cluster, the dissimilarity is only summed
#' between clusters, not within clusters. To do that efficiently, all pairs of
#' sites within the same clusters have their dissimilarity set to zero in
#' the dissimilarity matrix, and then the dissimilarity matrix is summed. The
#' `pc_distance` ratio is obtained by dividing the between-cluster sum of
#' dissimilarity by the total sum of dissimilarity.}
#' 
#' \item{`anosim`: This metric is the statistic used in Analysis of
#' Similarities, as suggested in \insertCite{Castro-Insua2018}{bioregion} (see
#' [vegan::anosim()][vegan::anosim]). It compares the between-cluster
#' dissimilarities to the within-cluster dissimilarities. It is based based on
#' the difference of mean ranks between groups and within groups with the
#' following formula:
#' \mjeqn{R = (r_B - r_W)/(N (N-1) / 4)}{R = (r_B - r_W)/(N (N-1) / 4)},
#' where \mjeqn{r_B}{r_B} and \mjeqn{r_W}{r_W} are the average ranks
#' between and within clusters respectively, and \mjeqn{N}{N} is the total
#' number of sites.
#' Note that the function does not estimate the significance here, it only
#' computes the statistic - for significance testing see
#' [vegan::anosim()][vegan::anosim]}.
#' 
#' \item{`avg_endemism`: this metric is the average percentage of
#' endemism in clusters as
#' recommended by \insertCite{Kreft2010}{bioregion}. Calculated as follows:
#' \mjeqn{End_{mean} = \frac{\sum_{i=1}^K E_i / S_i}{K}}{Pc_endemism_mean = sum(Ei / Si) / K}
#'  where \mjeqn{E_i}{Ei} is the number of endemic species in cluster i,
#' \mjeqn{S_i}{Si} is the number of
#' species in cluster i, and K the maximum number of clusters.
#' }
#' 
#' \item{`tot_endemism`: this metric is the total endemism across all clusters,
#' as recommended by \insertCite{Kreft2010}{bioregion}. Calculated as follows:
#' \mjeqn{End_{tot} = \frac{E}{C}}{Endemism_total = E/C}
#'
#' where \mjeqn{E}{E} is total the number of endemics (i.e., species found in
#' only one cluster) and \mjeqn{C}{C} is the number of non-endemic species.
#' }
#' }
#'
#' @return
#' a `list` of class `bioregion.partition.metrics` with two to three elements:
#' \itemize{
#' \item{`args`: input arguments
#' }
#' \item{`evaluation_df`: the data.frame containing `eval_metric`
#' for all explored numbers of clusters
#' }
#' \item{`endemism_results`: if endemism calculations were requested, a list
#' with the endemism results for each partition
#' }
#' }
#' @import data.table 
#' 
#' @references
#' \insertRef{Castro-Insua2018}{bioregion}
#'
#' \insertRef{Ficetola2017}{bioregion}
#'
#' \insertRef{Holt2013}{bioregion}
#'
#' \insertRef{Kreft2010}{bioregion}
#'
#' \insertRef{Langfelder2008}{bioregion}
#' 
#' @author
#' Boris Leroy (\email{leroy.boris@gmail.com}),
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) and
#' Pierre Denelle (\email{pierre.denelle@gmail.com})
#' 
#' @seealso [compare_partitions]
#' @examples
#' comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
#' 20, 25)
#' rownames(comat) <- paste0("Site",1:20)
#' colnames(comat) <- paste0("Species",1:25)
#' 
#' comnet <- mat_to_net(comat)
#' 
#' dissim <- dissimilarity(comat, metric = "all")
#' 
#' # User-defined number of clusters
#' tree1 <- hclu_hierarclust(dissim, n_clust = 2:20, index = "Simpson")
#' tree1
#' 
#' a <- partition_metrics(tree1, dissimilarity = dissim, net = comnet,
#'                        site_col = "Node1", species_col = "Node2",
#'                        eval_metric = c("tot_endemism", "avg_endemism",
#'                                       "pc_distance", "anosim"))
#' a
#'
#' @importFrom stats as.dist
#' @importFrom data.table chmatch as.data.table
#'
#'@export

partition_metrics <- function(
    cluster_object, dissimilarity = NULL,
    dissimilarity_index = NULL, 
    net = NULL,
    site_col = 1, species_col = 2,
    eval_metric = c("pc_distance", "anosim", "avg_endemism", "tot_endemism")){
  dissimilarity_based_metrics <- c("pc_distance", "anosim")
  compo_based_metrics <- c("avg_endemism", "tot_endemism")
  
  # 1. Does input object contains partitions? ---------------------------------
  if (inherits(cluster_object, "bioregion.clusters")) {
    if (inherits(cluster_object$clusters, "data.frame")) {
      has.clusters <- TRUE # To remove? Does not seem relevant anymore
    } else {
      if (cluster_object$name == "hierarchical_clustering") {
        stop("No clusters have been generated for your hierarchical tree,
        please extract clusters from the tree before using partition_metrics()
        See ?hclu_hierarclust or ?cut_tree")
      } else {
        stop(
          "cluster_object does not have the expected type of 'clusters' slot")
      }
    }
  } else {
    stop("This function is designed to work  on bioregion.clusters objects
         (outputs from clustering functions)")
    # Add here the possibility to work on data.frame / matrices of clusters
    # directly
  } 
  
  if (is.null(dissimilarity)) {
    has.dissimilarity <- FALSE
    if(any(eval_metric %in% dissimilarity_based_metrics)){
      warning(paste0("No dissimilarity matrix provided, so metrics ",
                     paste(eval_metric[which(eval_metric %in%
                                               dissimilarity_based_metrics)],
                           collapse = ", "),
                     " will not be computed\n"))
      eval_metric <-
        eval_metric[-which(eval_metric %in% dissimilarity_based_metrics)]
    }
  } else if (inherits(dissimilarity, "bioregion.pairwise.metric")) {
    if (attr(dissimilarity, "type") == "dissimilarity") {
      if(is.null(dissimilarity_index)) {
        dissimilarity_index <- cluster_object$args$index
      } else if(!(dissimilarity_index %in% colnames(dissimilarity))) {
        stop("dissimilarity_index does not exist in the dissimilarity object.
             Did you misspecify the metric name?")
      }
      dist_object <- stats::as.dist(
        net_to_mat(dissimilarity[, c(colnames(dissimilarity)[1:2],
                                     dissimilarity_index)],
                   weight = TRUE, squared = TRUE, symmetrical = TRUE))
      has.dissimilarity <- TRUE
    } else{
      stop("dissimilarity must be an object containing dissimilarity indices
           from dissimilarity() or similarity_to_dissimilarity(), or an object
           of class dist")
    }
  } else if(!any(inherits(dissimilarity, "bioregion.pairwise.metric"),
                 inherits(dissimilarity, "dist"))){
    if(is.numeric(dissimilarity_index)){
      dissimilarity_index <- names(dissimilarity)[dissimilarity_index]
    }
    dist_object <- stats::as.dist(
      net_to_mat(dissimilarity[, c(colnames(dissimilarity)[1:2],
                                   dissimilarity_index)],
                 weight = TRUE, squared = TRUE, symmetrical = TRUE))
    has.dissimilarity <- TRUE
    
    if(!(dissimilarity_index %in% colnames(dissimilarity))){
      stop("dissimilarity is not a bioregion.pairwise.metric object, a
           dissimilarity matrix (class dist) or a data.frame with at least 3
           columns (site1, site2, and your dissimilarity index)")
    }
  }
  
  if (is.null(net)) {
    has.contin <- FALSE
    if(any(eval_metric %in% compo_based_metrics)){
      warning(paste0("No species-site network provided, so metrics ",
                     paste(eval_metric[which(eval_metric %in%
                                               compo_based_metrics)],
                           collapse = ", "), " will not be computed\n"))
      eval_metric <- eval_metric[-which(eval_metric %in% compo_based_metrics)]
    }
  } else {
    if(has.clusters){
      if(any(!(cluster_object$clusters$ID %in% c(net[, site_col],
                                                 net[, species_col])))) {
        stop("Some elements of the cluster table (column ID) cannot be found in
             the network")
      }
    } 
    if(is.numeric(site_col)){
      site_col <- names(net)[site_col]
    } 
    if(is.numeric(species_col)){
      species_col <- names(net)[species_col]
    } 
    
    # Next line is to use fast match with data.table, it needs characters
    net[, c(site_col, species_col)] <- lapply(net[, c(site_col, species_col)],
                                              as.character)
    has.contin <- TRUE
  }
  
  if(!length(eval_metric)){
    stop("No evaluation metric can be computed because of missing arguments.
         Check arguments dissimilarity and sp_site_table")
  }
  
  nb_sites <- cluster_object$inputs$nb_sites
  
  # 2. Calculate metrics ------------------------------------------------------
  
  if(has.dissimilarity) {
    dist_mat <- as.matrix(dist_object)
    dist_sum_total <- sum(dist_mat) # Calculation for metric "pc_distance"
    
    # Create a vector of positions in the dissimilarity matrix indicating to
    # what row/column each element corresponds
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
  evaluation_df <- data.frame(matrix(
    nrow = ncol(cluster_object$clusters) - 1,
    ncol = 2 + length(eval_metric),
    dimnames = list(colnames(
      cluster_object$clusters)[2:ncol(cluster_object$clusters)],
      c("K", "n_clusters", eval_metric))))
  
  evaluation_df$K <-
    colnames(cluster_object$clusters)[2:ncol(cluster_object$clusters)]
  
  evaluation_df$n_clusters <- apply(
    cluster_object$clusters[, 2:(ncol(cluster_object$clusters)), drop = FALSE],
    2,
    function (x) length(unique(x)))
  
  evaluation_df <- evaluation_df[order(evaluation_df$n_clusters), ]
  
  ## Check correspondence 
  # net_long <- data.frame(
  #   ID = c(unique(net[, site_col]),
  #          unique(net[, species_col])),
  #   nodetype = c(rep("site", length(unique(net[, site_col]))),
  #                rep("species", length(unique(net[, species_col])))))
  # net_long <- data.frame(
  #   net_long,
  #   cluster_object$clusters[match(net_long$ID,
  #                                 cluster_object$clusters$ID), -1])
  
  if(has.dissimilarity & any(c("pc_distance", "anosim") %in% eval_metric)){
    message("Computing similarity-based metrics...")
    # The next line will create, for each element of the dissimilarity matrix,
    # a vector indicating whether if each dissimilarity is within or between
    # clusters
    dissimilarity <- data.frame(
      dissimilarity,
      cluster_object$clusters[data.table::chmatch(dissimilarity$Site1,
                                                  cluster_object$clusters$ID), 
                              cluster_object$cluster_info$partition_name,
                              drop = FALSE],
      cluster_object$clusters[data.table::chmatch(dissimilarity$Site2,
                                                  cluster_object$clusters$ID), 
                              cluster_object$cluster_info$partition_name,
                              drop = FALSE])
    
    dissimilarity[, cluster_object$cluster_info$partition_name] <- 
      dissimilarity[, cluster_object$cluster_info$partition_name] ==
      dissimilarity[, paste0(cluster_object$cluster_info$partition_name, ".1")]
    
    dissimilarity <-
      dissimilarity[,
                    -which(colnames(dissimilarity) %in%
                             paste0(cluster_object$cluster_info$partition_name,
                                    ".1"))]
    
    if("pc_distance" %in% eval_metric) {
      evaluation_df$pc_distance <-
        vapply(cluster_object$cluster_info$partition_name,
               FUN = function(x, dist., index.) {
                 sum(dist.[!dist.[, x], index.]) / sum(dist.[, index.])
               },
               FUN.VALUE = numeric(1),
               dist. = dissimilarity, index. = dissimilarity_index)
      
      message("  - pc_distance OK")
    }
    
    if("anosim" %in% eval_metric){
      dissimilarity$ranks <- rank(dissimilarity[, dissimilarity_index])
      denom <- nb_sites * (nb_sites - 1) / 4
      
      # Fast calculation of the anosim for all clusters
      evaluation_df$anosim <- vapply(
        cluster_object$cluster_info$partition_name,
        FUN = function(x, dist., denom.) {
          # Testing if there is only one cluster
          if(all(dist.[, x])){
            NA # If only one cluster we cannot calculate anosim
          } else {
            -diff(tapply(dist.$ranks,
                         dist.[, x],
                         mean)) / denom.
          }
        },
        FUN.VALUE = numeric(1),
        dist. = dissimilarity, denom. = denom)
      message("  - anosim OK")
    }
  }
  
  if(has.contin & any(c("avg_endemism", "tot_endemism") %in% eval_metric)){
    message("Computing composition-based metrics...")
    net <- data.frame(
      net, 
      cluster_object$clusters[data.table::chmatch(net[, site_col],
                                                  cluster_object$clusters$ID),
                              -1])
    
    # Correcting column names when there is only one clustering
    if("cluster_object.clusters.data.table..chmatch.net...site_col..." %in%
       colnames(net)){
      colnames(net)[colnames(net) ==
                      "cluster_object.clusters.data.table..chmatch.net...site_col..."] <-
        colnames(cluster_object$clusters)[2]
    }

    # Visible binding for global variable
    N <- endemism <- end_richness <- pc_endemism <- NULL 
    
    # Fast calculation of endemism per cluster
    endemism_results <- lapply(
      cluster_object$cluster_info$partition_name,
      function(x, network., species_col.) {
        # Create species per cluster network in data.table format for faster
        # calculations
        species_cluster <- data.table::as.data.table(
          unique(network.[, c(species_col., x)]))
        # Calculate richness per cluster
        rich_clusters <- species_cluster[, .N, by = x]
        # Calculate species occurrence across clusters
        occ_sp <- species_cluster[, .N, by = species_col.]
        # Add new column with endemism status
        occ_sp[, "endemism" := N == 1] 
        # Then add endemism status to species_cluster table
        species_cluster$endemism <-
          occ_sp$endemism[data.table::chmatch(species_cluster[[species_col.]],
                                              occ_sp[[species_col.]])] 
        # Calculate richness of endemics per cluster
        end_clusters <- species_cluster[endemism == TRUE, .N, by = x]
        # Merge total & endemism richness tables
        rich_clusters[, end_richness := end_clusters[
          data.table::chmatch(rich_clusters[[x]],
                              end_clusters[[x]]), N]]
        # Replace NAs (i.e. no endemics) by zeros
        rich_clusters[is.na(rich_clusters)] <- 0 
        # Calculate percentage of endemism
        rich_clusters[, pc_endemism := end_richness/N] 
        return(rich_clusters)
      }, network. = net, species_col. = species_col
    )
    names(endemism_results) <- cluster_object$cluster_info$partition_name
    
    # Average endemism per cluster
    if("avg_endemism" %in% eval_metric){
      evaluation_df$avg_endemism <- vapply(
        cluster_object$cluster_info$partition_name,
        FUN = function(x, end_list) {
          mean(end_list[[x]]$pc_endemism)
        },
        FUN.VALUE = numeric(1),
        end_list = endemism_results
      )
      message("  - avg_endemism OK")
    }
    # Total endemism
    if("tot_endemism" %in% eval_metric){
      nb_sp <- length(unique(net[, species_col])) 
      evaluation_df$tot_endemism <- vapply(
        cluster_object$cluster_info$partition_name,
        FUN = function(x, end_list) {
          sum(end_list[[x]]$end_richness)
        },
        FUN.VALUE = numeric(1),
        end_list = endemism_results
      ) / nb_sp
      message("  - tot_endemism OK")
    }
  }
  
  outputs <- list(args = list(eval_metric = eval_metric),
                  evaluation_df = evaluation_df)
  
  if(has.contin & any(c("avg_endemism", "tot_endemism") %in% eval_metric)){
    outputs$endemism_results <- endemism_results
  }
  
  class(outputs) <- append("bioregion.partition.metrics", class(outputs))
  return(outputs)
}
