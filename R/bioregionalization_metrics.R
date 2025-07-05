#' Calculate metrics for one or several bioregionalizations
#' 
#' This function calculates metrics for one or several bioregionalizations, 
#' typically based on outputs from `netclu_`, `hclu_`, or `nhclu_` functions. 
#' Some metrics may require users to provide either a similarity or dissimilarity 
#' matrix, or the initial species-site table.
#'
#' @param bioregionalization A `bioregion.clusters` object.
#' 
#' @param eval_metric A `character` vector or a single `character` string 
#' indicating the metric(s) to be calculated to assess the effect of different 
#' numbers of clusters. Available options are `"pc_distance"`, `"anosim"`,
#' `"avg_endemism"`, or `"tot_endemism"`. If `"all"` is specified, all metrics
#'  will be calculated.
#'  
#' @param dissimilarity A `dist` object or a `bioregion.pairwise` 
#' object (output from [similarity_to_dissimilarity()]). Required if 
#' `eval_metric` includes `"pc_distance"` and `tree` is not a
#' `bioregion.hierar.tree` object.
#' 
#' @param dissimilarity_index A `character` string indicating the dissimilarity
#' (beta-diversity) index to use if dissimilarity is a `data.frame` with
#' multiple dissimilarity indices.
#' 
#' @param net The site-species network (i.e., bipartite network). Should be
#' provided as a `data.frame` if `eval_metric` includes `"avg_endemism"` or 
#' `"tot_endemism"`.
#' 
#' @param site_col The name or index of the column representing site nodes 
#' (i.e., primary nodes). Should be provided if `eval_metric` includes 
#' `"avg_endemism"` or `"tot_endemism"`.
#' 
#' @param species_col The name or index of the column representing species nodes 
#' (i.e., feature nodes). Should be provided if `eval_metric` includes
#' `"avg_endemism"` or `"tot_endemism"`.
#' 
#' @return A `list` of class `bioregion.bioregionalization.metrics` with two to three elements:
#' \itemize{
#' \item{`args`: Input arguments.}
#' \item{`evaluation_df`: A `data.frame` containing the `eval_metric`
#' values for all explored numbers of clusters.}
#' \item{`endemism_results`: If endemism calculations are requested, a list
#' with the endemism results for each bioregionalization.}
#' }
#'
#' @details
#' **Evaluation metrics:**
#' 
#' \itemize{
#' 
#' \item{`pc_distance`: This metric, as used by Holt et al. (2013), is the 
#' ratio of the between-cluster sum of dissimilarities (beta-diversity) to the 
#' total sum of dissimilarities for the full dissimilarity matrix. It is calculated 
#' in two steps: 
#' - Compute the total sum of dissimilarities by summing all elements of the 
#' dissimilarity matrix.
#' - Compute the between-cluster sum of dissimilarities by setting within-cluster 
#' dissimilarities to zero and summing the matrix. 
#' The `pc_distance` ratio is obtained by dividing the between-cluster sum of 
#' dissimilarities by the total sum of dissimilarities.}
#' 
#' \item{`anosim`: This metric is the statistic used in the Analysis of 
#' Similarities, as described in Castro-Insua et al. (2018). It compares 
#' between-cluster and within-cluster dissimilarities. The statistic is computed as: 
#' R = (r_B - r_W) / (N (N-1) / 4), 
#' where r_B and r_W are the average ranks of between-cluster and within-cluster 
#' dissimilarities, respectively, and N is the total number of sites. 
#' Note: This function does not estimate significance; for significance testing, 
#' use [vegan::anosim()][vegan::anosim].}
#' 
#' \item{`avg_endemism`: This metric is the average percentage of 
#' endemism in clusters, as recommended by Kreft & Jetz (2010). It is calculated as: 
#' End_mean = sum_i (E_i / S_i) / K, 
#' where E_i is the number of endemic species in cluster i, S_i is the number of 
#' species in cluster i, and K is the total number of clusters.}
#' 
#' \item{`tot_endemism`: This metric is the total endemism across all clusters, 
#' as recommended by Kreft & Jetz (2010). It is calculated as: 
#' End_tot = E / C, 
#' where E is the total number of endemic species (i.e., species found in only one 
#' cluster) and C is the number of non-endemic species.}
#' }
#' 
#' @references
#' Castro-Insua A, Gómez-Rodríguez C & Baselga A (2018) Dissimilarity measures 
#' affected by richness differences yield biased delimitations of biogeographic 
#' realms. \emph{Nature Communications} 9, 9-11.
#'
#' Holt BG, Lessard J, Borregaard MK, Fritz SA, Araújo MB, Dimitrov D, Fabre P, 
#' Graham CH, Graves GR, Jønsson Ka, Nogués-Bravo D, Wang Z, Whittaker RJ, 
#' Fjeldså J & Rahbek C (2013) An update of Wallace's zoogeographic regions of 
#' the world. \emph{Science} 339, 74-78.
#'
#' Kreft H & Jetz W (2010) A framework for delineating biogeographical regions
#' based on species distributions. \emph{Journal of Biogeography} 37, 2029-2053.
#' 
#' @seealso 
#' For more details illustrated with a practical example, 
#' see the vignette: 
#' \url{https://biorgeo.github.io/bioregion/articles/a4_1_hierarchical_clustering.html#optimaln}.
#' 
#' Associated functions: 
#' [compare_bioregionalizations] [find_optimal_n]
#' 
#' @author
#' Boris Leroy (\email{leroy.boris@gmail.com}) \cr
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) \cr
#' Pierre Denelle (\email{pierre.denelle@gmail.com})
#' 
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
#' tree1 <- hclu_hierarclust(dissim, 
#'                           n_clust = 10:15, 
#'                           index = "Simpson")
#' tree1
#' 
#' a <- bioregionalization_metrics(tree1, 
#'                                 dissimilarity = dissim, 
#'                                 net = comnet,
#'                                 site_col = "Node1", 
#'                                 species_col = "Node2",
#'                                 eval_metric = c("tot_endemism", 
#'                                                 "avg_endemism",
#'                                                 "pc_distance", 
#'                                                 "anosim"))
#' a
#'
#' @importFrom stats as.dist
#' @importFrom data.table chmatch as.data.table
#' @import data.table 
#'
#'@export
bioregionalization_metrics <- function(bioregionalization, 
                                       dissimilarity = NULL,
                                       dissimilarity_index = NULL, 
                                       net = NULL,
                                       site_col = 1, 
                                       species_col = 2,
                                       eval_metric = "all"){
  
  dissimilarity_based_metrics <- c("pc_distance", 
                                    "anosim")
  compo_based_metrics <- c("avg_endemism", 
                           "tot_endemism")
  
  # Control bioregionalization
  if (inherits(bioregionalization, "bioregion.clusters")) {
    if (inherits(bioregionalization$clusters, "data.frame")) {
      has.clusters <- TRUE # To remove? Does not seem relevant anymore
    } else {
      if (bioregionalization$name == "hclu_hierarclust") {
        stop(paste0("No clusters have been generated for your hierarchical ",
                    "tree, please extract clusters from the tree before using ",
                    "bioregionalization_metrics().\n",
                    "See ?hclu_hierarclust or ?cut_tree"), 
             call. = FALSE)
      } else {
        stop(paste0("bioregionalization does not have the expected type of ",
                    "'clusters' slot"), 
             call. = FALSE)
      }
    }
  } else {
    stop(paste0("This function is designed to work on bioregion.clusters ",
                "objects (outputs from clustering functions)"), 
         call. = FALSE)
    # Add here the possibility to work on data.frame / matrices of clusters
    # directly
  } 
  
  # Control eval_metrics
  controls(args = eval_metric, data = NULL, type = "character_vector")
  if ("all" %in% eval_metric) {
    eval_metric <- c(dissimilarity_based_metrics, 
                     compo_based_metrics)
  }
  if (length(intersect(c(dissimilarity_based_metrics, compo_based_metrics), 
                       eval_metric)) != length(eval_metric)) {
    stop(paste0("One or several metric(s) chosen are not", 
                " available.\n",
                "Please choose from the following:\n",
                "pc_distance, anosim, avg_endemism or tot_endemism."),
         call. = FALSE)
  }

  # Control dissimilarity_index
  if(!is.null(dissimilarity_index)){
    controls(args = dissimilarity_index, type = "character")
  }
  
  # Control dissimilarity
  if (is.null(dissimilarity)) {
    has.dissimilarity <- FALSE
    if(any(eval_metric %in% dissimilarity_based_metrics)){
      warning(paste0("No dissimilarity oject provided, so metrics ",
                     paste(eval_metric[which(eval_metric %in%
                                             dissimilarity_based_metrics)],
                           collapse = ", "),
                     " will not be computed.\n"))
      eval_metric <-
        eval_metric[-which(eval_metric %in% dissimilarity_based_metrics)]
    }
  } else if (inherits(dissimilarity, "bioregion.pairwise")) {
    if (attr(dissimilarity, "type") == "dissimilarity") {
      if(is.null(dissimilarity_index)) {
        if("index" %in% names(bioregionalization$args)) {
          # If an index was already used for bioregionalization, use it here
          dissimilarity_index <- bioregionalization$args$index 
        } else {
          # else choose the 3rd column and warn the user if there are more than
          # 3 columns
          if(ncol(dissimilarity) > 3) {
            warning(paste0("You did not specify the dissimilarity ",
                           "index to use in dissimilarity.",
                           " Defaulting to the third column of",
                           " the dissimilarity object: ",
                           colnames(dissimilarity)[3]))
          }
          dissimilarity_index <- colnames(dissimilarity)[3]
        }
        
      } else if(!(dissimilarity_index %in% colnames(dissimilarity))) {
        stop(paste0("dissimilarity_index does not exist in the dissimilarity ",
                    "object. Did you misspecify the metric name?"),
             call. = FALSE)
      }
      dist_object <- stats::as.dist(
        net_to_mat(dissimilarity[, c(colnames(dissimilarity)[1:2],
                                     dissimilarity_index)],
                   weight = TRUE, squared = TRUE, symmetrical = TRUE))
      has.dissimilarity <- TRUE
    } else{
      stop(paste0("dissimilarity must be an object containing dissimilarity ",
                  "indices from dissimilarity() or ",
                  "similarity_to_dissimilarity(), or an object of class dist."),
           call. = FALSE)
    }
  } else if(!any(inherits(dissimilarity, "bioregion.pairwise"),
                 inherits(dissimilarity, "dist"))){
    #if(is.numeric(dissimilarity_index)){
    #  dissimilarity_index <- names(dissimilarity)[dissimilarity_index]
    #}
    #dist_object <- stats::as.dist(
    #  net_to_mat(dissimilarity[, c(colnames(dissimilarity)[1:2],
    #                               dissimilarity_index)],
    #             weight = TRUE, squared = TRUE, symmetrical = TRUE))
    #has.dissimilarity <- TRUE
    #
    #if(!(dissimilarity_index %in% colnames(dissimilarity))){
    #  stop(paste0("dissimilarity is not a bioregion.pairwise object, ",
    #              "a dissimilarity matrix (class dist) or a data.frame with ",
    #              "at least 3 columns (site1, site2, and your dissimilarity ",
    #              "index)"),
    #       call. = FALSE)
    #}
    
    stop(paste0("dissimilarity must be a bioregion.pairwise object or ",
                "a dissimilarity matrix (class dist)."),
        call. = FALSE)
  }
  
  if(has.dissimilarity){
    if(attr(dist_object, "Size") != bioregionalization$inputs$nb_sites){
      stop(paste0("bioregionalization and dissimilarity have different ",
                  "number of sites."),
           call. = FALSE)
    }
  }
  
  if (is.null(net)) {
    has.contin <- FALSE
    if(any(eval_metric %in% compo_based_metrics)){
      warning(paste0("No site-species network provided, so metrics ",
                     paste(eval_metric[which(eval_metric %in%
                                               compo_based_metrics)],
                           collapse = ", "), 
                     " will not be computed\n"))
      eval_metric <- eval_metric[-which(eval_metric %in% compo_based_metrics)]
    }
  } else {
    if(site_col == species_col){
      stop("site_col and species_col should not be the same.", call. = FALSE)
    }
    controls(args = site_col, data = net, type = "input_net_bip_col")
    controls(args = species_col, data = net, type = "input_net_bip_col")
    if(has.clusters){
      if(any(!(bioregionalization$clusters$ID %in% c(net[, site_col],
                                                 net[, species_col])))) {
        stop(paste0("Some elements of the cluster table (column ID) cannot be ",
                    "found in net."),
             call. = FALSE)
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
    stop(paste0("No evaluation metric can be computed because of missing ",
                "arguments. Check arguments dissimilarity and sp_site_table"), 
                call. = FALSE)
  }
  
  nb_sites <- bioregionalization$inputs$nb_sites
  
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
    nrow = ncol(bioregionalization$clusters) - 1,
    ncol = 2 + length(eval_metric),
    dimnames = list(colnames(
      bioregionalization$clusters)[2:ncol(bioregionalization$clusters)],
      c("K", "n_clusters", eval_metric))))
  
  evaluation_df$K <-
    colnames(bioregionalization$clusters)[2:ncol(bioregionalization$clusters)]
  
  evaluation_df$n_clusters <- apply(
    bioregionalization$clusters[, 2:(ncol(bioregionalization$clusters)), drop = FALSE],
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
  #   bioregionalization$clusters[match(net_long$ID,
  #                                 bioregionalization$clusters$ID), -1])
  
  if(has.dissimilarity & any(c("pc_distance", "anosim") %in% eval_metric)){
    message("Computing similarity-based metrics...")
    # The next line will create, for each element of the dissimilarity matrix,
    # a vector indicating whether if each dissimilarity is within or between
    # clusters
    dissimilarity <- data.frame(
      dissimilarity,
      bioregionalization$clusters[data.table::chmatch(dissimilarity$Site1,
                                                  bioregionalization$clusters$ID), 
                              bioregionalization$cluster_info$partition_name,
                              drop = FALSE],
      bioregionalization$clusters[data.table::chmatch(dissimilarity$Site2,
                                                  bioregionalization$clusters$ID), 
                              bioregionalization$cluster_info$partition_name,
                              drop = FALSE])
    
    dissimilarity[, bioregionalization$cluster_info$partition_name] <- 
      dissimilarity[, bioregionalization$cluster_info$partition_name] ==
      dissimilarity[, paste0(bioregionalization$cluster_info$partition_name, ".1")]
    
    dissimilarity <-
      dissimilarity[,
                    -which(colnames(dissimilarity) %in%
                             paste0(bioregionalization$cluster_info$partition_name,
                                    ".1"))]
    
    if("pc_distance" %in% eval_metric) {
      evaluation_df$pc_distance <-
        vapply(bioregionalization$cluster_info$partition_name,
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
        bioregionalization$cluster_info$partition_name,
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
      bioregionalization$clusters[data.table::chmatch(net[, site_col],
                                                  bioregionalization$clusters$ID),
                              -1])
    
    # Correcting column names when there is only one clustering
    if("bioregionalization.clusters.data.table..chmatch.net...site_col..." %in%
       colnames(net)){
      colnames(net)[colnames(net) ==
                      "bioregionalization.clusters.data.table..chmatch.net...site_col..."] <-
        colnames(bioregionalization$clusters)[2]
    }

    # Visible binding for global variable
    N <- endemism <- end_richness <- pc_endemism <- NULL 
    
    # Fast calculation of endemism per cluster
    endemism_results <- lapply(
      bioregionalization$cluster_info$partition_name,
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
    names(endemism_results) <- bioregionalization$cluster_info$partition_name
    
    # Average endemism per cluster
    if("avg_endemism" %in% eval_metric){
      evaluation_df$avg_endemism <- vapply(
        bioregionalization$cluster_info$partition_name,
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
        bioregionalization$cluster_info$partition_name,
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
  
  class(outputs) <- append("bioregion.bioregionalization.metrics", 
                           class(outputs))
  return(outputs)
}
