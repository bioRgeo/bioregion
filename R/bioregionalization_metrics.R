#' Calculate metrics for a bioregionalization
#'
#' This function calculates metrics at the bioregionalization level. These 
#' evaluation metrics can be used to assess and select a partition based on the 
#' optimal number of clusters.
#'
#' @param bioregionalization A `bioregion.clusters` object.
#' 
#' @param eval_metrics A `character` vector or a single `character` string 
#' indicating the metric(s) to be calculated. Available options are 
#' `"prop_between_dissim"`, `"anosim"`, `"mean_endemics"` or `"tot_endemics"`. Use `"all"` 
#' to compute all available metrics. See Details for metric descriptions.
#' 
#' @param dissimilarity A site-by-site dissimilarity object from [dissimilarity()] 
#' or [dissimilarity_to_similarity()]. Required only for `"prop_between_dissim"` 
#' and `"anosim"`.
#' 
#' @param dissimilarity_index The name or number of the column to use as 
#' dissimilarity. By default, the third column name of `dissimilarity` is used.
#' 
#' @param comat A site-species `matrix` with sites as rows and species as
#' columns. Should be provided if `eval_metrics` includes `"avg_endemism"` or 
#' `"tot_endemism"`.
#'
#' @param anosim_permutations The number of permutations used to compute the
#' p-value associated with the ANOSIM statistic. Defaults to 1, in which
#' case no p-value is returned.
#' 
#' @param eval_metric Deprecated.
#' 
#' @param net Deprecated.
#' 
#' @param site_col Deprecated.
#' 
#' @param species_col Deprecated.
#' 
#' @return A \code{data.frame} containing the \code{eval_metrics} values for a 
#' bioregionalization and its partition(s).
#'
#' @details
#' **Evaluation metrics:**
#' 
#' \itemize{
#' 
#' \item{[**prop_between_dissim**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#prop_between_dissim): 
#' The proportion of total dissimilarity occurring between
#' bioregions, following Holt et al. (2013), calculated as the sum of
#' between-bioregion dissimilarities divided by the total sum of 
#' dissimilarities.}
#' 
#' \item{[**anosim**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#anosim): 
#' The Analysis of Similarities (ANOSIM) statistic, based on
#' [vegan::anosim()][vegan::anosim]. It measures the separation between
#' within- and between-bioregion dissimilarities. When `permutation > 1`, 
#' a p-value is calculated using permutation testing.}
#' 
#' \item{[**mean_endemics**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#mean_endemics): 
#' The mean proportion of endemic species across bioregions,
#' following Kreft & Jetz (2010). For each bioregion, the proportion of endemic
#' species is calculated and then averaged across bioregions.}
#'
#' \item{[**tot_endemics**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#tot_endemics): 
#' The proportion of endemic species across all bioregion. 
#' It is calculated as the total number of endemic species divided by the total 
#' number of species. Endemic species are those occurring in only one bioregion.}
#' 
#' }
#' 
#' @references
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
#' \url{https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#bioregionalization}.
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
#' dissim <- dissimilarity(comat, metric = "all")
#' 
#' # User-defined number of clusters
#' bioreg <- hclu_hierarclust(dissim, 
#'                            n_clust = 10:15, 
#'                            index = "Simpson",
#'                            verbose = FALSE)
#' 
#' met <- bioregionalization_metrics(bioreg,
#'                                   eval_metrics = "all",
#'                                   dissimilarity = dissim,
#'                                   comat = comat)
#' met
#'
#'@export
bioregionalization_metrics <- function(bioregionalization,
                                       eval_metrics = c("prop_between_dissim", "anosim"), 
                                       dissimilarity,
                                       dissimilarity_index = names(dissimilarity)[3], 
                                       comat = NULL,
                                       anosim_permutations = 1,
                                       eval_metric = NULL,
                                       net = NULL,
                                       site_col = NULL, 
                                       species_col = NULL){
  
  # Control deprecated
  if (!is.null(eval_metric)) {
    warning("eval_metric is deprecated.", 
            call. = FALSE)
  }
  if (!is.null(net)) {
    warning("net is deprecated.", 
            call. = FALSE)
  }
  if (!is.null(site_col)) {
    warning("site_col is deprecated.", 
            call. = FALSE)
  }
  if (!is.null(species_col)) {
    warning("species_col is deprecated.", 
            call. = FALSE)
  }
  
  # Convert metrics names
  eval_metrics <- convert_metric_names(eval_metrics)
  
  # Control bioregionalization
  controls(args = NULL, 
           data = bioregionalization, 
           type ="input_bioregionalization")
  
  # Extract node_type
  b_node_type <- bioregionalization$inputs$node_type 
  if(b_node_type == "species"){
    stop(paste0("No bioregion are assigned to the site in bioregionalization."), 
         call. = FALSE)
  }
  b_site <- bioregionalization$clusters[attr(bioregionalization$clusters, 
                                             "node_type") == "site", 1]
  
  # Controls metrics
  dissim_metrics <- c("prop_between_dissim", 
                      "anosim")
  comat_metrics <- c("mean_endemics", 
                     "tot_endemics")
  controls(args = eval_metrics, data = NULL, type = "character_vector")
  metrics <- eval_metrics
  if ("all" %in% eval_metrics) {
    metrics <- c(dissim_metrics, comat_metrics)
  }
  if (length(intersect(c(dissim_metrics, comat_metrics), metrics)) !=
      length(metrics)) {
    stop(paste0("One or several evaluation metrics chosen are not", 
                " available.\n",
                "Please choose from the following:\n",
                "prop_between_dissim, anosim, mean_endemics and tot_endemics"),
         call. = FALSE)
  }
  
  # Check if comat and/or similarity are needed
  dissim_needed <- FALSE
  if(length(intersect(dissim_metrics, metrics))>0){
    dissim_needed <- TRUE
  }
  comat_needed <- FALSE
  if(length(intersect(comat_metrics, metrics))>0){
    comat_needed <- TRUE
  }
  
  if(is.null(dissimilarity) & is.null(comat)){
    stop(paste0("At least dissimilarity or comat should be provided."),
         call. = FALSE)
  }
  
  if(is.null(dissimilarity) & dissim_needed){
    warning(paste0("Some metrics (", 
                   paste(intersect(metrics, dissim_metrics), collapse = ", "),
                   ") will be skipped because no dissimilarity is provided."),
            call. = FALSE)
    metrics <- setdiff(metrics, dissim_metrics)
    dissim_needed <- FALSE
  }
  
  if(is.null(comat) & comat_needed){
    warning(paste0("Some metrics (", 
                   paste(intersect(metrics, comat_metrics), collapse = ", "),
                   ") will be skipped because no co-occurrence matrix is provided."),
            call. = FALSE)
    metrics <- setdiff(metrics, comat_metrics)
    comat_needed <- FALSE
  }
  
  # Stop if no metrics
  if(length(metrics) == 0){
    stop(paste0("At least one metric with the appropriate input ", 
                "should be specified."),
         call. = FALSE)
  }
  
  # Control permutations and update metrics
  if("anosim" %in% metrics){
    controls(args = anosim_permutations, data = NULL, type = "positive_integer")
    if(anosim_permutations > 1){
      metrics <- append(metrics, "anosim_pval", 
                        after = which(metrics == "anosim"))
    }
  }
  
  # Control dissimilarity if needed
  if(dissim_needed){
    
    controls(args = NULL, 
             data = dissimilarity, 
             type = "input_conversion_dissimilarity")
    controls(args = dissimilarity_index, 
             data = dissimilarity, 
             type = "input_net_index")
    
    dissimilarity <- dissimilarity
    dissimilarity[,3] <- dissimilarity[, dissimilarity_index]
    dissimilarity <- dissimilarity[,1:3]
    dissimilarity <- net_to_mat(dissimilarity, 
                                weight = TRUE, 
                                squared = TRUE,
                                symmetrical = TRUE)
    
    dissim_site <- rownames(dissimilarity)
    
    # Check that dissim_site are in bioregionalization 
    missing_sites <- setdiff(b_site, dissim_site)
    if(length(missing_sites) > 0){
      stop(paste0("Some sites are not found in dissimilarity:\n",
                  "  Missing sites: ", paste(utils::head(missing_sites, 10), collapse = ", "),
                  if(length(missing_sites) > 10) paste0(" ... (", length(missing_sites) - 10, " more)") else "",
                  "\n  Please ensure that all sites in 'bioregionalization' have corresponding entries in 'dissimilarity'."),
           call. = FALSE)
    }
    dissimilarity <- dissimilarity[match(b_site, dissim_site), 
                                   match(b_site, dissim_site)]
    dissim_site <- rownames(dissimilarity)
    
  }  
  
  # Control comat if needed
  if(comat_needed){
    
    controls(args = NULL, data = comat, type = "input_matrix")
    minco <- min(comat)
    if (minco < 0) {
      stop("Negative value(s) detected in comat!", 
           call. = FALSE)
    }
    comat_site <- rownames(comat)
    
    # Check that comat_site are in bioregionalization 
    missing_sites <- setdiff(b_site, comat_site)
    if(length(missing_sites) > 0){
      stop(paste0("Some sites are not found in comat:\n",
                  "  Missing sites: ", paste(utils::head(missing_sites, 10), collapse = ", "),
                  if(length(missing_sites) > 10) paste0(" ... (", length(missing_sites) - 10, " more)") else "",
                  "\n  Please ensure that all sites in 'bioregionalization' have corresponding entries in 'comat'."),
           call. = FALSE)
    }
    comat <- comat[match(b_site, comat_site), ]
    comat_site <- rownames(comat)
    
  } 
  
  # Prepare outputs
  output <- bioregionalization$cluster_info[,1:2]
  matmet <- matrix(0, dim(output)[1], length(metrics))
  
  # Use bioregion_metrics to compute Endemics if necessary
  if("mean_endemics" %in% metrics |
     "tot_endemics" %in% metrics){
    
    endemics <- bioregion_metrics(bioregionalization, comat)
    if(is.data.frame(endemics)){ # Use a list if one bioregionalization
      endemics <- list(endemics)
    }
  
  }
  
  # Loop over partitions
  bioregionalization <- bioregionalization$clusters
  nb_partitions <- dim(bioregionalization)[2] - 1
  
  for(k in 1:nb_partitions){
    
    # Partition k
    partk <- bioregionalization[attr(bioregionalization, "node_type") == "site", 
                                (k+1)]
    
    # Initialize counter metrics
    l <- 0
    
    # prop_between_dissim
    if("prop_between_dissim" %in% metrics){
      
      dem <- sum(dissimilarity)
      num <- dissimilarity
      num[outer(partk, partk, "==")] <- 0
      num <- sum(num)
      
      l <- l+1
      matmet[k,l] <- num / dem
      
    }
    
    # anosim
    if("anosim" %in% metrics){
      
      if(output$n_clust[k] > 1){
        anosim <- vegan::anosim(dissimilarity, 
                                partk, 
                                permutations = anosim_permutations)
        stat <- anosim$statistic
        pval <- anosim$signif
      }else{
        stat <- NA
        pval <- NA
      }
      
      l <- l+1
      matmet[k,l] <- stat
      
      if(anosim_permutations > 1){
        l <- l+1
        matmet[k,l] <- pval
      }
      
    }
    
    # mean_endemics
    if("mean_endemics" %in% metrics){
      l <- l+1
      matmet[k,l] <- mean(endemics[[k]]$rich_endemics / endemics[[k]]$richness)
    }
    
    # tot_endemics
    if("tot_endemics" %in% metrics){
      l <- l+1
      matmet[k,l] <- sum(endemics[[k]]$rich_endemics) / 
                     sum(apply(comat, 2, sum) > 0)
    }
      
    
  }  
  
  # Return output
  output <- cbind(output, matmet)
  colnames(output) <- c("partition", "n_bioregions", metrics)
  class(output) <- append("bioregion.bioregionalization.metrics", 
                            class(output))
  return(output)
  
}