#' Compare cluster memberships among multiple bioregionalizations
#' 
#' This function computes pairwise comparisons for several
#' bioregionalizations, usually outputs from `netclu_`, `hclu_`, or `nhclu_`
#' functions. It also provides the confusion matrix from pairwise comparisons, 
#' enabling the user to compute additional comparison metrics.
#' 
#' @param cluster_object A `bioregion.clusters` object, a `data.frame`, or a 
#' list of `data.frame` objects containing multiple bioregionalizations. At 
#' least two bioregionalizations are required. If a list of `data.frame` is 
#' provided, all `data.frame` objects must have the same number of rows (i.e., 
#' the same items in the clustering for all bioregionalizations).
#' 
#' @param indices `NULL` or `character`. Indices to compute for the pairwise
#' comparison of bioregionalizations. Currently available metrics are `"rand"`
#' and `"jaccard"`.
#' 
#' @param cor_frequency A `boolean`. If `TRUE`, computes the correlation 
#' between each bioregionalization and the total frequency of co-membership of
#' items across all bioregionalizations. This is useful for identifying which
#' bioregionalization(s) is(are) most representative of all computed
#' bioregionalizations.
#' 
#' @param store_pairwise_membership A `boolean`. If `TRUE`, stores the pairwise 
#' membership of items in the output object.
#' 
#' @param store_confusion_matrix A `boolean`. If `TRUE`, stores the confusion 
#' matrices of pairwise bioregionalization comparisons in the output object.
#' 
#' @details 
#' This function operates in two main steps:
#' 
#' 1. Within each bioregionalization, the function compares all pairs of items 
#' and documents whether they are clustered together (`TRUE`) or separately 
#' (`FALSE`). For example, if site 1 and site 2 are clustered in the same 
#' cluster in bioregionalization 1, their pairwise membership `site1_site2` 
#' will be `TRUE`. This output is stored in the `pairwise_membership` slot if 
#' `store_pairwise_membership = TRUE`.
#' 
#' 2. Across all bioregionalizations, the function compares their pairwise 
#' memberships to determine similarity. For each pair of bioregionalizations, 
#' it computes a confusion matrix with the following elements:
#'  * `a`: Number of item pairs grouped in both bioregionalizations.
#'  * `b`: Number of item pairs grouped in the first but not in the second 
#' bioregionalization.
#'  * `c`: Number of item pairs grouped in the second but not in the first 
#' bioregionalization.
#'  * `d`: Number of item pairs not grouped in either bioregionalization.
#'  
#' The confusion matrix is stored in `confusion_matrix` if 
#' `store_confusion_matrix = TRUE`.
#' 
#' Based on these confusion matrices, various indices can be computed to 
#' measure agreement among bioregionalizations. The currently implemented 
#' indices are:
#' 
#' * **Rand index**: `(a + d) / (a + b + c + d)`
#'   Measures agreement by considering both grouped and ungrouped item pairs.
#' 
#' * **Jaccard index**: `a / (a + b + c)`
#'   Measures agreement based only on grouped item pairs.
#' 
#' These indices are complementary: the Jaccard index evaluates clustering 
#' similarity, while the Rand index considers both clustering and separation. 
#' For example, if two bioregionalizations never group the same pairs, their 
#' Jaccard index will be 0, but their Rand index may be > 0 due to ungrouped 
#' pairs.
#' 
#' Users can compute additional indices manually using the list of confusion 
#' matrices.
#' 
#' To identify which bioregionalization is most representative of the others, 
#' the function can compute the correlation between the pairwise membership of 
#' each bioregionalization and the total frequency of pairwise membership across 
#' all bioregionalizations. This is enabled by setting `cor_frequency = TRUE`.
#' 
#' @return 
#' A `list` containing 4 to 7 elements:
#' 
#' \enumerate{
#' \item{**args**: A `list` of user-provided arguments.}
#' \item{**inputs**: A `list` containing information on the input 
#' bioregionalizations, such as the number of items clustered.}
#' \item{**pairwise_membership** (optional): If `store_pairwise_membership = TRUE`, 
#' a `boolean matrix` where `TRUE` indicates two items are in the same cluster, 
#' and `FALSE` indicates they are not.}
#' \item{**freq_item_pw_membership**: A `numeric vector` containing the number of 
#' times each item pair is clustered together, corresponding to the sum of rows 
#' in `pairwise_membership`.}
#' \item{**bioregionalization_freq_cor** (optional): If `cor_frequency = TRUE`, 
#' a `numeric vector` of correlations between individual bioregionalizations 
#' and the total frequency of pairwise membership.}
#' \item{**confusion_matrix** (optional): If `store_confusion_matrix = TRUE`, 
#' a `list` of confusion matrices for each pair of bioregionalizations.}
#' \item{**bioregionalization_comparison**: A `data.frame` containing comparison 
#' results, where the first column indicates the bioregionalizations compared, 
#' and the remaining columns contain the requested `indices`.}}
#'  
#' @seealso 
#' For more details illustrated with a practical example, 
#' see the vignette: 
#' \url{https://biorgeo.github.io/bioregion/articles/a5_2_compare_bioregionalizations.html}.
#' 
#' Associated functions: 
#' [bioregionalization_metrics]
#'  
#' @author
#' Boris Leroy (\email{leroy.boris@gmail.com}) \cr
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) \cr
#' Pierre Denelle (\email{pierre.denelle@gmail.com})
#' 
#' @examples 
#' # A simple case with four bioregionalizations of four items
#' bioregionalizations <- data.frame(matrix(nr = 4, 
#'                                          nc = 4, 
#'                                          c(1,2,1,1,1,2,2,1,2,1,3,1,2,1,4,2),
#'                                          byrow = TRUE))
#' bioregionalizations
#' 
#' compare_bioregionalizations(bioregionalizations)
#' 
#' # Find out which bioregionalizations are most representative
#' compare_bioregionalizations(bioregionalizations,
#'                             cor_frequency = TRUE)
#'                                 
#' @export
compare_bioregionalizations <- function(cluster_object,
                                        indices = c("rand", "jaccard"),
                                        cor_frequency = FALSE,
                                        store_pairwise_membership = TRUE,
                                        store_confusion_matrix = TRUE){
  
  # input can be of format bioregion.clusters
  if (inherits(cluster_object, "bioregion.clusters")) {
    if (inherits(cluster_object$clusters, "data.frame")) {
      has.clusters <- TRUE
      clusters <- cluster_object$clusters
    } else {
      if (cluster_object$name == "hclu_hierarclust") {
        stop(paste0("No clusters have been generated for your hierarchical ",
                    "tree, please extract clusters from the tree before using ",
                    "bioregionalization_metrics().\n",
                    "See ?hclu_hierarclust or ?cut_tree"), 
             call. = FALSE)
      } else {
        stop(paste0("cluster_object does not have the expected type of ",
                    "'clusters' slot"), 
             call. = FALSE)
      }
    }
    # or data.frame (next row)
  } else if (inherits(cluster_object, "data.frame")) {
    controls(data = cluster_object, type = "input_data_frame")
    clusters <- cluster_object
    # or list of data.frames
  } else if (inherits(cluster_object, "list")) {
    # test that all elements of the list are data.frames
    test <- lapply(cluster_object, FUN = function(x) {
      try(controls(data = x, type = "input_data_frame"), silent = TRUE)})
    if (any(vapply(test, function(x) inherits(x, "try-error"),
                   FUN.VALUE = logical(1)))) {
      stop(paste0("All elements in cluster_object should be ",
                  "of data.frame format"), 
                  call. = FALSE)
      # ensure that the number of rows are identical among df of the list
    } else if (length(unique(vapply(cluster_object,
                                    nrow,
                                    FUN.VALUE = numeric(1)))) == 1){
      clusters <- data.frame(cluster_object)
    }
    # if none of the above then stop
  }  else {
    stop(paste0("This function is designed to work either on ",
                "bioregion.clusters objects, on a data.frame or on a list of ",
                "data.frames"), 
         call. = FALSE)
  }
  
  if(ncol(clusters) == 1) {
    stop(paste0("This function is designed to be applied on multiple ",
                "bioregionalizations. Your cluster_object only has a single ",
                "bioregionalization (one column)."), 
         call. = FALSE)
  }
  
  if(!is.null(indices)) {
    controls(args = indices, type = "character_vector")
    if(!isTRUE(unique(indices %in% c("rand", "jaccard")))) {
      stop(paste0("Please choose indices from the following:\n",
                  "rand or jaccard."), 
           call. = FALSE)
    }
  }
  
  controls(args = cor_frequency, type = "boolean")
  controls(args = store_pairwise_membership, type = "boolean")
  controls(args = store_confusion_matrix, type = "boolean")

  message(Sys.time(), 
          " - Computing pairwise membership comparisons for each ",
          "bioregionalization...\n")
  # if(!is.null(sample_items)) {
  #   if(ncol(clusters) * (nrow(clusters) * (nrow(clusters) - 1)) / 2 > 10e6) {
  #     message("       /!\\\ NOTE: Very high number of comparisons ",
  #             "(", nrow(clusters), " rows in clusters --> ",
  #             (nrow(clusters) * (nrow(clusters) - 1)) / 2, 
  #             " pairwise comparisons per bioregionalization --> ",
  #             ncol(clusters) * (nrow(clusters) * (nrow(clusters) - 1)) / 2, 
  #             " comparisons in total across all bioregionalizations).\n      ",
  #             "If the computation ",
  #             "time is too long, consider using argument sample_items and",
  #             " set it to a reasonable number (e.g. 5000).")
  #   }
  # }
  # 
  
  # WITHIN bioregionalization - Pairwise membership comparison ----------------
  
  # For each bioregionalization, compare all pairs of items to detect if they
  # are in the same cluster or not
  item_pw_mb <- get_pairwise_membership(
    clusters)
  
  # BETWEEN bioregionalizationS - Pairwise bioregionalization comparison ------
  
  message(Sys.time(), 
          " - Comparing memberships among bioregionalizations...\n")
  
  # Prepare the pairwise bioregionalization comparisons
  partnames <- pw_cluster_comps <- t(utils::combn(seq_len(ncol(clusters)), 2))
  
  
  # Compute confusion matrices for all pairs of bioregionalizations
  all_conf_matrices <- lapply(seq_len(nrow(pw_cluster_comps)), function(x) {
    get_confusion_matrix(item_pw_mb[, pw_cluster_comps[x, 1]],
                         item_pw_mb[, pw_cluster_comps[x, 2]])
  })
  
  # partnames gets the names of each bioregionalization 
  partnames[] <- colnames(clusters)[pw_cluster_comps]
  names(all_conf_matrices) <- apply(partnames, 1, paste, collapse = "%")
  
  # Item membership frequency -------------------------------------------------
  item_pw_mb_freq <- rowSums(item_pw_mb)
  
  # Prepare data.frame comparing bioregionalizations --------------------------
  
  partcomp_indices <- data.frame(
    bioregionalization_comparison = names(all_conf_matrices))
  
  # Rand index ----------------------------------------------------------------
  if("rand" %in% indices) {
    message(Sys.time(), 
            " - Computing Rand index...\n")
    
    partcomp_indices$rand <- vapply(all_conf_matrices,
                                    rand_index,
                                    numeric(1))
  }
  
  # Jaccard index -------------------------------------------------------------
  if("jaccard" %in% indices) {
    message(Sys.time(), 
            " - Computing Jaccard index...\n")
    partcomp_indices$jaccard <- vapply(all_conf_matrices,
                                       jaccard_index,
                                       numeric(1))
  }
  
  # Point biserial correlation  -----------------------------------------------
  if(cor_frequency) {
    message(
      Sys.time(), 
      " - Computing the correlation between each bioregionalization and the",
      " vector of frequency of pairwise membership...\n")
    
    bioregionalization_freq_cor <- as.vector(suppressWarnings(
      cor(item_pw_mb, item_pw_mb_freq)))
    bioregionalization_freq_cor[is.na(bioregionalization_freq_cor)] <- 0
    names(bioregionalization_freq_cor) <- colnames(clusters)
  }
  
  # Store outputs -----------------------------------------------------------
  outputs <- list(
    args = list(indices = indices,
                cor_frequency = cor_frequency,
                store_pairwise_membership = store_pairwise_membership,
                store_confusion_matrix = store_confusion_matrix))
  
  outputs$inputs <- c(number_items = nrow(clusters),
                      number_bioregionalizations = ncol(clusters))
  
  if(store_pairwise_membership) {
    outputs$pairwise_membership <- item_pw_mb
  }
  
  outputs$freq_item_pw_membership <- item_pw_mb_freq
  
  if(cor_frequency) {
    outputs$bioregionalization_freq_cor <-   bioregionalization_freq_cor
  }
  
  if(store_confusion_matrix) {
    outputs$confusion_matrix <- all_conf_matrices
  }
  
  outputs$bioregionalization_comparison <- partcomp_indices
  
  class(outputs) <- append("bioregion.bioregionalization.comparison", class(outputs))
  return(outputs)
}

# Get pairwise membership -----------------------------------------------------
get_pairwise_membership <- function(input_clusters) {
  
  # Number of rows (items) and columns (iterations or clustering solutions)
  n_items <- nrow(input_clusters)
  n_bioregionalizations <- ncol(input_clusters)
  
  # Initialize matrix to store pairwise membership results
  pw_membership <- matrix(FALSE, nrow = choose(n_items, 2),
                          ncol = n_bioregionalizations)
  
  # Get all unique pairwise comparisons indices
  pairwise_comps <- t(utils::combn(seq_len(n_items), 2))
  
  # Compute pairwise membership
  for (col in seq_len(n_bioregionalizations)) {
    # Extract cluster memberships for this bioregionalization
    cluster_col <- input_clusters[, col]
    
    # Compare memberships directly for each pair of items
    pw_membership[, col] <- cluster_col[pairwise_comps[, 1]] == 
      cluster_col[pairwise_comps[, 2]]
  }
  
  # Set row names based on the pairwise comparisons (optional)
  rownames(pw_membership) <- apply(pairwise_comps, 1, paste, collapse = "_")
  colnames(pw_membership) <- colnames(input_clusters)
  
  return(pw_membership)
}

# Confusion matrix ------------------------------------------------------------
get_confusion_matrix <- function(clust1, clust2) {
  return(c(a = length(which(clust1 & clust2)), # TP
           b = length(which(clust1 & !clust2)),
           c = length(which(!clust1 & clust2)),
           d = length(which(!clust1 & !clust2)))) # TN
}

# Rand index ------------------------------------------------------------------
rand_index <- function(confusion_matrix) {
  RI <- (confusion_matrix["a"] +
           confusion_matrix["d"]) /
    (confusion_matrix["a"] + confusion_matrix["b"] + confusion_matrix["c"] +
       confusion_matrix["d"])
  names(RI) <- "rand"
  return(RI)
}

# Jaccard index ---------------------------------------------------------------
jaccard_index <- function(confusion_matrix) {
  JI <- confusion_matrix["a"] /
    (confusion_matrix["a"] + confusion_matrix["b"] + confusion_matrix["c"])
  names(JI) <- "jaccard"
  return(JI)
}
