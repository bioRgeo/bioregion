#' Compare cluster memberships among multiple bioregionalizations
#' 
#' This function aims at computing pairwise comparisons for several
#' bioregionalizations, usually an output from `netclu_`, `hclu_` or `nhclu_`
#' functions. It also provides the confusion matrix from pairwise comparisons, 
#' so that the user can compute additional comparison metrics.
#' 
#' @param cluster_object a `bioregion.clusters` object or a `data.frame` or a 
#' list of `data.frame` containing multiple bioregionalizations. At least two
#' bioregionalizations are required. If a list of `data.frame` is provided,
#' they should all have the same number of rows (i.e., same items in the
#' clustering for all bioregionalizations). 
#' 
#' @param indices `NULL` or `character`. Indices to compute for the pairwise
#' comparison of bioregionalizations. Current available metrics are `"rand"`
#' and `"jaccard"`.
#' 
#' @param cor_frequency a `boolean`. If `TRUE`, then computes the correlation 
#' between each bioregionalization and the total frequency of co-membership of
#' items across all bioregionalizations. Useful to identify which
#' bioregionalization(s) is(are) most representative of all the computed
#' bioregionalizations.
#' 
#' @param store_pairwise_membership a `boolean`. If `TRUE`, the pairwise 
#' membership of items is stored in the output object.
#' 
#' @param store_confusion_matrix a `boolean`. If `TRUE`, the confusion matrices 
#' of pairwise bioregionalization comparisons are stored in the output object.
#' 
#' @details 
#' This function proceeds in two main steps:
#' 
#' 1. The first step is done within each bioregionalization. It will compare
#' all pairs of items and document if they are clustered together (`TRUE`) or
#' separately (`FALSE`) in each bioregionalization. For example, if site 1 and
#' site 2 are clustered in the same cluster in bioregionalization 1, then the
#' pairwise membership site1_site2 will be `TRUE`. The output of this first
#' step is stored in the slot `pairwise_membership` if
#' `store_pairwise_membership = TRUE`.
#' 
#' 2. The second step compares all pairs of bioregionalizations by analysing if
#' their pairwise memberships are similar or not. To do so, for each pair of 
#' bioregionalizations, the function computes a confusion matrix with four
#' elements:
#'  * `a` number of pairs of items grouped in bioregionalization 1 and in
#'  bioregionalization 2
#'  * `b` number of pairs of items grouped in bioregionalization 1 but not in
#'  bioregionalization 2
#'  * `c` number of pairs of items not grouped in bioregionalization 1 but
#'  grouped in bioregionalization 2
#'  * `d` number of pairs of items not grouped in both bioregionalization 1
#'  & 2
#'  
#' The confusion matrix is stored in `confusion_matrix` if 
#' `store_confusion_matrix = TRUE`.
#' 
#' 
#' Based on the confusion matrices, we can compute a range of indices to 
#' indicate the agreement among bioregionalizations. As of now, we have
#' implemented:
#' * `Rand index` (a + d)/(a + b + c + d)
#' The Rand index measures agreement among bioregionalizations by accounting
#' for both the pairs of sites that are grouped, but also the pairs of sites
#' that are not grouped.
#' * `Jaccard index` a/(a + b + c)
#' The Jaccard index measures agreement among bioregionalizations by only
#' accounting for pairs of sites that are grouped.
#' 
#' These two metrics are complementary, because the Jaccard index will tell
#' if bioregionalizations are similar in their clustering structure, whereas
#' the Rand index will tell if bioregionalizations are similar not only in the
#' pairs of items clustered together, but also in terms of the pairs of sites
#' that are not clustered together. For example, take two bioregionalizations
#' which never group together the same pairs of sites. Their Jaccard index will
#' be 0, whereas the Rand index can be > 0 due to the sites that are not
#' grouped together. 
#' 
#' Additional indices can be manually computed by the users on the basis of the
#' list of confusion matrices.
#' 
#' In some cases, users may be interested in finding which of the
#' bioregionalizations is most representative of all bioregionalizations. To
#' find it out, we can compare the pairwise membership of each
#' bioregionalization with the total frequency of pairwise membership across
#' all bioregionalizations. This correlation can be requested with
#' `cor_frequency = TRUE`.
#' 
#' @seealso [bioregionalization_metrics]
#' 
#' @return 
#' A `list` with 4 to 7 elements:
#'  * `args`: arguments provided by the user
#'  * `inputs`: information on the input bioregionalizations, such as the
#'  number of items being clustered
#'  * (facultative) `pairwise_membership`: only if 
#'  `store_pairwise_membership = TRUE`. This 
#'  element contains the pairwise memberships of all items for each 
#'  bioregionalization, in the form of a `boolean matrix` where `TRUE` means
#'  that two items are in the same cluster, and `FALSE` means that two items 
#'  are not in the same cluster
#'  * `freq_item_pw_membership`: A `numeric vector` 
#'  containing the number of times each pair of items are clustered together.
#'  It corresponds to the sum of rows of the table in `pairwise_membership`
#'  * (facultative) `bioregionalization_freq_cor`:  only if
#'  `cor_frequency = TRUE`.
#'  A `numeric vector` indicating the correlation between individual
#'  bioregionalizations and the total frequency of pairwise membership across
#'  all bioregionalizations. It corresponds to the correlation between
#'  individual columns in `pairwise_membership` and `freq_item_pw_membership`
#'  * (facultative) `confusion_matrix`: only if
#'  `store_confusion_matrix = TRUE`.
#'  A `list` containing all confusion matrices between each pair of
#'  bioregionalizations.
#'  * `bioregionalization_comparison`: a `data.frame` containing the results
#'  of the comparison of bioregionalizations, where the first column indicates
#'  which bioregionalizations are compared, and the next columns correspond to
#'  the requested `indices`.
#'  
#' @author
#' Boris Leroy (\email{leroy.boris@gmail.com}) \cr
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) \cr
#' Pierre Denelle (\email{pierre.denelle@gmail.com})
#' 
#' @examples 
#' # A simple case with four bioregionalizations of four items
#' bioregionalizations <- data.frame(matrix(nr = 4, nc = 4, 
#'                                 c(1,2,1,1,1,2,2,1,2,1,3,1,2,1,4,2),
#'                                 byrow = TRUE))
#' bioregionalizations
#' compare_bioregionalizations(bioregionalizations)
#' 
#' # Find out which bioregionalizations are most representative
#' compare_bioregionalizations(bioregionalizations,
#'                    cor_frequency = TRUE)
#'                                 
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
      if (cluster_object$name == "hierarchical_clustering") {
        stop("No clusters have been generated for your hierarchical tree,
        please extract clusters from the tree before using bioregionalization_metrics()
        See ?hclu_hierarclust or ?cut_tree")
      } else {
        stop(
          "cluster_object does not have the expected type of 'clusters' slot")
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
      stop("All elements in cluster_object should be of data.frame format")
      # ensure that the number of rows are identical among df of the list
    } else if (length(unique(vapply(cluster_object,
                                    nrow,
                                    FUN.VALUE = numeric(1)))) == 1){
      clusters <- data.frame(cluster_object)
    }
    # if none of the above then stop
  }  else {
    stop("This function is designed to work either on bioregion.clusters
    objects, on a data.frame or on a list of data.frames")
  }
  
  if(ncol(clusters) == 1) {
    stop(
      "This function is designed to be applied on multiple bioregionalizations.",
      "Your cluster_object only has a single bioregionalization (one column)")
  }
  
  if(!is.null(indices)) {
    controls(args = indices, type = "character_vector")
    
    if(!isTRUE(unique(indices %in% c("rand", "jaccard")))){
      stop("Please choose algorithm among the followings values:
    rand or jaccard", call. = FALSE)
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
