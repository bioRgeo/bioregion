#' Compare cluster memberships among multiple partitions
#' 
#' This function aims at computing pairwise comparisons for several
#' partitions, usually on outputs from `netclu_`, `hclu_` or `nhclu_` functions.
#' It also provides the confusion matrix from pairwise comparisons, so that 
#' the user can compute additional comparison metrics.
#' 
#' @param cluster_object a `bioregion.clusters` object or a `data.frame` or a 
#' list of `data.frame` containing multiple partitions. At least two partitions
#' are required. If a list of `data.frame` is provided, they should all have
#' the same number of rows (i.e., same items in the clustering for all
#' partitions). 
#' @param sample_comparisons `NULL` or a positive integer. Reduce computation
#' time by sampling a number of pairwise comparisons in cluster membership
#' of items. Useful if the number of items clustered is high. Suggested
#' values 5000 or 10000.
#' @param indices `NULL` or `character`. Indices to compute for the pairwise
#' comparison of partitions. Current available metrics are `"rand"` and 
#' `"jaccard"`
#' @param cor_frequency a boolean. If `TRUE`, then computes the correlation 
#' between each partition and the total frequency of co-membership of items 
#' across all partitions. Useful to identify which partition(s) is(are) most
#' representative of all the computed partitions.
#' @param store_pairwise_membership a boolean. If `TRUE`, the pairwise 
#' membership of items is stored in the output object.
#' @param store_confusion_matrix a boolean. If `TRUE`, the confusion matrices 
#' of pairwise partition comparisons are stored in the output object.
#' 
#' @details 
#' \loadmathjax
#' This function proceeds in two main steps:
#' 1. The first step is done within each partition. It will compare all pairs of
#' items and document if they are clustered together (`TRUE`) or separately 
#' (`FALSE`) in each partition. For example, if site 1 and site 2 are clustered
#' in the same cluster in partition 1, then the pairwise membership site1_site2
#' will be `TRUE`. The output of this first step is stored in the slot
#' `pairwise_membership` if `store_pairwise_membership = TRUE`.
#' 2. The second step compares all pairs of partitions by analysing if their
#' pairwise memberships are similar or not. To do so, for each pair of 
#' partitions, the function computes a confusion matrix with four elements:
#'  * _a_: number of pairs of items grouped in partition 1 and in partition 2
#'  * _b_: number of pairs of items grouped in partition 1 but not in partition
#'  2
#'  * _c_: number of pairs of items not grouped in partition 1 but grouped in
#'  partition 2
#'  * _d_: number of pairs of items not grouped in both partition 1 & 2
#'  
#' The confusion matrix is stored in `confusion_matrix` if 
#' `store_confusion_matrix = TRUE`.
#' 
#' 
#' Based on the confusion matrices, we can compute a range of indices to 
#' indicate the agreement among partitions. As of now, we have implemented:
#' * _Rand index_ 
#' \mjeqn{(a + d)/(a + b + c + d)}{(a + d)/(a + b + c + d)}
#' The Rand index measures agreement among partitions by accounting for both
#' the pairs of sites that are grouped, but also the pairs of sites that are
#' not grouped.
#' * _Jaccard index_
#' \mjeqn{(a)/(a + b + c)}{(a)/(a + b + c)}
#' The Jaccard index measures agreement among partitions by only accounting
#' for pairs of sites that are grouped - it is
#' 
#' These two metrics are complementary, because the Jaccard index will tell
#' if partitions are similar in their clustering structure, whereas the 
#' Rand index will tell if partitions are similar not only in the pairs of
#' items clustered together, but also in terms of the pairs of sites that are
#' not clustered together. For example, take two partitions which
#' never group together the same pairs of sites. Their Jaccard index will be 0,
#' whereas the Rand index can be > 0 due to the sites that are not grouped 
#' together. 
#' 
#' Additional indices can be manually computed by the users on the basis of the
#' list of confusion matrices.
#' 
#' 
#' In some cases, users may be interested in finding which of the partitions
#' is most representative of all partitions. To find it out, we can
#' compare the pairwise membership of each partition with the total frequency
#' of pairwise membership across all partitions. This correlation can be
#' requested with `cor_frequency = TRUE`
#' @seealso [partition_metrics]
#' @return 
#' A `list` with 4 to 7 elements:
#'  * `args`: arguments provided by the user
#'  * `inputs`: information on the input partitions, such as the number of items
#'  being clustered
#'  * (facultative) `pairwise_membership`: only if 
#'  `store_pairwise_membership = TRUE`. This 
#'  element contains the pairwise memberships of all items for each 
#'  partition, in the form of a `boolean matrix` where `TRUE` means that 
#'  two items are in the same cluster, and `FALSE` means that two items 
#'  are not in the same cluster
#'  * `freq_item_pw_membership`: A `numeric vector` 
#'  containing the number of times each pair of items are clustered
#'  together. It corresponds to the sum of rows of the table in 
#'  `pairwise_membership`
#'  * (facultative) `partition_freq_cor`:  only if `cor_frequency = TRUE`.
#'   A `numeric vector`
#'  indicating the correlation between individual partitions and the total
#'  frequency of pairwise membership across all partitions. It corresponds to
#'  the correlation between individual columns in `pairwise_membership` and 
#'  `freq_item_pw_membership`
#'  * (facultative) `confusion_matrix`: only if `store_confusion_matrix = TRUE`.
#'   A `list` 
#'  containing all confusion matrices between each pair of partitions.
#'  * `partition_comparison`: a `data.frame` containing the results of the
#'  comparison of partitions, where the first column indicates which partitions 
#'  are compared, and the next columns correspond to the requested `indices`.
#'  
#' @author
#' Boris Leroy (\email{leroy.boris@gmail.com}),
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) and
#' Pierre Denelle (\email{pierre.denelle@gmail.com})
#' 
#' @examples 
#' # A simple case with four partitions of four items
#' partitions <- data.frame(matrix(nr = 4, nc = 4, 
#'                                 c(1,2,1,1,1,2,2,1,2,1,3,1,2,1,4,2),
#'                                 byrow = TRUE))
#' partitions
#' compare_partitions(partitions)
#' 
#' # Find out which partitions are most representative
#' compare_partitions(partitions,
#'                    cor_frequency = TRUE)
#'                                 
#' 
#' @export

compare_partitions <- function(cluster_object,
                               sample_comparisons = NULL,
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
        please extract clusters from the tree before using partition_metrics()
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
    stop("This function is designed to be applied on multiple partitions.",
         "Your cluster_object only has a single partition (one column)")
  }
  
  if(!is.null(sample_comparisons)) {
    controls(data = sample_comparisons, type = "integer")
  }
  
  if(!is.null(indices)) {
    indices <- controls(args = indices, type = "character_vector")
  }
  
  controls(args = cor_frequency, type = "boolean")
  controls(args = store_pairwise_membership, type = "boolean")
  controls(args = store_confusion_matrix, type = "boolean")
  
  
  message(Sys.time(), 
          " - Computing pairwise membership comparisons for each",
          "partition...\n")
  if(ncol(clusters) * (nrow(clusters) * (nrow(clusters) - 1)) / 2 > 10e6) {
    message("       /!\\\ NOTE: Very high number of comparisons ",
            "(", nrow(clusters), " rows in clusters --> ",
            (nrow(clusters) * (nrow(clusters) - 1)) / 2, 
            " pairwise comparisons per partition --> ",
            ncol(clusters) * (nrow(clusters) * (nrow(clusters) - 1)) / 2, 
            " comparisons in total across all partitions).\n      ",
            "If the computation ",
            "time is too long, consider using argument sample_comparisons and",
            " set it to a reasonable number (e.g. 5000).")
  }
  
  
  # WITHIN PARTITION - Pairwise membership comparison -----------------------
  
  # For each partition, compare all pairs of items to detect if they are in the
  # same cluster or not
  item_pw_mb <- get_pairwise_membership(
    clusters,
    sample_pw_comparisons = sample_comparisons)
  
  # BETWEEN PARTITIONS - Pairwise partition comparison ----------------------
  
  message(Sys.time(), 
          " - Comparing memberships among partitions...\n")
  
  # Prepare the pairwise partition comparisons
  partnames <- pw_cluster_comps <- t(utils::combn(seq_len(ncol(clusters)), 2))
  
  
  # Compute confusion matrices for all pairs of partitions
  all_conf_matrices <- lapply(seq_len(nrow(pw_cluster_comps)), function(x) {
    get_confusion_matrix(item_pw_mb[, pw_cluster_comps[x, 1]],
                         item_pw_mb[, pw_cluster_comps[x, 2]])
  })
  
  # partnames gets the names of each partition 
  partnames[] <- colnames(clusters)[pw_cluster_comps]
  names(all_conf_matrices) <- apply(partnames, 1, paste, collapse = "%")
  
  # Item membership frequency -----------------------------------------------
  
  item_pw_mb_freq <- rowSums(item_pw_mb)
  
  # Prepare data.frame comparing partitions ---------------------------------
  
  partcomp_indices <- data.frame(
    partition_comparison = names(all_conf_matrices))
  
  
  # Rand index --------------------------------------------------------------
  
  if("rand" %in% indices) {
    message(Sys.time(), 
            " - Computing Rand index...\n")
    
    partcomp_indices$rand <- vapply(all_conf_matrices,
                                    rand_index,
                                    numeric(1))
  }
  
  
  # Jaccard index -----------------------------------------------------------
  
  if("jaccard" %in% indices) {
    message(Sys.time(), 
            " - Computing Jaccard index...\n")
    partcomp_indices$jaccard <- vapply(all_conf_matrices,
                                       jaccard_index,
                                       numeric(1))
  }
  
  
  
  # Point biserial correlation  ---------------------------------------------
  
  if(cor_frequency) {
    message(Sys.time(), 
            " - Computing the correlation between each partition and the",
            " vector of frequency of pairwise membership...\n")
    
    partition_freq_cor <- as.vector(suppressWarnings(
      cor(item_pw_mb, item_pw_mb_freq)))
    partition_freq_cor[is.na(partition_freq_cor)] <- 0
    names(partition_freq_cor) <- colnames(clusters)
  }
  
  
  # Store outputs -----------------------------------------------------------
  
  outputs <- list(
    args = list(sample_comparisons = sample_comparisons,
                indices = indices,
                cor_frequency = cor_frequency,
                store_pairwise_membership = store_pairwise_membership,
                store_confusion_matrix = store_confusion_matrix))
  
  outputs$inputs <- c(number_items = nrow(clusters),
                      number_partitions = ncol(clusters))
  
  if(store_pairwise_membership) {
    outputs$pairwise_membership <- item_pw_mb
  }
  
  outputs$freq_item_pw_membership <- item_pw_mb_freq
  
  if(cor_frequency) {
    outputs$partition_freq_cor <-   partition_freq_cor
  }
  
  if(store_confusion_matrix) {
    outputs$confusion_matrix <- all_conf_matrices
  }
  
  outputs$partition_comparison <- partcomp_indices
  
  class(outputs) <- append("bioregion.partition.comparison", class(outputs))
  return(outputs)
}


# Get pairwise membership -------------------------------------------------


get_pairwise_membership <- function(input_clusters,
                                    sample_pw_comparisons = NULL) {
  
  # Get unique combinations of items to compare
  pairwise_comps <- t(utils::combn(seq_len(nrow(input_clusters)), 2))
  
  if(!is.null(sample_pw_comparisons)) {
    pairwise_comps <- pairwise_comps[sample(sample_pw_comparisons), ]
  }
  
  # Vectorised computation of group membership
  pw_membership <- vapply(seq_len(ncol(input_clusters)), function(x) {
    # First vapply to work on each columns
    vapply(seq_len(nrow(pairwise_comps)),
           # Second vapply to work on each pairwise combination of items
           function(j, col){
             # Test if item 1 = item 2
             input_clusters[pairwise_comps[j, 1], col] ==
               input_clusters[pairwise_comps[j, 2], col]
           },
           logical(1), # Expect 1 logical in output for each parwise comparison
           col = x) # Specify on which column of input_clusters we are working
  },
  logical(nrow(pairwise_comps))) # Expect a number of logicals in output equal
  # to the number of pairwise comparisons to do
  
  # return(list(pairwise_comparisons = pairwise_comps,
  #             pairwise_membership = pw_membership))
  
  rownames(pw_membership) <- apply(pairwise_comps, 1, paste, collapse = "_")
  colnames(pw_membership) <- colnames(input_clusters)
  
  return(pw_membership)
}


# Confusion matrix --------------------------------------------------------

get_confusion_matrix <- function(clust1, clust2) {
  return(c(a = length(which(clust1 & clust2)), # TP
           b = length(which(clust1 & !clust2)),
           c = length(which(!clust1 & clust2)),
           d = length(which(!clust1 & !clust2)))) # TN
}


# Rand index --------------------------------------------------------------


rand_index <- function(confusion_matrix) {
  RI <- (confusion_matrix["a"] +
           confusion_matrix["d"]) /
    (confusion_matrix["a"] + confusion_matrix["b"] + confusion_matrix["c"] +
       confusion_matrix["d"])
  names(RI) <- "rand"
  return(RI)
}

# Jaccard index -----------------------------------------------------------


jaccard_index <- function(confusion_matrix) {
  JI <- confusion_matrix["a"] /
    (confusion_matrix["a"] + confusion_matrix["b"] + confusion_matrix["c"])
  names(JI) <- "jaccard"
  return(JI)
}