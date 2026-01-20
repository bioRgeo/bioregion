#' @export
#' @method str bioregion.clusters
str.bioregion.clusters <- function(object, ...) {
  args <- list(...)
  if (is.null(args$max.level)) {
    args$max.level <- 2
  }
  NextMethod("str", object = object, max.level = args$max.level)
}


#' @export
#' @method summary bioregion.clusters
summary.bioregion.clusters <- function(object,
                                       n_bioregionalizations = 3,
                                       n_top_clusters = 10,
                                       ...) {
  cat("\n")
  cat("Summary of clustering results\n")
  cat("=============================\n\n")

  # Algorithm information
  cat("Algorithm: ", object$name, "\n")
  cat("Number of sites: ", object$inputs$nb_sites, "\n")

  # Check if bipartite and get node type info
  is_bipartite <- object$inputs$bipartite
  if (is_bipartite) {
    if (!is.null(attr(object$clusters, "node_type"))) {
      n_sites <- sum(attr(object$clusters, "node_type") == "site")
      n_species <- sum(attr(object$clusters, "node_type") == "species")
      cat("  - Site nodes: ", n_sites, "\n")
      cat("  - Species nodes: ", n_species, "\n")
    }
  }

  # Check if hierarchical
  is_hierarchical <- object$inputs$hierarchical
  
  # Check if clusters have been computed
  if (is.data.frame(object$clusters)) {
    n_bioregionalizations_available <- ncol(object$clusters) - 1
  } else {
    n_bioregionalizations_available <- 0
  }

  cat("Number of bioregionalizations: ", n_bioregionalizations_available, "\n")
  if (is_hierarchical) {
    cat("Hierarchical clustering: Yes\n")
  } else {
    cat("Hierarchical clustering: No\n")
  }
  cat("\n")
  
  # If no clusters have been computed, show message and return
  if (n_bioregionalizations_available == 0) {
    cat("No bioregionalizations have been computed yet.\n")
    if (object$name == "hclu_hierarclust") {
      cat("Use cut_tree() to cut the tree and obtain clusters.\n")
    }
    cat("\n")
    invisible(object)
    return(invisible(object))
  }

  # Limit number of bioregionalizations to display
  n_to_show <- min(n_bioregionalizations, n_bioregionalizations_available)

  # Display summary for each bioregionalization
  for (i in 1:n_to_show) {
    bioregionalization_name <- names(object$clusters)[i + 1]
    bioregionalization_data <- object$clusters[, i + 1]

    cat("Bioregionalization ", i, ": ", bioregionalization_name, "\n", sep = "")
    cat(strrep("-", nchar(paste0("Bioregionalization ", i, ": ", bioregionalization_name))), "\n")

    # Calculate cluster sizes
    cluster_sizes <- table(bioregionalization_data[!is.na(bioregionalization_data)])
    cluster_sizes <- sort(cluster_sizes, decreasing = TRUE)

    n_clusters <- length(cluster_sizes)
    cat("Total clusters: ", n_clusters, "\n")

    # Show top clusters
    n_show <- min(n_top_clusters, n_clusters)
    cat("Top ", n_show, " clusters by size:\n", sep = "")

    for (j in 1:n_show) {
      cluster_id <- names(cluster_sizes)[j]
      cluster_size <- cluster_sizes[j]

      # If bipartite, show breakdown by node type
      if (is_bipartite && !is.null(attr(object$clusters, "node_type"))) {
        node_types <- attr(object$clusters, "node_type")[
          bioregionalization_data == cluster_id & !is.na(bioregionalization_data)
        ]
        n_sites_clust <- sum(node_types == "site")
        n_species_clust <- sum(node_types == "species")
        cat(sprintf(
          "  Cluster %s: %d items (%d sites, %d species)\n",
          cluster_id, cluster_size, n_sites_clust, n_species_clust
        ))
      } else {
        cat(sprintf("  Cluster %s: %d items\n", cluster_id, cluster_size))
      }
    }

    if (n_clusters > n_show) {
      cat("  ... and ", n_clusters - n_show, " more cluster(s)\n", sep = "")
    }
    cat("\n")
  }

  # If hierarchical, show structure
  if (is_hierarchical && n_bioregionalizations_available > 1) {
    cat("Hierarchical structure\n")
    cat("======================\n\n")

    # Limit to first n_to_show levels
    n_levels <- min(n_to_show, n_bioregionalizations_available)

    # Build full hierarchy tree
    # Start with level 1 clusters (top level)
    level1_col <- object$clusters[, 2] # First bioregionalization column
    level1_clusters <- sort(unique(level1_col[!is.na(level1_col)]))
    level1_clusters <- level1_clusters[1:min(
      n_top_clusters,
      length(level1_clusters)
    )]

    # Display hierarchy for each top-level cluster
    for (i in seq_along(level1_clusters)) {
      is_last_top <- (i == length(level1_clusters))
      display_hierarchy(object, level1_clusters[i], 1, "",
        is_last = is_last_top,
        max_level = n_levels,
        n_top_clusters = n_top_clusters
      )
      if (i < length(level1_clusters)) cat("\n")
    }

    # Show message if there are more level 1 clusters
    if (length(level1_clusters) <
      length(unique(level1_col[!is.na(level1_col)]))) {
      cat("\n... and ",
        length(unique(level1_col[!is.na(level1_col)])) - length(level1_clusters),
        " more top-level cluster(s)\n",
        sep = ""
      )
    }
    cat("\n")
  }

  invisible(object)
}


# Recursive function to display hierarchy
display_hierarchy <- function(object, cluster_id, level, prefix = "",
                              is_last = TRUE,
                              max_level, n_top_clusters) {
  if (level > max_level) {
    return()
  }

  # Get current level data
  current_col <- object$clusters[, level + 1]
  cluster_size <- sum(current_col == cluster_id & !is.na(current_col))

  # Display current cluster
  if (level == 1) {
    cat(cluster_id, " (n=", cluster_size, ")\n", sep = "")
  } else {
    connector <- if (is_last) "\u2514\u2500" else "\u251c\u2500"
    cat(prefix, connector, cluster_id, " (n=", cluster_size, ")\n", sep = "")
  }

  # If not at max level, find and display children
  if (level < max_level) {
    next_col <- object$clusters[, level + 2]
    children <- unique(next_col[current_col == cluster_id &
      !is.na(next_col)])
    children <- sort(children)

    if (length(children) > 0) {
      # Limit number of children to show
      n_children_show <- min(n_top_clusters, length(children))

      for (j in 1:n_children_show) {
        child <- children[j]
        # Check if this is the last child
        is_last_child <- (j == n_children_show &&
          length(children) <= n_top_clusters)

        # Update prefix for next level
        # The prefix is what goes BEFORE the connector
        # For level 2 (children of level 1), we want NO prefix before connector
        # For level 3+, we want the parent's prefix + continuation
        if (level == 1) {
          # Children of level 1 start with empty prefix
          # Their own children will get vertical bar or spaces based on is_last_child
          new_prefix <- ""
        } else if (level == 2) {
          # Children of level 2 need to continue the line from their parent
          # Use the parent's prefix + vertical line or space
          new_prefix <- if (is_last) "  " else "\u2502 "
        } else {
          # At deeper levels, inherit parent's full prefix
          new_prefix <- paste0(
            prefix,
            if (is_last) "  " else "\u2502 "
          )
        }

        display_hierarchy(
          object, child, level + 1, new_prefix, is_last_child,
          max_level, n_top_clusters
        )
      }

      # Show truncation message if needed
      if (length(children) > n_children_show) {
        connector <- "\u2514\u2500"
        if (level == 1) {
          cat(connector, "... and ", length(children) - n_children_show,
            " more\n",
            sep = ""
          )
        } else if (level == 2) {
          cat(if (is_last) "  " else "\u2502 ", connector,
            "... and ", length(children) - n_children_show,
            " more\n",
            sep = ""
          )
        } else {
          cat(prefix, if (is_last) "  " else "\u2502 ", connector,
            "... and ", length(children) - n_children_show,
            " more\n",
            sep = ""
          )
        }
      }
    }
  }
}


#' @export
#' @method print bioregion.clusters
print.bioregion.clusters <- function(x, ...) {
  # algorithm name -----
  cat("Clustering results for algorithm : ")
  cat(x$name, "\n")
  if (x$name == "hclu_hierarclust") {
    cat("\t(hierarchical clustering based on a dissimilarity matrix)\n")
  }

  # dataset characteristics -----
  cat(" - Number of sites: ", x$inputs$nb_sites, "\n")

  # methodological details -----
  if (x$name %in% c(
    "hclu_hierarclust",
    "hclu_diana"
  )) {
    cat(
      " - Name of dissimilarity metric: ",
      ifelse(is.null(x$args$index),
        "Undefined",
        x$args$index
      ), "\n"
    )
    if (x$name == "hclu_hierarclust") {
      cat(" - Tree construction method: ", x$args$method, "\n")
      cat(
        " - Randomization of the dissimilarity matrix: ",
        ifelse(x$args$randomize, paste0(
          "yes, number of trials ",
          x$args$n_runs
        ), "no"), "\n"
      )
      cat(
        " - Method to compute the final tree: ",
        ifelse(x$args$optimal_tree_method == "best",
          "Tree with the best cophenetic correlation coefficient",
          ifelse(x$args$optimal_tree_method == "iterative_consensus_tree",
            "Iterative hierarchical consensus tree",
            paste0(
              "Consensus tree with p = ",
              x$args$consensus_p
            )
          )
        ), "\n"
      )
    }
    cat(
      " - Cophenetic correlation coefficient: ",
      round(x$algorithm$final.tree.coph.cor, 3), "\n"
    )
  }


  # number of clusters -----
  if (inherits(x$clusters, "data.frame")) {
    # Further methodological details if hclust
    if (x$name == "hclu_hierarclust") {
      if (!is.null(x$args$n_clust)) {
        cat(
          " - Number of clusters requested by the user: ",
          ifelse(length(x$args$n_clust) > 10,
            paste0(
              paste(x$args$n_clust[1:10], collapse = " "),
              " ... (with ",
              length(x$args$n_clust) - 10, " more values)"
            ),
            x$args$n_clust
          ), "\n"
        )
      }
      if (!is.null(x$args$cut_height)) {
        cat(
          " - Heights of cut requested by the user: ",
          ifelse(length(x$args$cut_height) > 10,
            paste0(
              paste(round(x$args$cut_height, 3)[1:10],
                collapse = " "
              ),
              " ... (with ",
              length(x$args$cut_height) - 10, " more values)"
            ),
            paste(round(x$args$cut_height, 3), collapse = " ")
          ), "\n"
        )
      }
      if (x$args$dynamic_tree_cut) {
        cat(paste0(
          " - Dynamic tree cut method chosen: '", x$args$dynamic_method,
          "', with minimum cluster size ", x$args$dynamic_minClusterSize,
          "\n"
        ))
      }
    }

    cat("Clustering results:\n")
    cat(
      " - Number of partitions: ",
      ncol(x$clusters) - 1, "\n"
    )

    if (ncol(x$clusters) > 2) {
      if (x$input$hierarchical) {
        cat(" - Partitions are hierarchical\n")
      } else {
        cat(" - Partitions are not hierarchical\n")
      }
    }

    nclust <- apply(
      x$clusters[, 2:ncol(x$clusters), drop = FALSE],
      2, function(y) length(unique(y))
    )

    cat(
      " - Number of clusters: ",
      ifelse(length(nclust) > 10,
        paste0(
          paste(nclust[1:10], collapse = " "),
          " ... (with ",
          length(nclust) - 10, " more values)"
        ),
        paste(nclust, collapse = " ")
      ),
      "\n"
    )

    if (x$name == "hclu_hierarclust") {
      if (x$args$find_h) {
        cat(
          " - Height of cut of the hierarchical tree:",
          ifelse(length(x$algorithm$output_cut_height) > 10,
            paste0(
              paste(round(x$algorithm$output_cut_height, 3)[1:10],
                collapse = " "
              ),
              " ... (with ",
              length(x$algorithm$output_cut_height) - 10,
              " more values)"
            ),
            paste(round(x$algorithm$output_cut_height, 3),
              collapse = " "
            )
          ), "\n"
        )
      } else {
        cat(" - Height of cut not searched for.", "\n")
      }
    }
    
    # Display color information if present
    if (!is.null(x$colors)) {
      cat(" - Color palette assigned:\n")
      
      # Process each partition
      for (part_name in names(x$colors)) {
        partition_colors <- x$colors[[part_name]]
        
        # Identify grey colors (RGB values all equal, excluding black)
        is_grey <- sapply(partition_colors$color, function(col) {
          if (col == "#000000") {
            return(FALSE)
          }
          r <- substr(col, 2, 3)
          g <- substr(col, 4, 5)
          b <- substr(col, 6, 7)
          return(r == g && g == b)
        })

        nb_vivid <- sum(partition_colors$color != "#000000" & !is_grey)
        nb_grey <- sum(is_grey)
        nb_insignificant <- sum(partition_colors$color == "#000000")

        cat("   * ", part_name, ": ")
        parts <- character(0)
        if (nb_vivid > 0) parts <- c(parts, paste(nb_vivid, "vivid colors"))
        if (nb_grey > 0) parts <- c(parts, paste(nb_grey, "grey shades"))
        if (nb_insignificant > 0) {
          parts <- c(
            parts,
            paste(nb_insignificant, "insignificant (black)")
          )
        }
        cat(paste(parts, collapse = ", "), "\n")
      }
    }
  } else {
    cat("Clustering procedure incomplete - no clusters yet\n")
  }
}


#' @export
#' @method plot bioregion.clusters
plot.bioregion.clusters <- function(x, ...) {
  if (x$name == ("hclu_hierarclust")) {
    args <- list(...)
    # Changing default arguments for hclust plot
    if (is.null(args$xlab)) {
      args$xlab <- ""
    }
    if (is.null(args$ylab)) {
      args$ylab <- paste0(x$args$index, " dissimilarity")
    }
    if (is.null(args$main)) {
      args$main <- ""
    }
    if (is.null(args$sub)) {
      args$sub <- ""
    }
    if (is.null(args$hang)) {
      args$hang <- -1
    }
    args$x <- x$algorithm$final.tree

    do.call(
      plot,
      args
    )
    if (!is.null(x$algorithm$output_cut_height)) {
      # abline(h = x$output_cut_height, lty = 3, col = "#756bb1")

      if (length(x$algorithm$output_cut_height) > 1) {
        if (length(x$algorithm$output_cut_height) > 3) {
          message(
            "Multiple cuts detected, plotting only the first three levels"
          )
        }

        cols <- c("#253494", "#2c7fb8", "#41b6c4")

        for (i in 1:min(3, length(x$algorithm$output_cut_height)))
        {
          stats::rect.hclust(x$algorithm$final.tree,
            h = x$algorithm$output_cut_height[i],
            border = cols[i]
          )
        }
      } else {
        stats::rect.hclust(x$algorithm$final.tree,
          h = x$algorithm$output_cut_height,
          border = "#377eb8"
        )
      }
    } else if (x$args$dynamic_tree_cut) {
      # Adding rectangles for dynamic tree cut
      vect_clust <- x$clusters[, 2]
      names(vect_clust) <- x$clusters[, 1]
      tot_l <- x$algorithm$output_n_clust + length(which(is.na(vect_clust)))

      vect_clust[is.na(vect_clust)] <- (x$algorithm$output_n_clust + 1):
      (x$algorithm$output_n_clust + length(which(is.na(vect_clust))))

      order_rect <- unique(vect_clust[x$algorithm$final.tree$order])

      true_cl <- which(order_rect %in% 1:x$algorithm$output_n_clust)

      stats::rect.hclust(x$final.tree,
        k = tot_l,
        which = true_cl,
        cluster = vect_clust,
        # to do: add border colours from a vector with a
        # distinct colour for each cluster
        border = "#377eb8"
      )
    }
  } else if (x$name == ("hclu_diana")) {
    args <- list(...)
    # Changing default arguments for hclust plot
    if (is.null(args$xlab)) {
      args$xlab <- ""
    }
    if (is.null(args$ylab)) {
      args$ylab <- paste0(x$args$index, " dissimilarity")
    }
    if (is.null(args$main)) {
      args$main <- ""
    }
    if (is.null(args$sub)) {
      args$sub <- ""
    }
    if (is.null(args$ask)) {
      args$ask <- FALSE
    }
    if (is.null(args$which.plots)) {
      args$which.plots <- 2
    }


    args$x <- x$algorithm$final.tree

    do.call(
      plot,
      args
    )
    if (!is.null(x$algorithm$output_cut_height)) {
      # abline(h = x$output_cut_height, lty = 3, col = "#756bb1")

      if (length(x$algorithm$output_cut_height) > 1) {
        if (length(x$algorithm$output_cut_height) > 3) {
          message(
            "Multiple cuts detected, plotting only the first three levels"
          )
        }

        cols <- c("#253494", "#2c7fb8", "#41b6c4")

        for (i in 1:min(3, length(x$algorithm$output_cut_height)))
        {
          stats::rect.hclust(x$algorithm$final.tree,
            h = x$algorithm$output_cut_height[i],
            border = cols[i]
          )
        }
      } else {
        stats::rect.hclust(x$algorithm$final.tree,
          h = x$algorithm$output_cut_height,
          border = "#377eb8"
        )
      }
    } else if (x$args$dynamic_tree_cut) {
      # Adding rectangles for dynamic tree cut
      vect_clust <- x$clusters[, 2]
      names(vect_clust) <- x$clusters[, 1]
      tot_l <- x$algorithm$output_n_clust + length(which(is.na(vect_clust)))

      vect_clust[is.na(vect_clust)] <- (x$algorithm$output_n_clust + 1):
      (x$algorithm$output_n_clust + length(which(is.na(vect_clust))))

      order_rect <- unique(vect_clust[x$algorithm$final.tree$order])

      true_cl <- which(order_rect %in% 1:x$algorithm$output_n_clust)

      stats::rect.hclust(x$final.tree,
        k = tot_l,
        which = true_cl,
        cluster = vect_clust,
        # to do: add border colours from a vector with a
        # distinct colour for each cluster
        border = "#377eb8"
      )
    }
  } else {
    stop("No plot method for this type of object")
  }
}

#' @export
#' @method print bioregion.partition.comparison
print.bioregion.partition.comparison <- function(x, ...) {
  cat("Partition comparison:\n")
  cat(" -", x$inputs["number_partitions"], "partitions compared\n")
  cat(" -", x$inputs["number_items"], "items in the clustering\n")

  if (!is.null(x$args$sample_lines)) {
    cat(
      " - ", x$args$sample_lines,
      " items used to compute comparisons among partitions\n"
    )
  }

  if (!is.null(x$args$indices)) {
    cat(" - Requested indices: ", x$args$indices, "\n")
    cat(" - Metric summary:\n")



    print(data.frame(
      sapply(
        x$partition_comparison[, x$args$indices],
        function(x) {
          c(
            min(x, na.rm = TRUE),
            mean(x, na.rm = TRUE),
            max(x, na.rm = TRUE)
          )
        }
      ),
      row.names = c("Min", "Mean", "Max")
    ))
  } else {
    cat(" - No metrics computed\n")
  }

  if (x$args$cor_frequency) {
    cat(
      " - Correlation between each partition and the total frequency of item",
      " pairwise membership computed:\n"
    )
    cat(
      "   # Range: ", round(min(x$partition_freq_cor), 3), " - ",
      round(max(x$partition_freq_cor), 3), "\n"
    )
    cat(
      "   # Partition(s) most representative (i.e., highest correlation): \n",
      paste(names(x$partition_freq_cor)[
        which(x$partition_freq_cor == max(x$partition_freq_cor))
      ], collapse = ", "),
      "\n Correlation = ", round(max(x$partition_freq_cor), 3), "\n"
    )
  }

  cat(
    " - Item pairwise membership", ifelse(x$args$store_pairwise_membership,
      "", "not"
    ),
    "stored in outputs\n"
  )
  cat(
    " - Confusion matrices of partition comparisons",
    ifelse(x$args$store_confusion_matrix,
      "", "not"
    ),
    "stored in outputs\n"
  )
}

#' @export
#' @method print bioregion.bioregionalization.metrics
print.bioregion.bioregionalization.metrics <- function(x, ...) {
  cat("Partition metrics:\n")
  cat(" -", nrow(x$evaluation_df), " partition(s) evaluated\n")
  cat(
    " - Range of clusters explored: from ", min(x$evaluation_df$n_clusters),
    " to ",
    max(x$evaluation_df$n_clusters), "\n"
  )
  cat(" - Requested metric(s): ", x$args$eval_metric, "\n")
  cat(" - Metric summary:\n")

  print(data.frame(
    sapply(
      x$evaluation_df[x$args$eval_metric],
      function(x) {
        c(
          min(x, na.rm = TRUE),
          mean(x, na.rm = TRUE),
          max(x, na.rm = TRUE)
        )
      }
    ),
    row.names = c("Min", "Mean", "Max")
  ))

  cat("\nAccess the data.frame of metrics with your_object$evaluation_df\n")
  if ("endemism_results" %in% names(x)) {
    cat("Details of endemism % for each bioregionalization are available in
        your_object$endemism_results\n")
  }
}

#' @export
#' @method print bioregion.optimal.n
print.bioregion.optimal.n <- function(x, ...) {
  cat("Search for an optimal number of clusters:\n")
  cat(" -", nrow(x$evaluation_df), " partition(s) evaluated\n")
  cat(
    " - Range of clusters explored: from ", min(x$evaluation_df$n_clusters),
    " to ",
    max(x$evaluation_df$n_clusters), "\n"
  )
  cat(" - Evaluated metric(s): ", x$args$metrics_to_use, "\n")

  cat("\nPotential optimal partition(s):\n")
  cat(
    " - Criterion chosen to optimise the number of clusters: ",
    x$args$criterion, "\n"
  )
  if (x$args$criterion %in% c("increasing_step", "decreasing_step")) ##
    {
      cat(
        "   (step quantile chosen: ", x$args$step_quantile,
        " (i.e., only the top", (1 - x$args$step_quantile) * 100,
        "% ",
        ifelse(x$args$criterion == "increasing_step", "increase", "decrease"),
        " in evaluation metrics",
        " are used as break points for the number of clusters)\n"
      )
    } else if (x$args$criterion == "cutoff") {
    cat("   --> cutoff(s) chosen: ", x$args$metric_cutoffs, "\n")
  }
  cat(" - Optimal partition(s) of clusters for each metric:\n")

  cat(paste(
    paste(names(x$optimal_nb_clusters),
      sapply(x$optimal_nb_clusters,
        paste,
        collapse = " "
      ),
      sep = " - "
    ),
    collapse = "\n"
  ))
  cat("\n")
}

#' @export
#' @method str bioregion.optimal.n
str.bioregion.optimal.n <- function(object, ...) {
  args <- list(...)
  if (is.null(args$max.level)) {
    args$max.level <- 2
  }
  NextMethod("str", object = object, max.level = args$max.level)
}


#' @export
#' @method print bioregion.pairwise
print.bioregion.pairwise <- function(x, ...) {
  metrics <- colnames(x)[-which(colnames(x) %in%
    c(
      "Site1", "Site2", "a", "b",
      "c", "A", "B", "C"
    ))]
  cat(paste0(
    "Data.frame of ",
    ifelse(attr(x, "type") == "similarity",
      "similarity",
      "dissimilarity"
    ),
    " between sites\n"
  ))
  cat(" - Total number of sites: ", attr(x, "nb_sites"), "\n")
  cat(" - Total number of species: ", attr(x, "nb_species"), "\n")
  cat(
    " - Number of rows: ",
    (attr(x, "nb_sites") * (attr(x, "nb_sites") - 1)) / 2, "\n"
  )
  # Warning, next line can be wrong if users alter the object
  cat(
    " - Number of", ifelse(attr(x, "type") == "similarity",
      "similarity",
      "dissimilarity"
    ), "metrics: ",
    length(metrics), "\n"
  )
  cat("\n\n")
  print(as.data.frame(x))
}

#' @export
#' @method `[` bioregion.pairwise
`[.bioregion.pairwise` <- function(x, i, j, ..., drop = TRUE) {
  metric_type <- attributes(x)$type
  nb_sites <- attributes(x)$nb_sites
  nb_species <- attributes(x)$nb_species

  class(x) <- "data.frame"
  out <- x[i, j, ..., drop = drop]
  # We keep track of pw metric class & attribute only if the subset is not a vector
  if (inherits(out, "data.frame")) {
    # if(class(out) == "data.frame") {
    class(out) <- append("bioregion.pairwise", class(out))
    attributes(out)$type <- metric_type
    attributes(out)$nb_sites <- nb_sites
    attributes(out)$nb_species <- nb_species
  }
  out
}

#' @export
#' @method print bioregion.site.species.metrics
print.bioregion.site.species.metrics <- function(x, n_preview = 3, ...) {
  cat("Site and species metrics\n")
  cat("========================")
  cat("\n\n")
  
  # Input summary
  cat("Settings:\n")
  n_part <- attr(x, "n_partitions")
  cat(" - Number of partitions:", n_part, "\n")
  cluster_on <- attr(x, "cluster_on")
  cat(" - Clusters based on:", cluster_on, "\n")
  
  clust_dt <- attr(x, "clustering_data_type")
  if(!is.null(clust_dt) && !is.na(clust_dt)) {
    cat(" - Clustering data type:", clust_dt, "\n")
  }
  
  idx_dt <- attr(x, "index_data_type")
  if(!is.null(idx_dt) && !is.na(idx_dt)) {
    cat(" - Metric data type:", idx_dt, "\n")
  }
  cat("\n")
  
  # Computed metrics
  cat("Computed metrics:\n")
  bio_occ <- attr(x, "bioregion_metrics_occ")
  bio_abd <- attr(x, "bioregion_metrics_abd")
  bioreg_occ <- attr(x, "bioregionalization_metrics_occ")
  bioreg_abd <- attr(x, "bioregionalization_metrics_abd")
  sim_idx <- attr(x, "similarity_metrics")
  
  if(length(bio_occ) > 0)
    cat(" - Per-cluster metrics (occurrence):", paste(bio_occ, collapse = ", "), "\n")
  if(length(bio_abd) > 0)
    cat(" - Per-cluster metrics (abundance):", paste(bio_abd, collapse = ", "), "\n")
  if(length(bioreg_occ) > 0)
    cat(" - Summary metrics (occurrence):", paste(bioreg_occ, collapse = ", "), "\n")
  if(length(bioreg_abd) > 0)
    cat(" - Summary metrics (abundance):", paste(bioreg_abd, collapse = ", "), "\n")
  if(length(sim_idx) > 0)
    cat(" - Similarity-based metrics:", paste(sim_idx, collapse = ", "), "\n")
  cat("\n")
  
  print_df_preview <- function(df, name, n_rows) {
    cat("$", name, " (", nrow(df), " rows x ", ncol(df), " cols):\n", sep = "")
    
    df_preview <- utils::head(df, n_rows)
    num_cols <- sapply(df_preview, is.numeric)
    df_preview[num_cols] <- lapply(df_preview[num_cols], round, 3)
    print(df_preview, row.names = FALSE)
    if(nrow(df) > n_rows) cat("# ... with", nrow(df) - n_rows, "more rows\n")
    cat("\n")
  }
  
  if(n_preview > 0) {
    # If data preview, show first few rows of the first partition
    if(n_part == 1) {
      cat("Data preview:\n")
      for(comp in names(x)) {
        if(is.data.frame(x[[comp]])) {
          print_df_preview(x[[comp]], comp, n_preview)
        }
      }
    } else {
      cat("Data preview (first partition: ", names(x)[1], "):\n", sep = "")
      first_part <- x[[1]]
      for(comp in names(first_part)) {
        if(is.data.frame(first_part[[comp]])) {
          print_df_preview(first_part[[comp]], comp, n_preview)
        }
      }
      cat("Partitions:", paste(names(x), collapse = ", "), "\n\n")
    }
  } else {
    # if not data preview: show dimensions of tables
    cat("Available components:\n")
    if(n_part == 1) {
      for(comp in names(x)) {
        if(is.data.frame(x[[comp]])) {
          cat(" -", comp, ":", nrow(x[[comp]]), "rows x", ncol(x[[comp]]), "columns\n")
        }
      }
    } else {
      first_part <- x[[1]]
      for(comp in names(first_part)) {
        if(is.data.frame(first_part[[comp]])) {
          cat(" -", comp, ":", nrow(first_part[[comp]]), "rows x", 
              ncol(first_part[[comp]]), "columns (per partition)\n")
        }
      }
      cat("\nPartitions:", paste(names(x), collapse = ", "), "\n")
    }
    cat("\n")
  }
  
  cat("Access data with:\n")
  if(n_part == 1) {
    components <- names(x)
    components <- components[sapply(x, is.data.frame)]
    for(comp in components) {
      cat("  your_object$", comp, "\n", sep = "")
    }
  } else {
    part_name <- names(x)[1]
    components <- names(x[[1]])
    components <- components[sapply(x[[1]], is.data.frame)]
    for(comp in components) {
      cat("  your_object$", part_name, "$", comp, "\n", sep = "")
    }
  }
  
  invisible(x)
}

#' @export
#' @method str bioregion.site.species.metrics
str.bioregion.site.species.metrics <- function(object, ...) {
  cat("bioregion.site.species.metrics object\n")
  cat(" - Partitions:", attr(object, "n_partitions"), "\n")
  cat(" - Cluster based on:", attr(object, "cluster_on"), "\n")
  
  clust_dt <- attr(object, "clustering_data_type")
  if(!is.null(clust_dt) && !is.na(clust_dt))
    cat(" - Clustering data type:", clust_dt, "\n")
    
  idx_dt <- attr(object, "index_data_type")
  if(!is.null(idx_dt) && !is.na(idx_dt))
    cat(" - Metric data type:", idx_dt, "\n")
  
  bio_occ <- attr(object, "bioregion_metrics_occ")
  bio_abd <- attr(object, "bioregion_metrics_abd")
  bioreg_occ <- attr(object, "bioregionalization_metrics_occ")
  bioreg_abd <- attr(object, "bioregionalization_metrics_abd")
  sim_idx <- attr(object, "similarity_metrics")
  
  if(length(bio_occ) > 0)
    cat(" - Per-cluster metrics (occurrence):", paste(bio_occ, collapse = ", "), "\n")
  if(length(bio_abd) > 0)
    cat(" - Per-cluster metrics (abundance):", paste(bio_abd, collapse = ", "), "\n")
  if(length(bioreg_occ) > 0)
    cat(" - Summary metrics (occurrence):", paste(bioreg_occ, collapse = ", "), "\n")
  if(length(bioreg_abd) > 0)
    cat(" - Summary metrics (abundance):", paste(bioreg_abd, collapse = ", "), "\n")
  if(length(sim_idx) > 0)
    cat(" - Similarity-based metrics:", paste(sim_idx, collapse = ", "), "\n")
  cat("\n")
  
  args <- list(...)
  if(is.null(args$max.level)) args$max.level <- 2
  NextMethod("str", object = object, max.level = args$max.level)
}

# Helper function to summarize metric columns
.summarize_metrics_df <- function(df, exclude_cols = c("Species", "Site", 
                                                        "Bioregion", 
                                                        "Chorotype",
                                                        "Chorotypes",
                                                        "Assigned")) {
  metric_cols <- setdiff(names(df), exclude_cols)
  metric_cols <- metric_cols[sapply(df[metric_cols], is.numeric)]
  if(length(metric_cols) == 0) return(NULL)
  
  # Helper to safely compute stats, returning NA if no valid values
  safe_min <- function(x) {
    x <- x[!is.na(x)]
    if(length(x) == 0) NA_real_ else min(x)
  }
  safe_max <- function(x) {
    x <- x[!is.na(x)]
    if(length(x) == 0) NA_real_ else max(x)
  }
  safe_mean <- function(x) {
    x <- x[!is.na(x)]
    if(length(x) == 0) NA_real_ else mean(x)
  }
  
  stats <- data.frame(
    Metric = metric_cols,
    Min = sapply(df[metric_cols], safe_min),
    Mean = sapply(df[metric_cols], safe_mean),
    Max = sapply(df[metric_cols], safe_max),
    row.names = NULL
  )
  stats
}

#' @export
#' @method summary bioregion.site.species.metrics
summary.bioregion.site.species.metrics <- function(object,
                                                   n_partitions = 3,
                                                   n_top = 5,
                                                   show_top_contributors = TRUE,
                                                   ...) {
  cat("\nSummary of site and species metrics\n")
  cat("===================================\n\n")
  
  # Settings
  cat("Settings:\n")
  n_part <- attr(object, "n_partitions")
  cat(" - Number of partitions:", n_part, "\n")
  cluster_on <- attr(object, "cluster_on")
  cat(" - Cluster based on:", cluster_on, "\n")
  
  clust_dt <- attr(object, "clustering_data_type")
  if(!is.null(clust_dt) && !is.na(clust_dt))
    cat(" - Clustering data type:", clust_dt, "\n")
    
  idx_dt <- attr(object, "index_data_type")
  if(!is.null(idx_dt) && !is.na(idx_dt))
    cat(" - Metric data type:", idx_dt, "\n")
  cat("\n")
  
  n_to_show <- min(n_partitions, n_part)
  
  # Process each partition
  for(p in seq_len(n_to_show)) {
    # Get partition data
    if(n_part == 1) {
      part_data <- object
      part_name <- "Single partition"
    } else {
      part_data <- object[[p]]
      part_name <- names(object)[p]
    }
    
    # Count clusters if available
    n_bio <- NA
    n_choro <- NA
    if(!is.null(part_data$species_bioregions)) {
      n_bio <- length(unique(part_data$species_bioregions$Bioregion))
    } else if(!is.null(part_data$site_bioregions)) {
      n_bio <- length(unique(part_data$site_bioregions$Bioregion))
    }
    if(!is.null(part_data$site_chorotypes)) {
      n_choro <- length(unique(part_data$site_chorotypes$Chorotypes))
    }
    
    clust_str <- ""
    if(!is.na(n_bio) && !is.na(n_choro)) {
      clust_str <- paste0(" (", n_bio, " bioregions, ", n_choro, " chorotypes)")
    } else if(!is.na(n_bio)) {
      clust_str <- paste0(" (", n_bio, " bioregions)")
    } else if(!is.na(n_choro)) {
      clust_str <- paste0(" (", n_choro, " chorotypes)")
    }
    cat("Partition ", p, ": ", part_name, clust_str, "\n", sep = "")
    cat(strrep("-", nchar(paste0("Partition ", p, ": ", part_name, clust_str))), "\n\n")
    
    # Species-per-bioregion metrics
    if(!is.null(part_data$species_bioregions)) {
      cat("Species-per-bioregion metrics ($species_bioregions):\n")
      stats <- .summarize_metrics_df(part_data$species_bioregions)
      if(!is.null(stats)) {
        stats$Min <- round(stats$Min, 3)
        stats$Mean <- round(stats$Mean, 3)
        stats$Max <- round(stats$Max, 3)
        print(stats, row.names = FALSE)
      }
      cat("\n")
      
      # Top contributors by IndVal
      if(show_top_contributors) {
        indval_cols <- grep("IndVal", names(part_data$species_bioregions), 
                           value = TRUE)
        if(length(indval_cols) > 0) {
          indval_col <- indval_cols[1]
          df <- part_data$species_bioregions
          df <- df[order(-df[[indval_col]]), ]
          n_show <- min(n_top, nrow(df))
          cat("Top species by ", indval_col, ":\n", sep = "")
          for(i in seq_len(n_show)) {
            cat("  ", i, ". ", df$Species[i], " (Bioregion ", df$Bioregion[i], 
                "): ", round(df[[indval_col]][i], 3), "\n", sep = "")
          }
          cat("\n")
        }
      }
    }
    
    # Species summary metrics
    if(!is.null(part_data$species_bioregionalization)) {
      cat("Species summary metrics ($species_bioregionalization):\n")
      stats <- .summarize_metrics_df(part_data$species_bioregionalization)
      if(!is.null(stats)) {
        stats$Min <- round(stats$Min, 3)
        stats$Mean <- round(stats$Mean, 3)
        stats$Max <- round(stats$Max, 3)
        print(stats, row.names = FALSE)
      }
      cat("\n")
    }
    
    # Site-per-chorotype metrics
    if(!is.null(part_data$site_chorotypes)) {
      cat("Site-per-chorotype metrics ($site_chorotypes):\n")
      stats <- .summarize_metrics_df(part_data$site_chorotypes)
      if(!is.null(stats)) {
        stats$Min <- round(stats$Min, 3)
        stats$Mean <- round(stats$Mean, 3)
        stats$Max <- round(stats$Max, 3)
        print(stats, row.names = FALSE)
      }
      cat("\n")
      
      # Top contributors by IndVal for sites
      if(show_top_contributors) {
        indval_cols <- grep("IndVal", names(part_data$site_chorotypes), 
                           value = TRUE)
        if(length(indval_cols) > 0) {
          indval_col <- indval_cols[1]
          df <- part_data$site_chorotypes
          df <- df[order(-df[[indval_col]]), ]
          n_show <- min(n_top, nrow(df))
          cat("Top sites by ", indval_col, ":\n", sep = "")
          for(i in seq_len(n_show)) {
            cat("  ", i, ". ", df$Site[i], " (Chorotype ", df$Chorotypes[i], 
                "): ", round(df[[indval_col]][i], 3), "\n", sep = "")
          }
          cat("\n")
        }
      }
    }
    
    # Site chorological summary metrics
    if(!is.null(part_data$site_chorological)) {
      cat("Site chorological summary metrics ($site_chorological):\n")
      stats <- .summarize_metrics_df(part_data$site_chorological)
      if(!is.null(stats)) {
        stats$Min <- round(stats$Min, 3)
        stats$Mean <- round(stats$Mean, 3)
        stats$Max <- round(stats$Max, 3)
        print(stats, row.names = FALSE)
      }
      cat("\n")
    }
    
    # Site-per-bioregion metrics (diversity & similarity-based)
    if(!is.null(part_data$site_bioregions)) {
      cat("Site-per-bioregion metrics ($site_bioregions):\n")
      stats <- .summarize_metrics_df(part_data$site_bioregions)
      if(!is.null(stats)) {
        stats$Min <- round(stats$Min, 3)
        stats$Mean <- round(stats$Mean, 3)
        stats$Max <- round(stats$Max, 3)
        print(stats, row.names = FALSE)
      }
      cat("\n")
    }
    
    # Site summary metrics (Silhouette)
    if(!is.null(part_data$site_bioregionalization)) {
      cat("Site summary metrics ($site_bioregionalization):\n")
      stats <- .summarize_metrics_df(part_data$site_bioregionalization)
      if(!is.null(stats)) {
        stats$Min <- round(stats$Min, 3)
        stats$Mean <- round(stats$Mean, 3)
        stats$Max <- round(stats$Max, 3)
        print(stats, row.names = FALSE)
      }
      cat("\n")
    }
  }
  
  # Message if truncated
  if(n_part > n_to_show) {
    cat("... and ", n_part - n_to_show, " more partition(s).\n\n", sep = "")
  }
  
  invisible(object)
}

#' @export
#' @method as.dist bioregion.pairwise
as.dist.bioregion.pairwise <- function(m, diag = FALSE,
                                       upper = FALSE) {
  if (ncol(x) > 3) {
    message(
      "More than 3 columns in x: using the third column as the distance",
      "index"
    )
    x <- x[, 1:3]
  }
  matrix.dist <- net_to_mat(x,
    weight = TRUE, squared = TRUE, symmetrical = TRUE
  )
  matrix.dist <- stats::as.dist(x,
    diag = diag,
    upper = upper
  )
  return(matrix.dist)
}
