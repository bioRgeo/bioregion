#' @export
#' @method str bioregion.clusters
str.bioregion.clusters <- function(object, ...)
{
  args <- list(...)
  if(is.null(args$max.level))
  {
    args$max.level <- 2
  }
  NextMethod("str", object = object, max.level = args$max.level)
}


#' @export
#' @method print bioregion.clusters
print.bioregion.clusters <- function(x, ...)
{
  # algorithm name -----
  cat("Clustering results for algorithm : ")
  cat(x$name, "\n")
  if(x$name == "hierarchical_clustering") {
    cat("\t(hierarchical clustering based on a dissimilarity matrix)\n")
  }
  
  # dataset characteristics -----
  cat(" - Number of sites: ", x$inputs$nb_sites, "\n")
  
  # methodological details -----
  if(x$name == "hierarchical_clustering") {
    cat(" - Name of dissimilarity metric: ",
        ifelse(is.null(x$args$index),
               "Undefined",
               x$args$index), "\n")
    cat(" - Tree construction method: ", x$args$method, "\n")
    cat(" - Randomization of the dissimilarity matrix: ",
        ifelse(x$args$randomize, paste0("yes, number of trials ",
                                        x$args$n_runs), "no"), "\n")
    cat(" - Cophenetic correlation coefficient: ",
        round(x$algorithm$final.tree.coph.cor, 3), "\n")
  }
  

  # number of clusters -----
  if (inherits(x$clusters, "data.frame")) {
    
    # Further methodological details if hclust
    if(x$name == "hierarchical_clustering") {
      if(!is.null(x$args$n_clust))
      {
        cat(" - Number of clusters requested by the user: ",
            ifelse(length(x$args$n_clust) > 10,
                   paste0(paste(x$args$n_clust[1:10], collapse = " "),
                          " ... (with ",
                          length(x$args$n_clust) - 10, " more values)"),
                   x$args$n_clust), "\n")
      }
      if(!is.null(x$args$cut_height))
      {
        cat(" - Heights of cut requested by the user: ",
            ifelse(length(x$args$cut_height) > 10,
                   paste0(paste(round(x$args$cut_height, 3)[1:10],
                                collapse = " "),
                          " ... (with ",
                          length(x$args$cut_height) - 10, " more values)"),
                   paste(round(x$args$cut_height, 3), collapse = " ")), "\n")
      }
      if(x$args$dynamic_tree_cut)
      {
        cat(paste0(
          " - Dynamic tree cut method chosen: '", x$args$dynamic_method,
          "', with minimum cluster size ", x$args$dynamic_minClusterSize,
          "\n"))
      }
      
    }
    
    cat("Clustering results:\n")
    cat(" - Number of partitions: ",
        ncol(x$clusters) - 1, "\n")
    
    if(ncol(x$clusters) > 2) {
      if(x$input$hierarchical) {
        cat(" - Partitions are hierarchical\n")
      } else {
        cat(" - Partitions are not hierarchical\n")
      }
    }
    
    nclust <- apply(x$clusters[, 2:ncol(x$clusters), drop = FALSE],
                    2, function(y) length(unique(y)))
    
    cat(" - Number of clusters: ",
        ifelse(length(nclust) > 10,
               paste0(paste(nclust[1:10], collapse = " "),
                      " ... (with ",
                      length(nclust) - 10, " more values)"),
               paste(nclust, collapse = " ")),
        "\n")
    
    if(x$name == "hierarchical_clustering") {
      if(x$args$find_h)
      {
        cat(" - Height of cut of the hierarchical tree:",
            ifelse(length(x$algorithm$output_cut_height) > 10,
                   paste0(paste(round(x$algorithm$output_cut_height, 3)[1:10],
                                collapse = " "),
                          " ... (with ",
                          length(x$algorithm$output_cut_height) - 10,
                          " more values)"),
                   paste(round(x$algorithm$output_cut_height, 3),
                         collapse = " ")), "\n")
      } else
      {
        cat(" - Height of cut not searched for.", "\n")
      }
    }
  } else {
    cat("Clustering procedure incomplete - no clusters yet\n")
  }
}


#' @export
#' @method plot bioregion.clusters
plot.bioregion.clusters <- function(x, ...)
{
  if(x$name == ("hierarchical_clustering"))
  {
    args <- list(...)
    # Changing default arguments for hclust plot
    if(is.null(args$xlab))
    {
      args$xlab <- ""
    }
    if(is.null(args$ylab))
    {
      args$ylab <- paste0(x$args$index, " dissimilarity")
    }
    if(is.null(args$main))
    {
      args$main <- ""
    }
    if(is.null(args$sub))
    {
      args$sub <- ""
    }
    if(is.null(args$hang))
    {
      args$hang <- -1
    }
    args$x <- x$algorithm$final.tree
    
    do.call(plot,
            args)
    if(!is.null(x$algorithm$output_cut_height))
    {
      # abline(h = x$output_cut_height, lty = 3, col = "#756bb1")
      
      if(length(x$algorithm$output_cut_height) > 1)
      {
        if(length(x$algorithm$output_cut_height) > 3)
        {
          message(
            "Multiple cuts detected, plotting only the first three levels")
        }
        
        cols <- c("#253494", "#2c7fb8", "#41b6c4")
        
        for(i in 1:min(3, length(x$algorithm$output_cut_height)))
        {
          stats::rect.hclust(x$algorithm$final.tree,
                             h = x$algorithm$output_cut_height[i],
                             border = cols[i])
        }
        
      } else
      {
        stats::rect.hclust(x$algorithm$final.tree,
                           h = x$algorithm$output_cut_height,
                           border = "#377eb8")
      }
    } else if(x$args$dynamic_tree_cut)
    {
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
                         border = "#377eb8")
    }
  } else
  {
    stop("No plot method for this type of object")
  }
}

#' @export
#' @method print bioregion.partition.comparison
print.bioregion.partition.comparison <- function(x, ...)
{
  cat("Partition comparison:\n")
  cat(" -", x$inputs["number_partitions"], "partitions compared\n")
  cat(" -", x$inputs["number_items"], "items in the clustering\n")
  
  if(!is.null(x$args$sample_comparisons)) {
    cat(" - ", x$args$sample_comparisons, 
        "pairwise item comparisons sampled\n")
  }
  
  if(!is.null(x$args$indices)) {
    cat(" - Requested indices: ", x$args$indices, "\n")
    cat(" - Metric summary:\n")

    
    
    print(data.frame(sapply(x$partition_comparison[, x$args$indices],
                            function(x) {
                              c(min(x, na.rm = TRUE), 
                                mean(x, na.rm = TRUE), 
                                max(x, na.rm = TRUE))}),
                     row.names = c("Min", "Mean", "Max")))
  } else {
    cat(" - No metrics computed\n")
  }
  
  if(x$args$cor_frequency) {
    cat(" - Correlation between each partition and the total frequency of item",
        " pairwise membership computed:\n")
    cat("   # Range: ", round(min(x$partition_freq_cor), 3), " - ", 
        round(max(x$partition_freq_cor), 3), "\n")
    cat("   # Partition(s) most representative (i.e., highest correlation): \n", 
        paste(names(x$partition_freq_cor)[
          which(x$partition_freq_cor == max(x$partition_freq_cor))
        ], collapse = ", "),
        "\n Correlation = ", round(max(x$partition_freq_cor), 3), "\n")
  }
 
  cat(" - Item pairwise membership", ifelse(x$args$store_pairwise_membership,
                                             "", "not"),
      "stored in outputs\n")
  cat(" - Confusion matrices of partition comparisons",
      ifelse(x$args$store_confusion_matrix,
             "", "not"),
      "stored in outputs\n")
}

#' @export
#' @method print bioregion.partition.metrics
print.bioregion.partition.metrics <- function(x, ...)
{
  cat("Partition metrics:\n")
  cat(" -", nrow(x$evaluation_df), " partition(s) evaluated\n")
  cat(" - Range of clusters explored: from ", min(x$evaluation_df$n_clusters),
      " to ",
      max(x$evaluation_df$n_clusters), "\n")
  cat(" - Requested metric(s): ", x$args$eval_metric, "\n")
  cat(" - Metric summary:\n")
  
  print(data.frame(sapply(x$evaluation_df[x$args$eval_metric],
                          function(x) {
                            c(min(x, na.rm = TRUE), 
                              mean(x, na.rm = TRUE), 
                              max(x, na.rm = TRUE))}),
                   row.names = c("Min", "Mean", "Max")))
  
  cat("\nAccess the data.frame of metrics with your_object$evaluation_df\n")
  if("endemism_results" %in% names(x)) {
    cat("Details of endemism % for each partition are available in 
        your_object$endemism_results\n")
  }
}

#' @export
#' @method print bioregion.optimal.n
print.bioregion.optimal.n <- function(x, ...)
{
  cat("Search for an optimal number of clusters:\n")
  cat(" -", nrow(x$evaluation_df), " partition(s) evaluated\n")
  cat(" - Range of clusters explored: from ", min(x$evaluation_df$n_clusters),
      " to ",
      max(x$evaluation_df$n_clusters), "\n")
  cat(" - Evaluated metric(s): ", x$args$metrics_to_use, "\n")
  
  cat("\nPotential optimal partition(s):\n")
  cat(" - Criterion chosen to optimise the number of clusters: ",
      x$args$criterion, "\n")
  if(x$args$criterion %in% c("increasing_step", "decreasing_step")) ##
  {
    cat("   (step quantile chosen: ", x$args$step_quantile,
        " (i.e., only the top", (1 -  x$args$step_quantile) * 100,
        "% ",
        ifelse(x$args$criterion == "increasing_step", "increase", "decrease"),
        " in evaluation metrics",
        " are used as break points for the number of clusters)\n")
  } else if(x$args$criterion == "cutoff")
  {
    cat("   --> cutoff(s) chosen: ", x$args$metric_cutoffs, "\n" )
  }
  cat(" - Optimal partition(s) of clusters for each metric:\n")
  
  cat(paste(paste(names(x$optimal_nb_clusters),
                  sapply(x$optimal_nb_clusters,
                         paste, collapse = " "), sep = " - "),
            collapse = "\n"))
  cat("\n")
}

#' @export
#' @method str bioregion.optimal.n
str.bioregion.optimal.n <- function(object, ...)
{
  args <- list(...)
  if(is.null(args$max.level))
  {
    args$max.level <- 2
  }
  NextMethod("str", object = object, max.level = args$max.level)
}


#' @export
#' @method print bioregion.pairwise.metric
print.bioregion.pairwise.metric <- function(x, ...)
{
  metrics <- colnames(x)[-which(colnames(x) %in%
                                  c("Site1", "Site2", "a", "b",
                                    "c", "A", "B", "C"))]
  cat(paste0("Data.frame of ",
             ifelse(attr(x, "type") == "similarity",
                    "similarity",
                    "dissimilarity"),
             " between sites\n"))
  cat(" - Total number of sites: ", attr(x, "nb_sites"), "\n")
  cat(" - Total number of species: ", attr(x, "nb_species"), "\n")
  cat(" - Number of rows: ", 
      (attr(x, "nb_sites") * (attr(x, "nb_sites") - 1)) / 2, "\n")
  # Warning, next line can be wrong if users alter the object
  cat(" - Number of", ifelse(attr(x, "type") == "similarity",
                             "similarity",
                             "dissimilarity"), "metrics: ",
      length(metrics), "\n")
  cat("\n\n")
  print(as.data.frame(x))
}

#' @export
#' @method `[` bioregion.pairwise.metric
`[.bioregion.pairwise.metric` <- function(x, i, j, ..., drop = TRUE) {
  metric_type <- attributes(x)$type
  nb_sites <- attributes(x)$nb_sites
  nb_species <- attributes(x)$nb_species
  
  class(x) <- "data.frame"
  out <- x[i, j, ..., drop = drop]
  # We keep track of pw metric class & attribute only if the subset is not a vector
  if(class(out) == "data.frame") {
    class(out) <- append("bioregion.pairwise.metric", class(out))
    attributes(out)$type <- metric_type
    attributes(out)$nb_sites <- nb_sites
    attributes(out)$nb_species <- nb_species
  }
  out
}

#' @export
#' @method as.dist bioregion.pairwise.metric
as.dist.bioregion.pairwise.metric <- function(m, diag = FALSE, 
                                              upper = FALSE)
{
  if(ncol(x) > 3) {
    message("More than 3 columns in x: using the third column as the distance",
            "index")
    x <- x[, 1:3]
  }
  matrix.dist <- net_to_mat(x,
                            weight = TRUE, squared = TRUE, symmetrical = TRUE)
  matrix.dist <- stats::as.dist(x, diag = diag, 
                                upper = upper)
  return(matrix.dist)
}
