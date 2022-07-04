#' @export
#' @method str bioRgeo.clusters
str.bioRgeo.clusters <- function(object, ...)
{
  args <- list(...)
  if(is.null(args$max.level))
  {
    args$max.level <- 2
  }
  NextMethod("str", object = object, max.level = args$max.level)
}


#' @export
#' @method print bioRgeo.clusters
print.bioRgeo.clusters <- function(x, ...)
{
  # algorithm name -----
  cat("Clustering results for algorithm : ")
  cat(x$name, "\n")
  if(x$name == "hierarchical_clustering") {
    cat("\t(hierarchical clustering based on a distance matrix)\n")
  }

  # dataset characteristics -----
  cat(" - Number of sites: ", x$inputs$nb_sites, "\n")

  # methodological details -----
  if(x$name == "hierarchical_clustering") {
    cat(" - Name of distance metric: ",
        ifelse(is.null(x$args$index),
                       "Undefined",
                       x$args$index), "\n")
    cat(" - Tree construction method: ", x$args$method, "\n")
    cat(" - Randomization of the distance matrix: ",
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
            x$args$n_clust, "\n")
      }
      if(!is.null(x$args$cut_height))
      {
        cat(" - Heights of cut requested by the user: ",
            round(x$args$cut_height, 3), "\n")
      }
      if(x$args$dynamic_tree_cut)
      {
        cat(paste0(" - Dynamic tree cut method chosen: '", x$args$dynamic_method,
                   "', with minimum cluster size ", x$args$dynamic_minClusterSize, "\n"))
      }

    }

    cat("Clustering results:\n")
    cat(" - Number of partitions: ",
        ncol(x$clusters) - 1, "\n")
    cat(" - Number of clusters: ",
        apply(x$clusters[, 2:ncol(x$clusters), drop = FALSE],
              2, function(y) length(unique(y))),
        "\n")

    if(x$name == "hierarchical_clustering") {
      if(x$args$find_h)
      {
        cat(" - Height of cut of the hierarchical tree:",
            round(x$algorithm$output_cut_height, 3), "\n")
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
#' @method plot bioRgeo.clusters
plot.bioRgeo.clusters <- function(x, ...)
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
      args$ylab <- paste0(x$args$index, " distance")
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
    if(!is.null(x$output_cut_height))
    {
      # abline(h = x$output_cut_height, lty = 3, col = "#756bb1")

      if(length(x$output_cut_height) > 1)
      {
        message("Multiple cut detected, plotting only the first three levels")

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
                         # to do: add border colours from a vector with a distinct colour for each
                         # cluster
                         border = "#377eb8")
    }
  } else
  {
    stop("No plot method for this type of object")
  }
}


#' @export
#' @method print bioRgeo.partition.metrics
print.bioRgeo.partition.metrics <- function(x, ...)
{
  cat("Partition metrics:\n")
  cat(" -", nrow(x$evaluation_df), " partition(s) evaluated\n")
  cat(" - Range of clusters explored: from ", min(x$evaluation_df$n_clusters), " to ",
      max(x$evaluation_df$n_clusters), "\n")
  cat(" - Requested metric(s): ", x$args$eval_metric, "\n")
  cat(" - Metric summary:\n")

  print(data.frame(sapply(x$evaluation_df[x$args$eval_metric],
               function(x) {
                 c(min(x), mean(x), max(x))}),
             row.names = c("Min", "Mean", "Max")))

  if(x$args$partition_optimisation){
    cat("\nPotential optimal partition(s):\n")
    if(length(x$args$eval_metric) > 1)
    {
      cat(" - Metric used for optimisation: ", x$args$eval_metric[1], "\n")
    }
    cat(" - Criterion chosen to optimise the number of clusters: ",
        x$args$criterion, "\n")
    if(x$args$criterion == "step")
    {
      cat("   (step quantile chosen: ", x$args$step_quantile,
          " (i.e., only the top", (1 -  x$args$step_quantile) * 100,
          "% increase in ", x$args$eval_metric[1],
          " are used as break points for the number of clusters)\n")
    } else if(x$args$criterion == "cutoff")
    {
      cat("   --> cutoff(s) chosen: ", x$args$metric_cutoffs, "\n" )
    }
    cat(" - Optimal partition(s) of clusters:\n")
    cat("\n", x$evaluation_df$K[x$evaluation_df$optimal_nclust], "\n\n")
    cat(" - Respective optimal number(s) of clusters:\n")
    cat("\n", x$optimal_nb_clusters, "\n")
  }

  cat("Access the data.frame of metrics with your_object$evaluation_df")
}

#' @export
#' @method str bioRgeo.hierar.tree
str.bioRgeo.hierar.tree <- function(object, ...)
{
  args <- list(...)
  if(is.null(args$max.level))
  {
    args$max.level <- 1
  }
  NextMethod("str", object = object, max.level = args$max.level)
}



#' @export
#' @method print bioRgeo.hierar.tree
print.bioRgeo.hierar.tree <- function(x, ...)
{
  cat("Hierarchical tree based on distances between sites\n")
  cat(" - Number of sites: ", attr(x$dist.matrix, "Size"), "\n")
  cat(" - Tree construction method: ", x$args$method, "\n")
  cat(" - Name of distance column used: ", x$args$index, "\n")
  cat(" - Randomization of the distance matrix: ",
      ifelse(x$args$randomize, paste0("yes, number of trials ",
                                      x$args$n_runs), "no"), "\n")
  cat(" - Cophenetic correlation coefficient: ",
      round(x$final.tree.coph.cor, 3), "\n")


  if(!is.null(x$clusters))
  {
    cat("Clusters obtained from this tree:\n")
    if(!is.null(x$args$n_clust))
    {
      cat(" - Number of clusters requested by the user: ",
          x$args$n_clust, "\n")
      cat(" - Number of clusters in output: ",
          paste(x$output_n_clust, sep = ", "), "\n")

      if(x$args$find_h)
      {
        cat(" - Height of cut of the hierarchical tree:",
            round(x$output_cut_height, 3), "\n")
      } else
      {
        cat(" - Height of cut not searched for.", "\n")
      }
    } else if(!is.null(x$args$cut_height))
    {
      if(length(x$args$cut_height) == 1)
      {
        cat(" - Height of cut requested by the user: ",
            round(x$args$cut_height, 3), "\n")
        cat(" - Number of clusters in output: ", x$output_n_clust, "\n")
      } else
      {
        cat(" - Heights of cut requested by the user: ",
            round(x$args$cut_height, 3), "\n")
        cat(" - Number of clusters in output: \n")
        for(cuts in x$args$cut_height)
        {
          cat("   > ", cuts, ": ", x$output_n_clust[paste0("h_", cuts)], "\n")
        }
      }
    } else if(x$args$dynamic_tree_cut)
    {
      cat(paste0(" - Dynamic tree cut method chosen: '", x$args$dynamic_method,
          "', with minimum cluster size ", x$args$dynamic_minClusterSize, "\n"))
      cat(" - Number of clusters in output: ", x$output_n_clust, "\n")
    }

  }
}


#' @export
#' @method print bioRgeo.pairwise.metric
print.bioRgeo.pairwise.metric <- function(x, ...)
{
  metrics <- colnames(x)[-which(colnames(x) %in%
                                  c("Site1", "Site2", "a", "b",
                                    "c", "A", "B", "C"))]
  cat(paste0("Data.frame of ",
             ifelse(attr(x, "type") == "similarity",
                    "similarity",
                    "distance"),
                    " between sites\n"))
  cat(" - Total number of sites: ", length(unique(c(x$Site1, x$Site2))), "\n")
  cat(" - Number of rows: ", nrow(x), "\n")
  # Warning, next line can be wrong if users alter the object
  cat(" - Number of", ifelse(attr(x, "type") == "similarity",
                              "similarity",
                              "distance"), "metrics: ",
      length(metrics), "\n")
  cat("\n\n")
  print(as.data.frame(x))
}

#' #' @export
#' #' @method print bioRgeo.distance
#' print.bioRgeo.distance <- function(x, ...)
#' {
#'   metrics <- colnames(x)[-which(colnames(x) %in%
#'                                   c("Site1", "Site2", "a", "b",
#'                                     "c", "A", "B", "C"))]
#'   cat("Data.frame of distances between sites\n")
#'   cat(" - Total number of sites: ", length(unique(c(x$Site1, x$Site2))), "\n")
#'   cat(" - Number of rows: ", nrow(x), "\n")
#'   # Warning, next line can be wrong if users alter the object
#'   cat(" - Number of distance metrics: ", length(metrics), "\n")
#'   cat("\n\n")
#'   print(as.data.frame(x))
#' }
