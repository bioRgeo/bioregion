#' @export
#' @method print bioRgeo.nclust.tree
print.bioRgeo.nclust.tree <- function(x, ...)
{
  cat("Optimal number(s) of clusters search:\n")
  cat(" - Range of clusters explored: from ", x$args$k_min, " to ", x$args$k_max, "\n")
  cat(" - Requested metric(s): ", x$args$eval_metric, "\n")
  if(length(x$args$eval_metric) > 1)
  {
    cat(" - Metric used for optimisation: ", x$args$eval_metric[1], "\n")
  }
  cat(" - Criterion chosen to optimise the number of clusters: ", x$args$criterion, "\n")
  if(x$args$criterion == "gap")
  {
    cat("   --> step quantile chosen: ", x$args$gap_quantile,
        " (i.e., only the top", (1 -  x$args$gap_quantile) * 100,
        "% increase in ", x$args$eval_metric[1],
        " are used as break points for the number of clusters)\n")
  } else if(x$args$criterion == "cutoff")
  {
    cat("   --> cutoff(s) chosen: ", x$args$metric_cutoffs, "\n" )
  }
  cat("\nOptimal number(s) of clusters: \n")
  cat(x$optimal_nb_clusters)
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
#' @method plot bioRgeo.hierar.tree
plot.bioRgeo.hierar.tree <- function(x, ...)
{
  plot(x$final.tree, ...)
}

#' @export
#' @method print bioRgeo.hierar.tree
print.bioRgeo.hierar.tree <- function(x, ...)
{
  cat("Hierarchical tree based on distances between sites\n")
  cat(" - Number of sites: ", length(x$dist.matrix), "\n")
  cat(" - Tree construction method: ", x$args$method, "\n")
  cat(" - Name of distance column used: ", x$args$index, "\n")
  cat(" - Randomization of the distance matrix: ", ifelse(x$args$randomize, paste0("yes, number of trials ", x$args$n_runs), "no"), "\n")


  if(!is.null(x$clusters))
  {
    cat("Clusters obtained from this tree:\n")
    if(!is.null(x$args$n_clust))
    {
      cat(" - Number of clusters requested by the user: ", x$args$n_clust, "\n")
      cat(" - Number of clusters in output: ", paste(x$output_n_clust, sep = ", "), "\n")

      if(x$args$find_h)
      {
        cat(" - Height of cut of the hierarchical tree:", x$output_cut_height)
      } else
      {
        cat(" - Height of cut not searched for.")
      }
    }
    if(!is.null(x$args$cut_height))
    {
      if(length(x$args$cut_height) == 1)
      {
        cat(" - Height of cut requested by the user: ", x$args$cut_height, "\n")
        cat(" - Number of clusters in output: ", x$output_n_clust, "\n")
      } else
      {
        cat(" - Heights of cut requested by the user: ", x$args$cut_height, "\n")
        cat(" - Number of clusters in output: \n")
        for(cuts in x$args$cut_height)
        {
          cat("   > ", cuts, ": ", x$output_n_clust[paste0("h_", cuts)], "\n")
        }
      }
    }
  }
}


#' @export
#' @method print bioRgeo.similarity
print.bioRgeo.similarity <- function(x, ...)
{
  metrics <- colnames(x)[-which(colnames(x) %in% c("Site1", "Site2", "a", "b", "c", "A", "B", "C"))]
  cat("Data.frame of similarity between sites\n")
  cat(" - Total number of sites: ", length(unique(c(x$Site1, x$Site2))), "\n")
  cat(" - Number of rows: ", nrow(x), "\n")
  # Warning, next line can be wrong if users alter the object
  cat(" - Number of similarity metrics: ", length(metrics), "\n")
  cat("\n\n")
  print(as.data.frame(x))
}

#' @export
#' @method print bioRgeo.distance
print.bioRgeo.distance <- function(x, ...)
{
  metrics <- colnames(x)[-which(colnames(x) %in% c("Site1", "Site2", "a", "b", "c", "A", "B", "C"))]
  cat("Data.frame of distances between sites\n")
  cat(" - Total number of sites: ", length(unique(c(x$Site1, x$Site2))), "\n")
  cat(" - Number of rows: ", nrow(x), "\n")
  # Warning, next line can be wrong if users alter the object
  cat(" - Number of distance metrics: ", length(metrics), "\n")
  cat("\n\n")
  print(as.data.frame(x))
}
