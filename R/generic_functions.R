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
      cat(" - Number of clusters in output: ", max(x$clusters$cluster), "\n")
      if(x$args$find_h)
      {
        cat(" - Height of cut of the hierarchical tree:", attr(x$clusters, "cut_height"))
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
        cat(" - Number of clusters in output: ", max(x$clusters$cluster), "\n")
      } else
      {
        cat(" - Heights of cut requested by the user: ", x$args$cut_height, "\n")
        cat(" - Number of clusters in output: \n")
        for(cuts in x$args$cut_height)
        {
          cat("   > ", cuts, ": ", max(x$clusters[, paste0("cluster.", cuts)]), "\n")
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
