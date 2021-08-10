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
