#' Convert betapart dissimilarity to bioregion dissimilarity (DEPRECATED)
#'
#' This function converts dissimilarity results produced by the betapart package
#' (and packages using betapart, such as phyloregion) into a dissimilarity 
#' object compatible with the bioregion package. This function only converts 
#' object types to make them compatible with bioregion; it does not modify the 
#' beta-diversity values. This function allows the inclusion of phylogenetic 
#' beta diversity to compute bioregions with bioregion.
#'
#' @param betapart_result An object produced by the betapart package (e.g., 
#' using the `beta.pair` function).
#' 
#' @return 
#' A dissimilarity object of class `bioregion.pairwise`, 
#' compatible with the bioregion package.
#' 
#' @seealso 
#' This function is deprecated, use [as_bioregion_pairwise] instead.
#' 
#' @author
#' Boris Leroy (\email{leroy.boris@gmail.com}) \cr
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) \cr
#' Pierre Denelle (\email{pierre.denelle@gmail.com})
#' 
#' @examples
#' comat <- matrix(sample(0:1000, size = 50, replace = TRUE,
#' prob = 1 / 1:1001), 5, 10)
#' rownames(comat) <- paste0("Site", 1:5)
#' colnames(comat) <- paste0("Species", 1:10)
#' 
#' \dontrun{
#' beta_div <- betapart::beta.pair.abund(comat)
#' betapart_to_bioregion(beta_div)
#' }
#' 
#' @export
betapart_to_bioregion <- function(betapart_result) {
  
  .Deprecated("as_bioregion_pairwise")
  
  if (!inherits(betapart_result, "list")) {
    stop("betapart_result must be a valid object from the betapart package.",
         call. = FALSE)
  }
 
  nb_sites <- nrow(betapart_result[[1]])
  
  site_labels <- attr(betapart_result[[1]], "Labels")
  metric_names <- names(betapart_result)
  
  result_df <-  mat_to_net(as.matrix(betapart_result[[1]]), 
                           weight = TRUE, 
                           remove_zeroes = FALSE, 
                           include_lower = FALSE, 
                           include_diag = FALSE)
  
  # Loop over other metrics and extract pairwise values
  for (metric_name in metric_names[2:length(metric_names)]) {
    metric_values <-  mat_to_net(as.matrix(betapart_result[[metric_name]]), 
                                 weight = TRUE, 
                                 remove_zeroes = FALSE, 
                                 include_lower = FALSE, 
                                 include_diag = FALSE)
    result_df[[metric_name]] <- metric_values$Weight
  }
  
  attr(result_df, "type") <- "dissimilarity"
  attr(result_df, "nb_sites") <- nb_sites
  attr(result_df, "nb_species") <- NA
  
  class(result_df) <- append("bioregion.pairwise", class(result_df))
  
  return(result_df)
}
