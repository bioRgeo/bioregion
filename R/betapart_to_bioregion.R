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
#'   beta_div <- betapart::beta.pair.abund(comat)
#'   betapart_to_bioregion(beta_div)
#' }
#' 
#' @export
betapart_to_bioregion <- function(betapart_result) {
  
  stop("This function is deprecated, please use as_bioregion_pairwise instead.",
       call. = FALSE)

}
