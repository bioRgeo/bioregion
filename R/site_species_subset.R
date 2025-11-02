#' Extract a subset of sites or species from a `bioregion.clusters` object
#'
#' This function extracts a subset of nodes based on their type (`"site"` or 
#' `"species"`) from a `bioregion.clusters` object, which contains both types of 
#' nodes (sites and species).
#'
#' @param clusters An object of class `bioregion.clusters`.
#' 
#' @param node_type A `character` string indicating the type of nodes to 
#' extract. Possible values are `"site"` or `"species"`. The default is 
#' `"site"`.
#'
#' @return 
#' An object of class `bioregion.clusters` containing only the specified 
#' node type (sites or species).
#' 
#' @note 
#' Some `bioregion.clusters` objects may contain both types of nodes (sites and 
#' species). This information is available in the `$inputs$node_type` slot. 
#' 
#' This function allows you to extract a specific type of node 
#' (either sites or species) from any `bioregion.clusters` object that 
#' includes both.
#'
#' @author
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) \cr
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) \cr
#' Boris Leroy (\email{leroy.boris@gmail.com})
#'
#' @examples
#' net <- data.frame(
#'   Site = c(rep("A", 2), rep("B", 3), rep("C", 2)),
#'   Species = c("a", "b", "a", "c", "d", "b", "d"),
#'   Weight = c(10, 100, 1, 20, 50, 10, 20)
#' )
#'
#' clusters <- netclu_louvain(net, lang = "igraph", bipartite = TRUE)
#' 
#' clusters_sites <- site_species_subset(clusters, node_type = "site")
#'
#' @export
site_species_subset <- function(clusters, 
                                node_type = "site") {

  # Control clusters
  controls(args = NULL, 
           data = clusters, 
           type ="input_bioregionalization")
  
  # Control node_type
  controls(args = node_type, data = NULL, type = "character")
  if (!(node_type %in% c("site", "species"))) {
    stop(paste0("Please choose node_type from the following:\n",
                "sites and species"), 
                call. = FALSE)
  }
  
  # Control node_type in clusters
  if(clusters$inputs$node_type != "both"){
    stop("clusters must contain both types of node.",
         call. = FALSE
    )
  }
  
  # Get type
  if(node_type == "site"){
    clusters$clusters <- clusters$clusters[
      attributes(clusters$clusters)$node_type == "site", ]
  }
  if(node_type == "species"){
    clusters$clusters <- clusters$clusters[
      attributes(clusters$clusters)$node_type == "species", ]
  }
  
  # Update return_node_type
  clusters$inputs$node_type <- node_type
  #if(!is.null(clusters$args$return_node_type)){
  #  clusters$args$return_node_type <- node_type
  #}
  
  # Return output
  return(clusters)
}
