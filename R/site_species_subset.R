#' Extract a subset of sites or species from a bioregion.clusters object
#'
#' This function extracts a subset of nodes according to its type (sites or 
#' species) from a bioregion.clusters object containing both types of 
#' nodes (sites and species).
#'
#' @param clusters an object of class `bioregion.clusters`.
#'
#' @param node_type a `character` indicating what types of nodes
#' ("site" or "species") should be extracted
#' (`node_type = "site"` by default).
#'
#' @return An object of class `bioregion.clusters` with a given node type (sites 
#' or species).
#' 
#' @note 
#' The network clustering functions (prefix `netclu_`) may return both types of 
#' nodes (sites and species) when applied on bipartite networks 
#' (argument `bipartite`). In this case, the type of nodes returned in the 
#' output can be chosen with the argument `return_node_type`. This function 
#' allows to retrieve a particular type of nodes (sites or species) from the 
#' output and modify the return_node_type accordingly.
#'
#' @author
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}),
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) and
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
site_species_subset <- function(clusters, node_type = "site") {

  # Control node_type
  controls(args = node_type, data = NULL, type = "character")
  if (!(node_type %in% c("site", "species"))) {
    stop("Please choose node_type among the followings values:
sites and species", call. = FALSE)
  }
  
  # Control input 
  if (!inherits(clusters, "bioregion.clusters")) {
    stop("clusters must be a bioregion.clusters object.",
         call. = FALSE
    )
  }
  
  func <- clusters$name
  if(substr(func, 1,7) != "netclu_"){
    stop("clusters must be an output of a 'netclu_' function.",
         call. = FALSE
    )
  }
  
  bip <- FALSE
  if(func == "netclu_beckett"){
    bip <- TRUE
  } else if(func == "netclu_infomap"){
    if(clusters$args$bipartite | clusters$args$bipartite_version){
      bip <- TRUE
    }
  }else{
    if(clusters$args$bipartite){
      bip <- TRUE
    }
  }
  if(!bip){
    stop("clusters must be based on a bipartite network.",
         call. = FALSE
    )
  }

  if(clusters$args$return_node_type != "both"){
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
  clusters$args$return_node_type <- node_type

  # Return output
  return(clusters)
}
