#' Community structure detection in weighted bipartite network via modularity
#' optimization
#'
#' This function takes a bipartite weighted graph and computes modules by
#' applying Newman’s modularity measure in a bipartite weighted version to it.
#'
#' @param net a `data.frame` representing a bipartite network with the two
#' first columns as undirected links between pair of nodes and and the next
#' column(s) are the weight of the links.
#' 
#' @param weight a `boolean` indicating if the weights should be considered
#' if there are more than two columns (see Note).
#' 
#' @param index name or number of the column to use as weight. By default,
#' the third column name of `net` is used.
#' 
#' @param site_col name or number for the column of site nodes (i.e. primary
#' nodes).
#' 
#' @param species_col name or number for the column of species nodes (i.e.
#' feature nodes).
#' 
#' @param return_node_type a `character` indicating what types of nodes
#' ("sites", "species" or "both") should be returned in the output
#' (`return_node_type = "both"` by default).
#' 
#' @param forceLPA a `boolean` indicating if the even faster pure
#' LPA-algorithm of Beckett should be used? DIRT-LPA, the default, is less
#' likely to get trapped in a local minimum, but is slightly slower. Defaults
#' to FALSE.
#' 
#' @param algorithm_in_output a `boolean` indicating if the original output
#' of `computeModules` should be returned in the output (see Value).
#' Default to TRUE.
#' 
#' @details
#' This function is based on the modularity optimization algorithm provided by
#' Stephen Beckett \insertCite{Beckett2016}{bioregion} as implemented in the
#' [bipartite](https://cran.r-project.org/package=bipartite)
#' package ([computeModules][bipartite::computeModules]).
#'
#' @note
#' Beckett has been designed to deal with weighted bipartite networks. Note
#' that if `weight = FALSE`, a weight of 1 will be assigned to each pair of
#' nodes. Do not forget to indicate which of the first two columns is
#' dedicated to the site nodes (i.e. primary nodes) and species nodes (i.e.
#' feature nodes) using the arguments `site_col` and `species_col`. The type of
#' nodes returned in the output can be chosen with the argument
#' `return_node_type` equal to `"both"` to keep both types of nodes,`"sites"`
#' to preserve only the sites nodes and `"species"` to preserve only the
#' species nodes.
#'
#' @return
#' A `list` of class `bioregion.clusters` with five slots:
#' \enumerate{
#' \item{**name**: `character string` containing the name of the algorithm}
#' \item{**args**: `list` of input arguments as provided by the user}
#' \item{**inputs**: `list` of characteristics of the clustering process}
#' \item{**algorithm**: `list` of all objects associated with the
#'  clustering procedure, such as original cluster objects (only if
#'  `algorithm_in_output = TRUE`)}
#' \item{**clusters**: `data.frame` containing the clustering results}}
#'
#' In the `algorithm` slot, if `algorithm_in_output = TRUE`, users can find an
#' object of class "moduleWeb", output of
#' [computeModules][bipartite::computeModules].
#'
#' @author
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}),
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) and
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' 
#' @seealso [netclu_infomap], [netclu_oslom]
#' 
#' @examples
#' net <- data.frame(
#'   Site = c(rep("A", 2), rep("B", 3), rep("C", 2)),
#'   Species = c("a", "b", "a", "c", "d", "b", "d"),
#'   Weight = c(10, 100, 1, 20, 50, 10, 20))
#'
#' com <- netclu_beckett(net)
#' 
#' @references
#' \insertRef{Beckett2016}{bioregion}
#' 
#' @importFrom bipartite computeModules
#' 
#' @export
netclu_beckett <- function(net,
                           weight = TRUE,
                           index = names(net)[3],
                           forceLPA = FALSE,
                           site_col = 1,
                           species_col = 2,
                           return_node_type = "both",
                           algorithm_in_output = TRUE){
  
  # Controls inputs
  controls(args = forceLPA, data = NULL, type = "boolean")
  controls(args = algorithm_in_output, data = NULL, type = "boolean")
  
  controls(args = NULL, data = net, type = "input_net")
  
  controls(args = NULL, data = net, type = "input_net_bip")
  if(site_col == species_col){
    stop("site_col and species_col should not be the same.", call. = FALSE)
  }
  controls(args = site_col, data = net, type = "input_net_bip_col")
  controls(args = species_col, data = net, type = "input_net_bip_col")
  controls(args = return_node_type, data = NULL, type = "character")
  if(!(return_node_type %in% c("both", "sites", "species"))) {
    stop("Please choose return_node_type among the followings values:
both, sites or species", call. = FALSE)}
  
  controls(args = weight, data = net, type = "input_net_weight")
  if (weight) {
    controls(args = index, data = net, type = "input_net_index")
    net[, 3] <- net[, index]
    net <- net[, 1:3]
    controls(args = NULL, data = net, type = "input_net_index_positive_value")
  }
  
  # Prepare input
  idprim <- as.character(net[, site_col])
  idprim <- idprim[!duplicated(idprim)]
  nbsites <- length(idprim)
  idfeat <- as.character(net[, species_col])
  idfeat <- idfeat[!duplicated(idfeat)]
  
  idnode <- c(idprim, idfeat)
  idnode <- data.frame(ID = 1:length(idnode), ID_NODE = idnode)
  
  netemp <- data.frame(
    node1 = idnode[match(net[, site_col], idnode[, 2]), 1],
    node2 = idnode[match(net[, species_col], idnode[, 2]), 1])
  
  if(weight) {
    netemp <- cbind(netemp, net[, 3])
    netemp <- netemp[netemp[, 3] > 0, ]
    colnames(netemp)[3] <- "weight"
  } else {
    netemp$weight <- 1
  }
  
  # Class preparation
  outputs <- list(name = "netclu_beckett")
  
  outputs$args <- list(weight = weight,
                       index = index,
                       site_col = site_col,
                       species_col = species_col,
                       return_node_type = return_node_type,
                       forceLPA = forceLPA,
                       algorithm_in_output = algorithm_in_output)
  
  outputs$inputs <- list(
    bipartite = TRUE,
    weight = weight,
    pairwise = FALSE,
    pairwise_metric = NA,
    dissimilarity = FALSE,
    nb_sites = nbsites,
    hierarchical = FALSE)
  
  outputs$algorithm <- list()
  
  # Transform netemp into a contingency table
  comat <- net_to_mat(netemp, weight = weight)
  
  # Run algo
  outalg <- bipartite::computeModules(comat, forceLPA = forceLPA)
  comtemp <- outalg@modules[-1, -c(1, 2)]
  comtemp <- cbind(c(as.numeric(rownames(comat)),
                     as.numeric(colnames(comat))),
                   apply(comtemp, 2, function(x) which(x > 0)))
  
  com <- data.frame(ID = idnode[, 2], Com = 0)
  com[match(comtemp[, 1], idnode[, 1]), 2] <- comtemp[, 2]
  
  # Rename and reorder columns
  com <- knbclu(com)
  
  # Add attributes and return_node_type
  attr(com, "node_type") <- rep("site", dim(com)[1])
  attributes(com)$node_type[!is.na(match(com[, 1], idfeat))] <- "species"
  if (return_node_type == "sites") {
    com <- com[attributes(com)$node_type == "site", ]
  }
  if (return_node_type == "species") {
    com <- com[attributes(com)$node_type == "species", ]
  }
  
  # Set algorithm in outputs
  if (!algorithm_in_output) {
    outalg <- NA
  }
  outputs$algorithm <- outalg
  
  # Set clusters and cluster_info in output
  outputs$clusters <- com
  outputs$cluster_info <- data.frame(
    partition_name = names(outputs$clusters)[2:length(outputs$clusters),
                                             drop = FALSE],
    n_clust = apply(
      outputs$clusters[, 2:length(outputs$clusters), drop = FALSE],
      2, function(x) length(unique(x))))
  
  # Return outputs
  class(outputs) <- append("bioregion.clusters", class(outputs))
  return(outputs)
}
