#' Finding communities using the Leiden algorithm
#'
#' This function finds communities in a (un)weighted undirected network based
#' on the Leiden algorithm of Traag, van Eck & Waltman.
#'
#' @param net the output object from [similarity()] or
#' [dissimilarity_to_similarity()]. If a `data.frame` is used, the first two
#' columns represent pairs of sites (or any pair of nodes), and the next
#' column(s) are the similarity indices.
#' 
#' @param weight a `boolean` indicating if the weights should be considered
#' if there are more than two columns.
#' 
#' @param index name or number of the column to use as weight. By default,
#' the third column name of `net` is used.
#'
#' @param objective_function Whether to use the Constant Potts Model (CPM) or
#' modularity. Must be either "CPM" or "modularity".
#'
#' @param resolution_parameter The resolution parameter to use. Higher
#' resolutions lead to more smaller communities, while lower resolutions lead
#' to fewer larger communities.
#'
#' @param beta Parameter affecting the randomness in the Leiden algorithm. This
#' affects only the refinement step of the algorithm.
#'
#' @param initial_membership If provided, the Leiden algorithm will try to
#' improve this provided membership. If no argument is provided, the aglorithm
#' simply starts from the singleton partition.
#'
#' @param n_iterations the number of iterations to iterate the Leiden
#' algorithm. Each iteration may improve the partition further.
#'
#' @param vertex_weights the vertex weights used in the Leiden algorithm. If
#' this is not provided, it will be automatically determined on the basis of
#' the objective_function. Please see the details of this function how to
#' interpret the vertex weights.
#'
#' @param bipartite a `boolean` indicating if the network is bipartite
#' (see Details).
#' 
#' @param site_col name or number for the column of site nodes
#' (i.e. primary nodes).
#' @param species_col name or number for the column of species nodes
#' (i.e. feature nodes).
#' 
#' @param return_node_type a `character` indicating what types of nodes
#' ("sites", "species" or "both") should be returned in the output
#' (`keep_nodes_type="both"` by default).
#' 
#' @param algorithm_in_output a `boolean` indicating if the original output
#' of `communities` should be returned in the output (see Value).
#'
#' @details
#' This function is based on the Leiden algorithm
#' \insertCite{Traag2019}{bioRgeo} as implemented in the
#' [igraph](https://cran.r-project.org/web/packages/igraph/index.html)
#' package ([cluster_leiden][igraph::cluster_leiden]).
#'
#' @note
#' Although this algorithm was not primarily designed to deal with bipartite
#' network, it is possible to consider the bipartite network as unipartite
#' network (`bipartite = TRUE`).
#'
#' Do not forget to indicate which of the first two columns is
#' dedicated to the site nodes (i.e. primary nodes) and species nodes (i.e.
#' feature nodes) using the arguments `site_col` and `species_col`.
#' The type of nodes returned in the output can be chosen with the argument
#' `return_node_type` equal to `"both"` to keep both types of nodes,
#' `"sites"` to preserve only the sites nodes and `"species"` to
#' preserve only the species nodes.
#'
#' @return
#' A `list` of class `bioRgeo.clusters` with five slots:
#' \enumerate{
#' \item{**name**: `character string` containing the name of the algorithm}
#' \item{**args**: `list` of input arguments as provided by the user}
#' \item{**inputs**: `list` of characteristics of the input dataset}
#' \item{**algorithm**: `list` of all objects associated with the
#'  clustering procedure, such as original cluster objects (only if
#'  `algorithm_in_output = TRUE`)}
#' \item{**clusters**: `data.frame` containing the clustering results}}
#'
#' In the `algorithm` slot, if `algorithm_in_output = TRUE`, users can
#' find an "communities" object, output of
#' [cluster_leiden][igraph::cluster_leiden].
#'
#' @author
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}),
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) and
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' 
#' @examples
#' \dontrun{
#' comat <- matrix(sample(1000, 50), 5, 10)
#' rownames(comat) <- paste0("Site", 1:5)
#' colnames(comat) <- paste0("Species", 1:10)
#'
#' net <- similarity(comat, metric = "Simpson")
#' com <- netclu_leiden(net)
#' 
#' net_bip <- mat_to_net(comat, weight = TRUE)
#' clust2 <- netclu_leiden(net_bip, bipartite = TRUE)
#' }
#' 
#' @references
#' \insertRef{Traag2019}{bioRgeo}
#' 
#' @importFrom igraph graph_from_data_frame cluster_leiden
#' 
#' @export

netclu_leiden <- function(net,
                          weight = TRUE,
                          index = names(net)[3],
                          objective_function = c("CPM", "modularity"),
                          resolution_parameter = 1,
                          beta = 0.01,
                          initial_membership = NULL,
                          n_iterations = 2,
                          vertex_weights = NULL,
                          bipartite = FALSE,
                          site_col = 1,
                          species_col = 2,
                          return_node_type = "both",
                          algorithm_in_output = TRUE) {
  
  # Control input net
  controls(args = NULL, data = net, type = "input_bioRgeo.pairwise.metric")
  controls(args = NULL, data = net, type = "input_net")
  
  # Control input weight & index
  controls(args = weight, data = net, type = "input_net_weight")
  if (weight) {
    controls(args = index, data = net, type = "input_net_index")
    net[, 3] <- net[, index]
    net <- net[, 1:3]
    controls(args = NULL, data = net, type = "input_net_index_value")
  }
  
  # Control input bipartite
  controls(args = bipartite, data = NULL, type = "boolean")
  isbip <- bipartite
  if (isbip) {
    controls(args = NULL, data = net, type = "input_net_bip")
    controls(args = site_col, data = net, type = "input_net_bip_col")
    controls(args = species_col, data = net, type = "input_net_bip_col")
    if (!(return_node_type %in% c("both", "sites", "species"))) {
      stop("Please choose return_node_type among the followings values:
both, sites and species", call. = FALSE)
    }
  }
  
  # Control input directed
  controls(args = NULL, data = net, type = "input_net_isdirected")
  
  # Control algorithm_in_output
  controls(args = algorithm_in_output, data = NULL, type = "boolean")
  
  # Controls for other arguments
  controls(args = objective_function, data = NULL, type = "character")
  if(!all(objective_function %in% c("CPM", "modularity"))){
    stop("objective_function must be either 'CPM' or 'modularity'.")
  }
  
  controls(args = resolution_parameter, data = NULL, type = "positive_integer")
  controls(args = beta, data = NULL, type = "numeric")
  # controls(args = initial_membership, data = NULL, type = "boolean")
  controls(args = n_iterations, data = NULL, type = "positive_integer")
  # controls(args = vertex_weights, data = NULL, type = "boolean")
  
  # Prepare input
  if (isbip) {
    idprim <- as.character(net[, site_col])
    idprim <- idprim[!duplicated(idprim)]
    nbsites <- length(idprim)
    idfeat <- as.character(net[, species_col])
    idfeat <- idfeat[!duplicated(idfeat)]
    
    idnode <- c(idprim, idfeat)
    idnode <- data.frame(ID = 1:length(idnode), ID_NODE = idnode)
    netemp <- data.frame(
      node1 = idnode[match(net[, site_col], idnode[, 2]), 1],
      node2 = idnode[match(net[, species_col], idnode[, 2]), 1]
    )
  } else {
    idnode1 <- as.character(net[, 1])
    idnode2 <- as.character(net[, 2])
    if (isbip) {
      message("The network seems to be bipartite! 
The bipartite argument should probably be set to TRUE.")
    }
    idnode <- c(idnode1, idnode2)
    idnode <- idnode[!duplicated(idnode)]
    nbsites <- length(idnode)
    idnode <- data.frame(ID = 1:length(idnode), ID_NODE = idnode)
    netemp <- data.frame(
      node1 = idnode[match(net[, 1], idnode[, 2]), 1],
      node2 = idnode[match(net[, 2], idnode[, 2]), 1]
    )
  }
  
  if (weight) {
    netemp <- cbind(netemp, net[, 3])
    netemp <- netemp[netemp[, 3] > 0, ]
    colnames(netemp)[3] <- "weight"
  }
  
  # Class preparation
  outputs <- list(name = "netclu_leiden")
  
  outputs$args <- list(
    weight = weight,
    index = index,
    bipartite = bipartite,
    site_col = site_col,
    species_col = species_col,
    return_node_type = return_node_type,
    algorithm_in_output = algorithm_in_output
  )
  
  outputs$inputs <- list(
    bipartite = isbip,
    weight = weight,
    pairwise = ifelse(isbip, FALSE, TRUE),
    pairwise_metric = ifelse(isbip, NA, index),
    dissimilarity = FALSE,
    nb_sites = nbsites
  )
  
  outputs$algorithm <- list()
  
  # Run algo
  net <- igraph::graph_from_data_frame(netemp, directed = FALSE)
  outalg <- igraph::cluster_leiden(
    graph = net,
    weight = weight,
    objective_function = objective_function,
    resolution_parameter = resolution_parameter,
    beta = beta,
    initial_membership = initial_membership,
    n_iterations = n_iterations,
    vertex_weights = vertex_weights)
  comtemp <- cbind(as.numeric(outalg$names), as.numeric(outalg$membership))
  
  com <- data.frame(ID = idnode[, 2], Com = 0)
  com[match(comtemp[, 1], idnode[, 1]), 2] <- comtemp[, 2]
  
  # Rename and reorder columns
  com <- knbclu(com)
  
  # Add attributes and return_node_type
  if (isbip) {
    attr(com, "node_type") <- rep("site", dim(com)[1])
    attributes(com)$node_type[!is.na(match(com[, 1], idfeat))] <- "species"
    if (return_node_type == "sites") {
      com <- com[attributes(com)$node_type == "site", ]
    }
    if (return_node_type == "species") {
      com <- com[attributes(com)$node_type == "species", ]
    }
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
                                             drop = FALSE
    ],
    n_clust = apply(
      outputs$clusters[, 2:length(outputs$clusters), drop = FALSE],
      2, function(x) length(unique(x))))
  
  # Return outputs
  class(outputs) <- append("bioRgeo.clusters", class(outputs))
  return(outputs)
}
