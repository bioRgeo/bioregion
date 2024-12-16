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
#' @param cut_weight a minimal weight value. If `weight` is TRUE, the links 
#' between sites with a weight strictly lower than this value will not be 
#' considered (O by default).
#' 
#' @param index name or number of the column to use as weight. By default,
#' the third column name of `net` is used.
#' 
#' @param seed for the random number generator (NULL for random by default).
#'
#' @param objective_function a string indicating the objective function to use,
#' the Constant Potts Model ("CPM") or "modularity" ("CPM" by default).
#'
#' @param resolution_parameter the resolution parameter to use. Higher
#' resolutions lead to more smaller communities, while lower resolutions lead
#' to fewer larger communities.
#'
#' @param beta parameter affecting the randomness in the Leiden algorithm. This
#' affects only the refinement step of the algorithm.
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
#' 
#' @param species_col name or number for the column of species nodes
#' (i.e. feature nodes).
#' 
#' @param return_node_type a `character` indicating what types of nodes
#' ("site", "species" or "both") should be returned in the output
#' (`return_node_type = "both"` by default).
#' 
#' @param algorithm_in_output a `boolean` indicating if the original output
#' of [cluster_leiden][igraph::cluster_leiden] should be returned in the 
#' output (`TRUE` by default, see Value).
#'
#' @details
#' This function is based on the Leiden algorithm
#' (Traag et al., 2019) as implemented in the
#' [igraph](https://cran.r-project.org/package=igraph)
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
#' `"site"` to preserve only the sites nodes and `"species"` to
#' preserve only the species nodes.
#'
#' @return
#' A `list` of class `bioregion.clusters` with five slots:
#' \enumerate{
#' \item{**name**: `character` containing the name of the algorithm}
#' \item{**args**: `list` of input arguments as provided by the user}
#' \item{**inputs**: `list` of characteristics of the clustering process}
#' \item{**algorithm**: `list` of all objects associated with the
#'  clustering procedure, such as original cluster objects (only if
#'  `algorithm_in_output = TRUE`)}
#' \item{**clusters**: `data.frame` containing the clustering results}}
#'
#' In the `algorithm` slot, if `algorithm_in_output = TRUE`, users can
#' find the output of
#' [cluster_leiden][igraph::cluster_leiden].
#'
#' @author
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) \cr
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) \cr
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' 
#' @examples
#' comat <- matrix(sample(1000, 50), 5, 10)
#' rownames(comat) <- paste0("Site", 1:5)
#' colnames(comat) <- paste0("Species", 1:10)
#'
#' net <- similarity(comat, metric = "Simpson")
#' com <- netclu_leiden(net)
#' 
#' net_bip <- mat_to_net(comat, weight = TRUE)
#' clust2 <- netclu_leiden(net_bip, bipartite = TRUE)
#' 
#' @references 
#' Traag VA, Waltman L & Van Eck NJ (2019) From Louvain to Leiden: guaranteeing 
#' well-connected communities. \emph{Scientific reports}, 9(1), 5233. 
#' 
#' @importFrom igraph graph_from_data_frame cluster_leiden
#' 
#' @export

netclu_leiden <- function(net,
                          weight = TRUE,
                          cut_weight = 0,
                          index = names(net)[3],
                          seed = NULL,
                          objective_function = "CPM",
                          resolution_parameter = 1,
                          beta = 0.01,
                          n_iterations = 2,
                          vertex_weights = NULL,
                          bipartite = FALSE,
                          site_col = 1,
                          species_col = 2,
                          return_node_type = "both",
                          algorithm_in_output = TRUE) {
  
  # Control input net (+ check similarity if not bipartite)
  controls(args = bipartite, data = NULL, type = "boolean")
  isbip <- bipartite
  if(!isbip){
    controls(args = NULL, data = net, type = "input_similarity")
  }
  controls(args = NULL, data = net, type = "input_net")
  
  # Convert tibble into dataframe
  if(inherits(net, "tbl_df")){
    net <- as.data.frame(net)
  }
  
  # Control input weight & index
  controls(args = weight, data = net, type = "input_net_weight")
  if (weight) {
    controls(args = cut_weight, data = net, type = "positive_numeric")
    controls(args = index, data = net, type = "input_net_index")
    net[, 3] <- net[, index]
    net <- net[, 1:3]
    controls(args = NULL, data = net, type = "input_net_index_positive_value")
  }
  
  # Control input bipartite
  if (isbip) {
    controls(args = NULL, data = net, type = "input_net_bip")
    if(site_col == species_col){
      stop("site_col and species_col should not be the same.", call. = FALSE)
    }
    controls(args = site_col, data = net, type = "input_net_bip_col")
    controls(args = species_col, data = net, type = "input_net_bip_col")
    controls(args = return_node_type, data = NULL, type = "character")
    if (!(return_node_type %in% c("both", "site", "species"))) {
      stop("Please choose return_node_type among the followings values:
both, sites or species", call. = FALSE)
    }
  }
  
  # Control input loop or directed
  controls(args = NULL, data = net, type = "input_net_isloop")
  controls(args = NULL, data = net, type = "input_net_isdirected")
  
  # Control algorithm_in_output
  controls(args = algorithm_in_output, data = NULL, type = "boolean")
  
  # Control parameters LEIDEN
  if(!is.null(seed)){
    controls(args = seed, data = NULL, type = "strict_positive_integer")
  }
  controls(args = objective_function, data = NULL, type = "character")
  if (!(objective_function %in% c("CPM", "modularity"))) {
    stop("Please choose objective_function among the following values: 
CPM or modularity", call. = FALSE)
  }
  controls(args = resolution_parameter, data = NULL, type = "positive_integer")
  controls(args = beta, data = NULL, type = "numeric")
  controls(args = n_iterations, data = NULL, type = "strict_positive_integer")
  if(!is.null(vertex_weights)){
    controls(args = vertex_weights, data = NULL, type = "numeric_vector")
  }
  
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
    netemp <- netemp[netemp[, 3] > cut_weight, ]
    colnames(netemp)[3] <- "weight"
  }
  
  # Class preparation
  outputs <- list(name = "netclu_leiden")
  
  outputs$args <- list(
    weight = weight,
    cut_weight = cut_weight,
    index = index,
    seed = seed,
    objective_function = objective_function,
    resolution_parameter = resolution_parameter,
    beta = beta,
    n_iterations = n_iterations,
    vertex_weights = vertex_weights,
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
    pairwise_metric = ifelse(!isbip & weight, 
                             ifelse(is.numeric(index), names(net)[3], index), 
                             NA),
    dissimilarity = FALSE,
    nb_sites = nbsites,
    hierarchical = FALSE
  )
  
  outputs$algorithm <- list()
  
  # Run algo (with seed)
  net <- igraph::graph_from_data_frame(netemp, directed = FALSE)
  if(is.null(seed)){
    outalg <- igraph::cluster_leiden(
      graph = net,
      objective_function = objective_function,
      resolution = resolution_parameter,
      beta = beta,
      n_iterations = n_iterations,
      vertex_weights = vertex_weights)
  }else{
    set.seed(seed)
    outalg <- igraph::cluster_leiden(
      graph = net,
      objective_function = objective_function,
      resolution = resolution_parameter,
      beta = beta,
      n_iterations = n_iterations,
      vertex_weights = vertex_weights)
    rm(.Random.seed, envir=globalenv())
  }

  comtemp <- cbind(as.numeric(outalg$names), as.numeric(outalg$membership))
  
  com <- data.frame(ID = idnode[, 2], Com = NA)
  com[match(comtemp[, 1], idnode[, 1]), 2] <- comtemp[, 2]
  
  # Rename and reorder columns
  com <- knbclu(com)
  
  # Add attributes and return_node_type
  if (isbip) {
    attr(com, "node_type") <- rep("site", dim(com)[1])
    attributes(com)$node_type[!is.na(match(com[, 1], idfeat))] <- "species"
    if (return_node_type == "site") {
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
      2, function(x) length(unique(x[!is.na(x)]))))
  
  # Return outputs
  class(outputs) <- append("bioregion.clusters", class(outputs))
  return(outputs)
}
