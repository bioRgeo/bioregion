#' Community structure detection via short random walks
#'
#' This function finds communities in a (un)weighted undirected network via
#' short random walks.
#'
#' @param net The output object from [similarity()] or
#' [dissimilarity_to_similarity()]. If a `data.frame` is used, the first two
#' columns represent pairs of sites (or any pair of nodes), and the next
#' column(s) are the similarity indices.
#' 
#' @param weight A `boolean` indicating if the weights should be considered
#' if there are more than two columns.
#' 
#' @param cut_weight A minimal weight value. If `weight` is TRUE, the links 
#' between sites with a weight strictly lower than this value will not be 
#' considered (0 by default).
#' 
#' @param index Name or number of the column to use as weight. By default,
#' the third column name of `net` is used.
#' 
#' @param steps The length of the random walks to perform.
#' 
#' @param bipartite A `boolean` indicating if the network is bipartite
#' (see Details).
#' 
#' @param site_col Name or number for the column of site nodes
#' (i.e. primary nodes).
#' 
#' @param species_col Name or number for the column of species nodes
#' (i.e. feature nodes).
#' 
#' @param return_node_type A `character` indicating what types of nodes
#' (`site`, `species`, or `both`) should be returned in the output
#' (`return_node_type = "both"` by default).
#' 
#' @param algorithm_in_output A `boolean` indicating if the original output
#' of [cluster_walktrap][igraph::cluster_walktrap] should be returned in the 
#' output (`TRUE` by default, see Value).
#' 
#' @return
#' A `list` of class `bioregion.clusters` with five slots:
#' \enumerate{
#' \item{**name**: A `character` containing the name of the algorithm.}
#' \item{**args**: A `list` of input arguments as provided by the user.}
#' \item{**inputs**: A `list` of characteristics of the clustering process.}
#' \item{**algorithm**: A `list` of all objects associated with the
#'  clustering procedure, such as original cluster objects (only if
#'  `algorithm_in_output = TRUE`).}
#' \item{**clusters**: A `data.frame` containing the clustering results.}}
#'
#' In the `algorithm` slot, if `algorithm_in_output = TRUE`, users can
#' find the output of
#' [cluster_walktrap][igraph::cluster_walktrap].
#'
#' @details
#' This function is based on random walks (Pons & Latapy, 2005)
#' as implemented in the [igraph](https://cran.r-project.org/package=igraph)
#' package ([cluster_walktrap][igraph::cluster_walktrap]).
#'
#' @note
#' Although this algorithm was not primarily designed to deal with bipartite
#' networks, it is possible to consider the bipartite network as unipartite
#' network (`bipartite = TRUE`).
#'
#' Do not forget to indicate which of the first two columns is
#' dedicated to the site nodes (i.e. primary nodes) and species nodes (i.e.
#' feature nodes) using the arguments `site_col` and `species_col`.
#' The type of nodes returned in the output can be chosen with the argument
#' `return_node_type` equal to `both` to keep both types of nodes,
#' `sites` to preserve only the site nodes, and `species` to
#' preserve only the species nodes.

#' 
#' @references 
#' Pons P & Latapy M (2005) Computing Communities in Large Networks 
#' Using Random Walks. In Yolum I, Güngör T, Gürgen F, Özturan C (eds.), 
#' \emph{Computer and Information Sciences - ISCIS 2005}, Lecture Notes in 
#' Computer Science, 284-293.
#' 
#' @seealso 
#' For more details illustrated with a practical example, 
#' see the vignette: 
#' \url{https://biorgeo.github.io/bioregion/articles/a4_3_network_clustering.html}.
#' 
#' Associated functions: 
#' [netclu_infomap] [netclu_louvain] [netclu_oslom]
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
#' com <- netclu_walktrap(net)
#' 
#' net_bip <- mat_to_net(comat, weight = TRUE)
#' clust2 <- netclu_walktrap(net_bip, bipartite = TRUE)
#' 
#' @importFrom igraph graph_from_data_frame cluster_walktrap
#' 
#' @export
netclu_walktrap <- function(net,
                            weight = TRUE,
                            cut_weight = 0,
                            index = names(net)[3],
                            steps = 4,
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
      stop(paste0("Please choose return_node_type from the following:\n",
                  "both, sites or species."), 
           call. = FALSE) 
    }
  }
  
  # Control input loop or directed
  controls(args = NULL, data = net, type = "input_net_isloop")
  controls(args = NULL, data = net, type = "input_net_isdirected")
  
  # Control algorithm_in_output
  controls(args = algorithm_in_output, data = NULL, type = "boolean")
  
  # Control parameters WALKTRAPS
  controls(args = steps, data = NULL, type = "strict_positive_integer")
  
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
  outputs <- list(name = "netclu_walktrap")
  
  outputs$args <- list(
    weight = weight,
    cut_weight = cut_weight,
    index = index,
    steps = steps,
    bipartite = bipartite,
    site_col = site_col,
    species_col = species_col,
    return_node_type = return_node_type,
    algorithm_in_output = algorithm_in_output)
  
  outputs$inputs <- list(
    bipartite = isbip,
    weight = weight,
    pairwise = ifelse(isbip, FALSE, TRUE),
    pairwise_metric = ifelse(!isbip & weight, 
                             ifelse(is.numeric(index), names(net)[3], index), 
                             NA),
    dissimilarity = FALSE,
    nb_sites = nbsites,
    hierarchical = FALSE)
  
  outputs$algorithm <- list()
  
  # Run algo
  net <- igraph::graph_from_data_frame(netemp, directed = FALSE)
  outalg <- igraph::cluster_walktrap(net, steps = steps)
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
