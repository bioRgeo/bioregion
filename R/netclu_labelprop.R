#' Finding communities based on propagating labels
#'
#' This function finds communities in a (un)weighted undirected network based on propagating labels.
#'
#' @param net a two- or three-column \code{data.frame} representing a network with the two first columns
#' as (un)directed links between pair of nodes and an optional third column indicating the weight of the link
#' @param weight a boolean indicating if the weights should be considered if there is a third column
#' @param bipartite a boolean indicating if the network is bipartite (see Details)
#' @param primary_col name or number for the column of primary nodes (i.e. site)
#' @param feature_col name or number for the column of feature nodes (i.e. species)
#' @param remove_feature a boolean indicating if the feature nodes should be removed from the outputs (TRUE by default)
#'
#' @export
#' @details
#' This function is based on propagating labels \insertCite{Raghavan2007}{bioRgeo}
#' as implemented in the \href{https://cran.r-project.org/web/packages/igraph/index.html}{igraph} package
#' (\link[igraph]{cluster_label_prop}).
#'
#' Although this algorithm was not primarily designed to deal with bipartite network, it is possible to consider
#' the bipartite network as unipartite network by using the arguments \code{bipartite}, \code{primary_col},
#' \code{feature_col} and \code{remove_feature}).
#'
#' @return A \code{data.frame} providing one community by node.
#'
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}),
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) and
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' @seealso \link{netclu_infomap}, \link{netclu_oslom}
#' @examples
#' comat <- matrix(sample(1000, 50), 5, 10)
#' rownames(comat) <- paste0("Site", 1:5)
#' colnames(comat) <- paste0("Species", 1:10)
#'
#' net <- similarity(comat, metric = "Simpson")
#' com <- netclu_labelprop(net)
#' @references
#' \insertRef{Raghavan2007}{bioRgeo}
#' @export
netclu_labelprop <- function(net, weight = TRUE,
                             bipartite = FALSE, primary_col = 1, feature_col = 2, remove_feature = TRUE) {

  # Controls input net
  if (!is.data.frame(net)) {
    stop("net must be a two- or three-columns data.frame")
  }

  if (dim(net)[2] != 2 & dim(net)[2] != 3) {
    stop("net must be a two- or three-columns data.frame")
  }

  sco <- sum(is.na(net))
  if (sco > 0) {
    stop("NA(s) detected in the data.frame")
  }

  # Controls parameters
  if (weight & dim(net)[2] == 2) {
    stop("net must be a three-columns data.frame if weight equal TRUE")
  }

  if (weight & dim(net)[2] == 3) {
    if (!is.numeric(net[, 3])) {
      stop("The third column of net must be numeric")
    }
  }

  if (!is.logical(weight)) {
    stop("weight must be a boolean")
  }

  # Controls bipartite arguments
  if (!is.logical(bipartite)) {
    stop("bipartite must be a boolean")
  }

  if (is.character(primary_col)) {
    if (!(primary_col %in% colnames(net))) {
      stop("primary_col should be a column name")
    }
  } else if (is.numeric(primary_col)) {
    if (primary_col <= 0) {
      stop("primary_col must be strictly positive")
    } else {
      if (primary_col %% 1 != 0) {
        stop("primary_col must be an integer")
      }
    }
  } else {
    stop("primary_col should be numeric or character")
  }

  if (is.character(feature_col)) {
    if (!(feature_col %in% colnames(net))) {
      stop("feature_col should be a column name")
    }
  } else if (is.numeric(feature_col)) {
    if (feature_col <= 0) {
      stop("feature_col must be strictly positive")
    } else {
      if (feature_col %% 1 != 0) {
        stop("feature_col must be an integer")
      }
    }
  } else {
    stop("feature_col should be numeric or character")
  }

  if (!is.logical(remove_feature)) {
    stop("remove_feature must be a boolean")
  }

  # Prepare input
  if (bipartite) {
    idprim <- as.character(net[, primary_col])
    idprim <- idprim[!duplicated(idprim)]
    idfeat <- as.character(net[, feature_col])
    idfeat <- idfeat[!duplicated(idfeat)]
    # Control bipartite
    if (length(intersect(idprim, idfeat)) > 0) {
      stop("If bipartite = TRUE primary and feature nodes should be different.")
    }
    idnode <- c(idprim, idfeat)
    idnode <- data.frame(ID = 1:length(idnode), ID_NODE = idnode)
    netemp <- data.frame(node1 = idnode[match(net[, primary_col], idnode[, 2]), 1], node2 = idnode[match(net[, feature_col], idnode[, 2]), 1])
  } else {
    idnode1 <- as.character(net[, 1])
    idnode2 <- as.character(net[, 2])
    if (length(intersect(idnode1, idnode2)) == 0) {
      stop("The network is bipartite! The bipartite argument should be set to TRUE.")
    }
    idnode <- c(idnode1, idnode2)
    idnode <- idnode[!duplicated(idnode)]
    idnode <- data.frame(ID = 1:length(idnode), ID_NODE = idnode)
    netemp <- data.frame(node1 = idnode[match(net[, 1], idnode[, 2]), 1], node2 = idnode[match(net[, 2], idnode[, 2]), 1])
  }

  if (weight) {
    netemp <- cbind(netemp, net[, 3])
    netemp <- netemp[netemp[, 3] > 0, ]
    colnames(netemp)[3] <- "weight"
  }

  # Run algo
  net <- igraph::graph_from_data_frame(netemp, directed = FALSE)
  comtemp <- igraph::cluster_label_prop(net)
  comtemp <- cbind(as.numeric(comtemp$names), as.numeric(comtemp$membership))

  com <- data.frame(ID = idnode[, 2], Com = 0)
  com[match(comtemp[, 1], idnode[, 1]), 2] <- comtemp[, 2]

  # Remove feature nodes
  if (bipartite & remove_feature) {
    com <- com[match(idprim, com[, 1]), ]
  }

  # Rename and reorder columns
  com[, 1] <- as.character(com[, 1])
  com <- knbclu(com)

  # Return output
  return(com)
  
}
