#' Community structure detection in weighted bipartite network via modularity optimisation
#'
#' This function takes a bipartite weighted graph and computes modules by applying Newman’s modularity
#' measure in a bipartite weighted version to it.
#'
#' @param net a two- or three-column \code{data.frame} representing a network with the two first columns
#' as undirected links between pair of nodes and an optional third column indicating the weight of the link
#' @param weight a boolean indicating if the weights should be considered if there is a third column
#' @param forceLPA a boolean indicating if the even faster pure LPA-algorithm of Beckett should be used? DIRT-
#' LPA, the default, is less likely to get trapped in a local minimum, but is slightly
#' slower. Defaults to FALSE.
#' @param primary_col name or number for the column of primary nodes (i.e. site)
#' @param feature_col name or number for the column of feature nodes (i.e. species)
#' @param remove_feature a boolean indicating if the feature nodes should be removed from the outputs (TRUE by default)
#' @export
#' @details
#' This function is based on the modularity optimization algorithm provided by Stephen Beckett \insertCite{Beckett2016}{bioRgeo}
#' as implemented in the \href{https://cran.r-project.org/web/packages/bipartite/index.html}{bipartite} package
#' (\link[bipartite]{computeModules}).
#'
#' @return A \code{data.frame} providing one community by node.
#'
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}),
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) and
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' @seealso \link{netclu_infomap}, \link{netclu_oslom}
#' @examples
#' net <- data.frame(
#'   Site = c(rep("A", 2), rep("B", 3), rep("C", 2)),
#'   Species = c("a", "b", "a", "c", "d", "b", "d"),
#'   Weight = c(10, 100, 1, 20, 50, 10, 20)
#' )
#'
#' com <- netclu_beckett(net)
#' @references
#' \insertRef{Beckett2016}{bioRgeo}
#' @export
netclu_beckett <- function(net, weight = TRUE, forceLPA = FALSE,
                           primary_col = 1, feature_col = 2, remove_feature = TRUE) {

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

  if (!is.logical(forceLPA)) {
    stop("forceLPA must be a boolean")
  }

  # Controls bipartite arguments
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
  idprim <- as.character(net[, primary_col])
  idprim <- idprim[!duplicated(idprim)]
  idfeat <- as.character(net[, feature_col])
  idfeat <- idfeat[!duplicated(idfeat)]

  if (length(intersect(idprim, idfeat)) > 0) { # Control bipartite
    stop("The network should be bipartite!")
  }

  idnode <- c(idprim, idfeat)
  idnode <- data.frame(ID = 1:length(idnode), ID_NODE = idnode)

  netemp <- data.frame(node1 = idnode[match(net[, primary_col], idnode[, 2]), 1], node2 = idnode[match(net[, feature_col], idnode[, 2]), 1])
  if (weight) {
    netemp <- cbind(netemp, net[, 3])
    netemp <- netemp[netemp[, 3] > 0, ]
    colnames(netemp)[3] <- "weight"
  }

  # Transform netemp into a contingency table
  comat <- net_to_mat(netemp, weight = weight)

  # Run algo
  comtemp <- bipartite::computeModules(comat, forceLPA = forceLPA)@modules[-1, -c(1, 2)]
  comtemp <- cbind(c(as.numeric(rownames(comat)), as.numeric(colnames(comat))), apply(comtemp, 2, function(x) which(x > 0)))

  com <- data.frame(ID = idnode[, 2], Com = 0)
  com[match(comtemp[, 1], idnode[, 1]), 2] <- comtemp[, 2]

  # Remove feature nodes
  if (remove_feature) {
    com <- com[match(idprim, com[, 1]), ]
  }

  # Rename and reorder columns
  com[, 1] <- as.character(com[, 1])
  com <- knbclu(com)

  # Return output
  return(com)
}
