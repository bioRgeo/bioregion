#' Infomap community finding
#'
#' This function finds communities in a (un)weighted (un)directed network based on the Infomap algorithm
#' (\url{https://github.com/mapequation/infomap}, version 1.6.0).
#'
#' @param net a two- or three-column \code{data.frame} representing a network with the two first columns
#' as (un)directed links between pair of nodes and an optional third column indicating the weight of the link
#' @param weight a boolean indicating if the weights should be considered if there is a third column
#' @param nbmod penalize solutions the more they differ from this number (0 by default for no preferred number of modules)
#' @param markovtime scales link flow to change the cost of moving between modules, higher values results
#' in fewer modules (default is 1)
#' @param seed for the random number generator
#' @param numtrials for the number of trials before picking up the best solution
#' @param twolevel a boolean indicating if the algorithm should optimize a two-level partition of the network
#' (default is multi-level)
#' @param directed a boolean indicating if the network is directed (from column 1 to column 2)
#' @param bipartite a boolean indicating if the network is bipartite (see Details)
#' @param bipartite_version a boolean indicating if the bipartite version of Infomap should be used (see Details)
#' @param primary_col name or number for the column of primary nodes (i.e. site)
#' @param feature_col name or number for the column of feature nodes (i.e. species)
#' @param remove_feature a boolean indicating if the feature nodes should be removed from the outputs (TRUE by default)
#' @param delete_temp a boolean indicating if the temporary folder should be removed (see Details)
#' @param path_temp a string indicating the path to the temporary folder (see Details)
#' @param binpath a string indicating the path to the bin folder (see \link{install_binaries} and Details)
#' @export
#' @details
#' Infomap is a network clustering algorithm based on the Map equation proposed in
#' \insertCite{Rosvall2008}{bioRgeo} that finds communities in (un)weighted and (un)directed networks.
#' Infomap has two ways to deal with bipartite networks. The first possibility is to consider the bipartite network
#' as unipartite network (arguments \code{bipartite}, \code{primary_col}, \code{feature_col} and \code{remove_feature}).
#' The second possibility is to set the \code{bipartite_version} argument to TRUE in order to
#' approximate a two-step random walker (see \url{https://www.mapequation.org/infomap/} for more information).
#'
#' This function is based on the 2.1.0 C++ version of Infomap (\url{https://github.com/mapequation/infomap/releases}).
#' This function needs executable files to run. They can be installed with \link{install_binaries}. If you set the path to
#' the folder that will host the bin folder manually while running \link{install_binaries} please make sure to set \code{binpath}
#' accordingly.
#'
#' The C++ version of Infomap generates temporary folders and/or files that are stored in the \code{path_temp} folder
#' (folder "infomap_temp" in the working directory by default). This temporary folder is removed by default
#' (\code{delete_temp = TRUE}).
#'
#' @return A \code{data.frame} providing one partition by hierarchical level.
#'
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}),
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) and
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' @seealso \link{netclu_louvain}, \link{netclu_oslom}
#' @examples
#' comat <- matrix(sample(1000, 50), 5, 10)
#' rownames(comat) <- paste0("Site", 1:5)
#' colnames(comat) <- paste0("Species", 1:10)
#'
#' net <- similarity(comat, metric = "Simpson")
#' # com=netclu_infomap(net) # run install_binaries() to use this function
#' @references
#' \insertRef{Rosvall2008}{bioRgeo}
#' @export
netclu_infomap <- function(net, weight = TRUE, nbmod = 0, markovtime = 1, seed = 1, numtrials = 1, twolevel = FALSE, directed = FALSE,
                           bipartite_version = FALSE, bipartite = FALSE, primary_col = 1, feature_col = 2,
                           remove_feature = TRUE, delete_temp = TRUE, path_temp = "infomap_temp", binpath = NULL) {

  # Remove warning for tidyr
  defaultW <- getOption("warn")
  options(warn = -1)

  # Set binpath
  if (is.null(binpath)) {
    # Identify bioRgeo directory on your computer
    biodir <- list.dirs(.libPaths(), recursive = FALSE)
    binpath <- biodir[grep("bioRgeo", biodir)]
  } else {
    # Control
    if (!is.character(binpath)) {
      stop("path must be a string")
    }
    if (!file.exists(binpath)) {
      stop(paste0("Impossible to access ", binpath))
    }
  }

  # Check OS
  os <- Sys.info()[["sysname"]]

  # Check if INFOMAP has successfully been installed
  if (!file.exists(paste0(binpath, "/bin/INFOMAP/check.txt"))) {
    stop("Infomap is not installed... Please have a look at https//biorgeo.github.io/bioRgeo/articles/bin.html for more details.")
  }

  # Control input net
  if (!is.data.frame(net)) {
    stop("net must be a two- or three-columns data.frame")
  }

  if (dim(net)[2] != 2 & dim(net)[2] != 3) {
    stop("net must be a two- or three-columns data.frame")
  }

  nbna <- sum(is.na(net))
  if (nbna > 0) {
    stop("NA(s) detected in the data.frame")
  }

  # Control parameters
  if (!is.logical(weight)) {
    stop("weight must be a boolean")
  }

  if (weight & dim(net)[2] == 2) {
    stop("net must be a three-columns data.frame if weight equal TRUE")
  }

  if (weight & dim(net)[2] == 3) {
    if (!is.numeric(net[, 3])) {
      stop("The third column of net must be numeric")
    } else {
      minet <- min(net[, 3])
      if (minet < 0) {
        stop("The third column of net should contains only positive real: negative value detected!")
      }
    }
  }

  if (!is.numeric(nbmod)) {
    stop("nbmod must be numeric")
  } else {
    if (nbmod < 0) {
      stop("nbmod must be positive")
    } else {
      if (nbmod %% 1 != 0) {
        stop("nbmod must be an integer")
      }
    }
  }

  if (!is.numeric(markovtime)) {
    stop("markovtime must be numeric")
  } else {
    if (markovtime <= 0) {
      stop("markovtime must be strictly higher than 0")
    }
  }

  if (!is.numeric(seed)) {
    stop("seed must be numeric")
  } else {
    if (seed <= 0) {
      stop("seed must be strictly higher than 0")
    } else {
      if (seed %% 1 != 0) {
        stop("nbmod must be an integer")
      }
    }
  }

  if (!is.numeric(numtrials)) {
    stop("numtrials must be numeric")
  } else {
    if (numtrials <= 0) {
      stop("numtrials must be strictly higher than 0")
    } else {
      if (numtrials %% 1 != 0) {
        stop("nbmod must be an integer")
      }
    }
  }

  if (!is.logical(twolevel)) {
    stop("twolevel must be a boolean")
  }

  if (!is.logical(directed)) {
    stop("directed must be a boolean")
  }

  if (!is.logical(delete_temp)) {
    stop("delete_temp must be a boolean")
  }

  if (!is.character(path_temp)) {
    stop("path_temp must be a string")
  }

  # Controls bipartite arguments
  if (!is.logical(bipartite)) {
    stop("bipartite must be a boolean")
  }

  if (!is.logical(bipartite_version)) {
    stop("bipartite_version must be a boolean")
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

  # Create temp folder
  dir.create(path_temp, showWarnings = FALSE, recursive = TRUE)
  if (!file.exists(path_temp)) {
    stop(paste0("Impossible to create directory ", path_temp))
  }

  # Prepare input for INFOMAP
  if (bipartite | bipartite_version) {
    idprim <- as.character(net[, primary_col])
    idprim <- idprim[!duplicated(idprim)]
    idfeat <- as.character(net[, feature_col])
    idfeat <- idfeat[!duplicated(idfeat)]
    # Control bipartite
    if (length(intersect(idprim, idfeat)) > 0) {
      stop("If bipartite = TRUE or bipartite_version = TRUE primary and feature nodes should be different.")
    }
    idprim <- data.frame(ID = 1:length(idprim), ID_NODE = idprim, Type = 1) # Primary nodes
    idfeat <- data.frame(ID = ((dim(idprim)[1] + 1):(dim(idprim)[1] + length(idfeat))), ID_NODE = idfeat, Type = 2) # Feature nodes
    N <- dim(idprim)[1] + 1 # First node id of the feature node type
    idnode <- rbind(idprim, idfeat)
    if (!bipartite_version) {
      idnode <- idnode[, 1:2]
    }
    netemp <- data.frame(node1 = idnode[match(net[, primary_col], idnode[, 2]), 1], node2 = idnode[match(net[, feature_col], idnode[, 2]), 1])
  } else {
    idnode1 <- as.character(net[, 1])
    idnode2 <- as.character(net[, 2])
    if (length(intersect(idnode1, idnode2)) == 0) {
      stop("The network is bipartite! The bipartite or bipartite_version argument should be set to TRUE.")
    }
    idnode <- c(idnode1, idnode2)
    idnode <- idnode[!duplicated(idnode)]
    idnode <- data.frame(ID = 1:length(idnode), ID_NODE = idnode)
    netemp <- data.frame(node1 = idnode[match(net[, 1], idnode[, 2]), 1], node2 = idnode[match(net[, 2], idnode[, 2]), 1])
  }

  if (weight) {
    netemp <- cbind(netemp, net[, 3])
    netemp <- netemp[netemp[, 3] > 0, ]
  }

  # Export input in INFOMAP folder
  if (bipartite_version) { # Add tag if bipartite
    cat(paste0("*Bipartite ", N), "\n", file = paste0(path_temp, "/net.txt"))
    utils::write.table(netemp, paste0(path_temp, "/net.txt"), append = TRUE, row.names = FALSE, col.names = FALSE, sep = " ")
  } else {
    utils::write.table(netemp, paste0(path_temp, "/net.txt"), row.names = FALSE, col.names = FALSE, sep = " ")
  }

  # Prepare command to run INFOMAP
  cmd <- paste0("--silent --seed ", seed, "--num-trials", numtrials, " --preferred-number-of-modules ", nbmod, " --markov-time ", markovtime)
  if (twolevel) {
    cmd <- paste0(cmd, " --two-level")
  }
  if (directed) {
    cmd <- paste0(cmd, " --flow-model directed")
  } else {
    cmd <- paste0(cmd, " --flow-model undirected")
  }

  cmd <- paste0(cmd, " ", path_temp, "/net.txt ", path_temp)


  if (os == "Linux") {
    cmd <- paste0(binpath, "/bin/INFOMAP/infomap_lin ", cmd)
  } else if (os == "Windows") {
    cmd <- paste0(binpath, "/bin/INFOMAP/infomap_win.exe ", cmd)
  } else if (os == "Darwin") {
    cmd <- paste0(binpath, "/bin/INFOMAP/infomap_mac ", cmd)
  } else {
    stop("Linux, Windows or Mac distributions only.")
  }

  # Run INFOMAP
  system(command = cmd)

  # Control: if the command line did not work
  if (!("net.tree" %in% list.files(paste0(path_temp)))) {
    stop("Command line was wrongly implemented. Infomap did not run.")
  }

  # Retrieve output from net.tree
  tree <- utils::read.table(paste0(path_temp, "/net.tree"))

  # Reformat tree [TO COMMENT]
  idinf <- as.numeric(tree[, 4]) # INFOMAP node ids

  # require(tidyr) # Extract the modules from tree
  df <- data.frame(x = as.character(tree[, 1]))
  cominf <- tidyr::separate(df, "x", as.character(1:100), sep = ":") # Max 100 levels
  cominf[is.na(cominf)] <- 0
  for (k in 1:dim(cominf)[2]) { # Transform in numeric
    cominf[, k] <- as.numeric(as.character(cominf[, k]))
  }
  cominf <- cominf[, apply(cominf, 2, sum) > 0] # Data frame with information contains in the first column of tree (0 when no info)
  nblev <- dim(cominf)[2] # The number of columns of cominf correspond to the number of levels + one dummy column

  for (k in 2:nblev) { # Set a value 0 for the dummy column
    cominf[cominf[, k] == 0, (k - 1)] <- 0
  }
  cominf[cominf[, nblev] > 0, nblev] <- 0

  for (k in 2:nblev) { # Extract a real partition for lower levels
    cominf[, k] <- as.numeric(factor(paste0(cominf[, k - 1], "_", cominf[, k])))
  }

  cominf <- cominf[, nblev:1] # Reverse column order

  com <- data.frame(ID = idnode[, 2], dum = 0) # Dummy level
  com[match(idinf, idnode[, 1]), 2] <- cominf[, 1]

  if (nblev >= 2) {
    com$Com <- 0
    com[match(idinf, idnode[, 1]), 3] <- cominf[, 2]
  }
  if (nblev >= 3) {
    com$HCom <- 0
    com[match(idinf, idnode[, 1]), 4] <- cominf[, 3]
  }
  if (nblev >= 4) {
    com$HHCom <- 0
    com[match(idinf, idnode[, 1]), 5] <- cominf[, 4]
  }
  if (nblev >= 5) {
    com$HHHCom <- 0
    com[match(idinf, idnode[, 1]), 6] <- cominf[, 5]
  }
  if (nblev >= 6) {
    com$HHHHCom <- 0
    com[match(idinf, idnode[, 1]), 7] <- cominf[, 6]
  }

  com <- com[, -2] # Remove dum

  # Remove temporary
  if (delete_temp) {
    unlink(paste0(path_temp), recursive = TRUE)
  }

  # Remove feature nodes
  if ((bipartite | bipartite_version) & remove_feature) {
    com <- com[match(idprim$ID_NODE, com[, 1]), ]
  }

  # Rename and reorder columns
  com[, 1] <- as.character(com[, 1])
  com <- knbclu(com)

  # Put the warning back
  options(warn = defaultW)

  # Return output
  return(com)
}
