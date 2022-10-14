#' Louvain community finding
#'
#' This function finds communities in a (un)weighted undirected network based on
#' the Louvain algorithm.
#'
#' @param net the output object from \code{\link{similarity}} or
#' \code{\link{dissimilarity_to_similarity}}.
#' If a \code{data.frame} is used, the first two columns represent pairs of
#' sites (or any pair of nodes), and the next column(s) are the similarity
#' indices.
#' @param weight a \code{boolean} indicating if the weights should be considered
#' if there are more than two columns.
#' @param index name or number of the column to use as weight. By default,
#' the third column name of \code{net} is used.
#' @param lang a string indicating what version of Louvain should be used
#' (igraph or Cpp, see Details).
#' @param q the quality function used to compute partition of the graph
#' (modularity is chosen by default, see Details).
#' @param c the parameter for the Owsinski-Zadrozny quality function
#' (between 0 and 1, 0.5 is chosen by default).
#' @param k the kappa_min value for the Shi-Malik quality function
#' (it must be > 0, 1 is chosen by default).
#' @param bipartite a boolean indicating if the network is bipartite
#' (see Details).
#' @param site_col name or number for the column of site nodes
#' (i.e. primary nodes).
#' @param species_col name or number for the column of species nodes
#' (i.e. feature nodes).
#' @param return_node_type a \code{character} indicating what types of nodes
#' ("sites", "species" or "both") should be returned in the output
#' (\code{keep_nodes_type="both"} by default).
#' @param delete_temp a \code{boolean} indicating if the temporary folder should
#' be removed (see Details).
#' @param path_temp a \code{character} indicating the path to the temporary
#' folder (see Details).
#' @param binpath a \code{character} indicating the path to the bin folder
#' (see \link{install_binaries} and Details).
#' @param algorithm_in_output a \code{boolean} indicating if the original output
#' of \code{communities} should be returned in the output (see Value).
#' Default to TRUE.
#' @export
#' @details
#' Louvain is a network community detection algorithm proposed in
#' \insertCite{Blondel2008}{bioRgeo}. This function proposed two implementations
#' of the function (parameter \code{lang}):
#' the \href{https://cran.r-project.org/web/packages/igraph/index.html}{igraph}
#' implementation (\link[igraph]{cluster_louvain}) and the C++ implementation
#' (\url{https://sourceforge.net/projects/louvain/}, version 0.3). The latest
#' offers the possibility to choose among several quality functions,
#' \code{q = 0} for the classical Newman-Girvan criterion (also called
#' "Modularity"), 1 for the Zahn-Condorcet criterion, 2 for the Owsinski-Zadrozny
#' criterion (you should specify the value of the parameter with the \code{c}
#' argument),
#' 3	for the Goldberg Density criterion, 4	for the A-weighted Condorcet criterion,
#' 5 for the Deviation to Indetermination criterion, 6 for the Deviation to
#' Uniformity criterion, 7 for the Profile Difference criterion, 8	for the
#' Shi-Malik criterion (you should specify the value of kappa_min with \code{k}
#'  argument)
#' and 9	for the Balanced Modularity criterion.
#'
#' The C++ version of Louvain is based on the version 0.3
#' (\url{https://sourceforge.net/projects/louvain/}). This function needs
#' executable binary files to run. They can be installed with
#' \link{install_binaries}. If you set the path to the folder that will host the
#' bin folder manually while running \link{install_binaries} please make sure to
#' set \code{binpath} accordingly.
#'
#' The C++ version of Louvain generates temporary folders and/or files that are
#' stored in the \code{path_temp} folder ("louvain_temp" with an unique timestamp
#' located in the working directory by default). This temporary folder is removed
#' by default (\code{delete_temp = TRUE}).
#'
#' @note
#' Although this algorithm was not primarily designed to deal with bipartite
#' network, it is possible to consider the bipartite network as unipartite
#' network (\code{bipartite = TRUE}).
#'
#' Do not forget to indicate which of the first two columns is
#' dedicated to the site nodes (i.e. primary nodes) and species nodes (i.e.
#' feature nodes) using the arguments \code{site_col} and \code{species_col}.
#' The type of nodes returned in the output can be chosen with the argument
#' \code{return_node_type} equal to \code{"both"} to keep both types of nodes,
#' \code{"sites"} to preserve only the sites nodes and \code{"species"} to
#' preserve only the species nodes.
#'
#' @return
#' A \code{list} of class \code{bioRgeo.clusters} with five slots:
#' \enumerate{
#' \item{\bold{name}: \code{character string} containing the name of the algorithm}
#' \item{\bold{args}: \code{list} of input arguments as provided by the user}
#' \item{\bold{inputs}: \code{list} of characteristics of the input dataset}
#' \item{\bold{algorithm}: \code{list} of all objects associated with the
#'  clustering procedure, such as original cluster objects (only if
#'  \code{algorithm_in_output = TRUE})}
#' \item{\bold{clusters}: \code{data.frame} containing the clustering results}}
#'
#' In the \code{algorithm} slot, if \code{algorithm_in_output = TRUE}, users can
#' find an "communities" object, output of \link[igraph]{cluster_louvain}
#' if \code{lang = "igraph"} and the following element if \code{lang = "Cpp"}:
#'
#' \itemize{
#' \item{\code{cmd}: the command line use to run Louvain}
#' \item{\code{version}: the Louvain version}
#' \item{\code{web}: Louvain's website}
#' }.
#'
#' @author
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}),
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) and
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' @seealso \link{install_binaries}, \link{netclu_infomap}, \link{netclu_oslom}
#' @examples
#' comat <- matrix(sample(1000, 50), 5, 10)
#' rownames(comat) <- paste0("Site", 1:5)
#' colnames(comat) <- paste0("Species", 1:10)
#'
#' net <- similarity(comat, metric = "Simpson")
#' com <- netclu_louvain(net, lang = "igraph")
#' @references
#' \insertRef{Blondel2008}{bioRgeo}
#' @export
netclu_louvain <- function(net,
                           weight = TRUE,
                           index = names(net)[3],
                           lang = "Cpp",
                           q = 0,
                           c = 0.5,
                           k = 1,
                           bipartite = FALSE,
                           site_col = 1,
                           species_col = 2,
                           return_node_type = "both",
                           delete_temp = TRUE,
                           path_temp = "louvain_temp",
                           binpath = NULL,
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

  # Control parameters LOUVAIN
  if (!(lang %in% c("Cpp", "igraph"))) {
    stop("Please choose lang among the following values:
Cpp, igraph", call. = FALSE)
  }
  controls(args = q, data = NULL, type = "positive_integer")
  controls(args = c, data = NULL, type = "strict_positive_numeric")
  if (c > 1) {
    stop("c must be in the interval (0,1)!", call. = FALSE)
  }
  controls(args = k, data = NULL, type = "strict_positive_numeric")
  controls(args = algorithm_in_output, data = NULL, type = "boolean")

  # Prepare input for LOUVAIN
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
  outputs <- list(name = "netclu_louvain")

  outputs$args <- list(
    weight = weight,
    index = index,
    lang = lang,
    q = q,
    c = c,
    k = k,
    bipartite = bipartite,
    site_col = site_col,
    species_col = species_col,
    return_node_type = return_node_type,
    delete_temp = delete_temp,
    path_temp = path_temp,
    binpath = binpath,
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

  # igraph
  if (lang == "igraph") {

    # Run algo
    net <- igraph::graph_from_data_frame(netemp, directed = FALSE)
    outalg <- igraph::cluster_louvain(net)
    comtemp <- cbind(as.numeric(outalg$names), as.numeric(outalg$membership))

    com <- data.frame(ID = idnode[, 2], Com = 0)
    com[match(comtemp[, 1], idnode[, 1]), 2] <- comtemp[, 2]

    # Set algorithm in outputs
    if (!algorithm_in_output) {
      outalg <- NA
    }
    outputs$algorithm <- outalg
  }

  # Cpp
  if (lang == "Cpp") {

    # Set binpath
    if (is.null(binpath)) {
      # Identify bioRgeo directory on your computer
      biodir <- list.dirs(.libPaths(), recursive = FALSE)
      binpath <- biodir[grep("bioRgeo", biodir)]
      if (length(binpath) > 1) {
        message("Several bioRgeo directories have been detected in your default package/library folder(s). 
The first one will be used by default.
Please use the binpath argument to manually set the path to the bin folder.")
        binpath <- binpath[1]
      }
    } else {
      # Control
      controls(args = binpath, data = NULL, type = "character")
      if (!file.exists(binpath)) {
        stop(paste0("Impossible to access ", binpath), call. = FALSE)
      }
    }

    # Check OS
    os <- Sys.info()[["sysname"]]

    # Check if LOUVAIN has successfully been installed
    if (!file.exists(paste0(binpath, "/bin/LOUVAIN/check.txt"))) {
      message("Louvain is not installed... 
Please have a look at https://biorgeo.github.io/bioRgeo/articles/a3_1_install_executable_binary_files.html for more details.")
    } else {

      # Control temp folder + create temp folder
      path_temp <- controls(args = path_temp, data = NULL, type = "character")
      controls(args = delete_temp, data = NULL, type = "boolean")
      if (path_temp == "louvain_temp") {
        path_temp <- paste0(path_temp, "_", round(as.numeric(as.POSIXct(Sys.time()))))
      } else {
        if (file.exists(path_temp)) {
          stop(paste0(path_temp, " already exists. Please rename it or remove it."),
            call. = FALSE
          )
        }
      }
      dir.create(path_temp, showWarnings = FALSE, recursive = TRUE)
      if (!file.exists(path_temp)) {
        stop(paste0("Impossible to create directory ", path_temp), call. = FALSE)
      }

      # Export input in LOUVAIN folder
      utils::write.table(netemp, paste0(path_temp, "/net.txt"),
        row.names = FALSE, col.names = FALSE, sep = " "
      )

      # Prepare command to run LOUVAIN
      # Convert net.txt with LOUVAIN
      if (weight) {
        cmd <- paste0(
          "-i ", path_temp, "/net.txt -o ", path_temp, "/net.bin -w ",
          path_temp, "/net.weights"
        )
      } else {
        cmd <- paste0("-i ", path_temp, "/net.txt -o ", path_temp, "/net.bin")
      }

      if (os == "Linux") {
        cmd <- paste0(binpath, "/bin/LOUVAIN/convert_lin ", cmd)
      } else if (os == "Windows") {
        cmd <- paste0(binpath, "/bin/LOUVAIN/convert_win.exe ", cmd)
      } else if (os == "Darwin") {
        cmd <- paste0(binpath, "/bin/LOUVAIN/convert_mac ", cmd)
      } else {
        stop("Linux, Windows or Mac distributions only.")
      }

      tree <- system(command = cmd)

      # Run LOUVAIN
      if (weight) {
        cmd <- paste0(
          path_temp, "/net.bin -l -1 -q ", q, " -c ", c, " -k ", k,
          " -w ", path_temp, "/net.weights"
        )
      } else {
        cmd <- paste0(path_temp, "/net.bin -l -1 -q ", q, " -c ", c, " -k ", k)
      }

      if (os == "Linux") {
        cmd <- paste0(
          binpath, "/bin/LOUVAIN/louvain_lin ", cmd, " > ",
          path_temp, "/net.tree"
        )
        system(command = cmd)
      } else if (os == "Windows") {
        cmd <- paste0(binpath, "/bin/LOUVAIN/louvain_win.exe ", cmd)
        tree <- system(command = cmd, intern = TRUE)
        cat(tree[1:(length(tree) - 1)],
          file = paste0(path_temp, "/net.tree"),
          sep = "\n"
        )
      } else if (os == "Darwin") {
        cmd <- paste0(
          binpath, "/bin/LOUVAIN/louvain_mac ", cmd, " > ",
          path_temp, "/net.tree"
        )
        system(command = cmd)
      } else {
        stop("Linux, Windows or Mac distributions only.")
      }

      # Control: if the command line did not work
      if (!("net.tree" %in% list.files(paste0(path_temp)))) {
        stop("Command line was wrongly implemented. Louvain did not run.")
      }

      # Retrieve output from net.tree
      tree <- utils::read.table(paste0(path_temp, "/net.tree"))

      id0 <- which(tree[, 1] == 0)
      tree <- tree[(id0[1] + 1):(id0[2] - 1), ]

      com <- data.frame(ID = idnode[, 2], Com = 0)
      com[match(tree[, 1], idnode[, 1]), 2] <- tree[, 2]

      # Remove temporary file
      if (delete_temp) {
        unlink(paste0(path_temp), recursive = TRUE)
      }

      # Set algorithm in outputs
      outputs$algorithm$cmd <- cmd
      outputs$algorithm$version <- "0.3"
      outputs$algorithm$web <- "https://sourceforge.net/projects/louvain/"
    }
  }

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

  # Set clusters and cluster_info in output
  outputs$clusters <- com
  outputs$cluster_info <- data.frame(
    partition_name = names(outputs$clusters)[2:length(outputs$clusters),
      drop = FALSE
    ],
    n_clust = apply(
      outputs$clusters[, 2:length(outputs$clusters), drop = FALSE],
      2, function(x) length(unique(x))
    )
  )

  # Return outputs
  class(outputs) <- append("bioRgeo.clusters", class(outputs))
  return(outputs)
}
