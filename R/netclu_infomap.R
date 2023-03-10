#' Infomap community finding
#'
#' This function finds communities in a (un)weighted (un)directed network based
#' on the Infomap algorithm (<https://github.com/mapequation/infomap>).
#'
#' @param net the output object from [similarity()] or
#' [dissimilarity_to_similarity()].
#' If a `data.frame` is used, the first two columns represent pairs of
#' sites (or any pair of nodes), and the next column(s) are the similarity
#' indices.
#' 
#' @param weight a `boolean` indicating if the weights should be considered
#' if there are more than two columns.
#' 
#' @param index name or number of the column to use as weight. By default,
#' the third column name of `net` is used.
#' 
#' @param nbmod penalize solutions the more they differ from this number (0 by
#' default for no preferred number of modules).
#' 
#' @param markovtime scales link flow to change the cost of moving between
#' modules, higher values results in fewer modules (default is 1).
#' 
#' @param seed for the random number generator (0 for random by default).
#' 
#' @param numtrials for the number of trials before picking up the best
#' solution.
#' 
#' @param twolevel a `boolean` indicating if the algorithm should optimize a
#' two-level partition of the network (default is multi-level).
#' 
#' @param show_hierarchy a `boolean` specifying if the hierarchy of community
#' should be identifiable in the outputs (FALSE by default).
#' 
#' @param directed a `boolean` indicating if the network is directed (from
#' column 1 to column 2).
#' 
#' @param bipartite a `boolean` indicating if the network is bipartite
#' (see Note).
#' 
#' @param bipartite_version a `boolean` indicating if the bipartite version of
#' Infomap should be used (see Note).
#' 
#' @param site_col name or number for the column of site nodes 
#' (i.e. primary nodes).
#' 
#' @param species_col name or number for the column of species nodes
#' (i.e. feature nodes).
#' 
#' @param return_node_type a `character` indicating what types of nodes
#' ("sites", "species" or "both") should be returned in the output
#' (`keep_nodes_type="both"` by default).
#' 
#' @param version a `character` indicating the Infomap version to use.
#' 
#' @param delete_temp a `boolean` indicating if the temporary folder should
#' be removed (see Details).
#' 
#' @param path_temp a `character` indicating the path to the temporary folder
#' (see Details).
#' 
#' @param binpath a `character` indicating the path to the bin folder
#' (see [install_binaries] and Details).
#' 
#' @details
#' Infomap is a network clustering algorithm based on the Map equation proposed
#' in \insertCite{Rosvall2008}{bioregion} that finds communities in (un)weighted
#' and (un)directed networks.
#'
#' This function is based on the C++ version of Infomap
#' (<https://github.com/mapequation/infomap/releases>).
#' This function needs binary files to run. They can be installed
#' with [install_binaries()]. If you set the path to the folder that will host
#' the bin folder manually while running [install_binaries] please make sure to
#' set `binpath` accordingly.
#'
#' The C++ version of Infomap generates temporary folders and/or files that are
#' stored in the `path_temp` folder ("infomap_temp" with an unique time stamp
#' located in the working directory by default). This temporary folder is
#' removed by default (`delete_temp = TRUE`).
#'
#' Several version of Infomap are available in the package. See
#' [install_binaries()] for more details.
#'
#' @note
#' Infomap has been designed to deal with bipartite networks. To use this
#' functionality set the `bipartite_version` argument to TRUE in order to
#' approximate a two-step random walker (see
#' <https://www.mapequation.org/infomap/> for more information). Note that
#' a bipartite network can also be considered as unipartite network
#' (`bipartite = TRUE`).
#'
#' In both cases do not forget to indicate which of the first two columns is
#' dedicated to the site nodes (i.e. primary nodes) and species nodes (i.e.
#' feature nodes) using the arguments `site_col` and `species_col`.
#' The type of nodes returned in the output can be chosen with the argument
#' `return_node_type` equal to `"both"` to keep both types of nodes, `"sites"`
#' to preserve only the sites nodes and `"species"` to preserve only the
#' species nodes.
#'
#' @return
#' A `list` of class `bioregion.clusters` with five slots:
#' \enumerate{
#' \item{**name**: `character string` containing the name of the algorithm}
#' \item{**args**: `list` of input arguments as provided by the user}
#' \item{**inputs**: `list` of characteristics of the input dataset}
#' \item{**algorithm**: `list` of all objects associated with the
#'  clustering procedure, such as original cluster objects}
#' \item{**clusters**: `data.frame` containing the clustering results}}
#'
#' In the `algorithm` slot, users can find the following elements:
#'
#' \itemize{
#' \item{`cmd`: the command line use to run Infomap}
#' \item{`version`: the Infomap version}
#' \item{`web`: Infomap's GitHub repository}
#' }
#'
#' @author
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}),
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) and
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' 
#' @seealso [install_binaries], [netclu_louvain], [netclu_oslom]
#' 
#' @examples
#' \dontrun{
#' comat <- matrix(sample(1000, 50), 5, 10)
#' rownames(comat) <- paste0("Site", 1:5)
#' colnames(comat) <- paste0("Species", 1:10)
#'
#' net <- similarity(comat, metric = "Simpson")
#' com <- netclu_infomap(net)
#' }
#' 
#' @references
#' \insertRef{Rosvall2008}{bioregion}
#' 
#' @export
netclu_infomap <- function(net,
                           weight = TRUE,
                           index = names(net)[3],
                           nbmod = 0,
                           markovtime = 1,
                           seed = 0,
                           numtrials = 1,
                           twolevel = FALSE,
                           show_hierarchy = FALSE,
                           directed = FALSE,
                           bipartite_version = FALSE,
                           bipartite = FALSE,
                           site_col = 1,
                           species_col = 2,
                           return_node_type = "both",
                           version = "2.6.0",
                           delete_temp = TRUE,
                           path_temp = "infomap_temp",
                           binpath = NULL) {
  
  # Remove warning for tidyr
  defaultW <- getOption("warn")
  options(warn = -1)
  
  # Set binpath
  if (is.null(binpath)) {
    # Identify bioregion directory on your computer
    biodir <- .libPaths()[1]
    binpath <- paste0(biodir,"/bioregion")
  } else {
    # Control
    controls(args = binpath, data = NULL, type = "character")
    if (!dir.exists(binpath)) {
      stop(paste0("Impossible to access ", binpath), call. = FALSE)
    }
  }
  
  # Control version
  controls(args = version, data = NULL, type = "character")
  # versiondispo <- list.files(paste0(binpath, "/bin/INFOMAP/"))
  # if (!(version %in% versiondispo)) {
  #  stop(paste0(
  #    "Please choose a version of Infomap already installed: ",
  #    paste(versiondispo, collapse = " ")
  #  ), call. = FALSE)
  # }
  
  # Check OS
  os <- Sys.info()[["sysname"]]
  
  # Check if INFOMAP has successfully been installed
  if (!file.exists(paste0(binpath, "/bin/INFOMAP/", version, "/check.txt"))) {
    message(paste0(
      "Infomap ", version, " is not installed... Please have a look at
    https//bioregion.github.io/bioregion/articles/a1_install_binary_files.html
    for more details.\n", 
      "It should be located in ", 
      paste0(binpath, "/bin/INFOMAP/", version, "/")))
  } else {
    
    # Control input net
    controls(args = NULL, data = net, type = "input_bioregion.pairwise.metric")
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
    controls(args = bipartite_version, data = NULL, type = "boolean")
    controls(args = bipartite, data = NULL, type = "boolean")
    isbip <- (bipartite | bipartite_version)
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
    if (!isbip) {
      controls(args = directed, data = net, type = "input_net_directed")
    } else {
      if (directed) {
        stop("directed cannot be set to TRUE if the network is bipartite!",
             call. = FALSE
        )
      }
    }
    
    # Control parameters INFOMAP
    controls(args = nbmod, data = NULL, type = "positive_integer")
    controls(args = markovtime, data = NULL, type = "strict_positive_numeric")
    controls(args = seed, data = NULL, type = "positive_integer")
    if (seed == 0) {
      seed <- round(as.numeric(as.POSIXct(Sys.time())))
    }
    controls(args = numtrials, data = NULL, type = "strict_positive_integer")
    controls(args = twolevel, data = NULL, type = "boolean")
    controls(args = show_hierarchy, data = NULL, type = "boolean")
    
    
    # Control temp folder + create temp folder
    path_temp <- controls(args = path_temp, data = NULL, type = "character")
    controls(args = delete_temp, data = NULL, type = "boolean")
    if (path_temp == "infomap_temp") {
      path_temp <- paste0(path_temp, "_", round(as.numeric(as.POSIXct(Sys.time()))))
    } else {
      if (dir.exists(path_temp)) {
        stop(paste0(path_temp, " already exists. Please rename it or remove it."),
             call. = FALSE
        )
      }
    }
    dir.create(path_temp, showWarnings = FALSE, recursive = TRUE)
    if (!dir.exists(path_temp)) {
      stop(paste0("Impossible to create directory ", path_temp), call. = FALSE)
    }
    
    # Prepare input for INFOMAP
    if (isbip) {
      idprim <- as.character(net[, site_col])
      idprim <- idprim[!duplicated(idprim)]
      nbsites <- length(idprim)
      idfeat <- as.character(net[, species_col])
      idfeat <- idfeat[!duplicated(idfeat)]
      
      idprim <- data.frame(
        ID = 1:length(idprim), # Primary nodes
        ID_NODE = idprim,
        Type = 1
      )
      idfeat <- data.frame(
        ID = ((dim(idprim)[1] + 1):(dim(idprim)[1] # Feature nodes
                                    + length(idfeat))),
        ID_NODE = idfeat, Type = 2
      )
      N <- dim(idprim)[1] + 1 # First node id of the feature node type
      idnode <- rbind(idprim, idfeat)
      if (!bipartite_version) {
        idnode <- idnode[, 1:2]
      }
      netemp <- data.frame(
        node1 = idnode[match(net[, site_col], idnode[, 2]), 1],
        node2 = idnode[match(net[, species_col], idnode[, 2]), 1])
    } else {
      idnode1 <- as.character(net[, 1])
      idnode2 <- as.character(net[, 2])
      if (isbip) {
        message("The network seems to be bipartite! 
The bipartite or bipartite_version argument should probably be set to TRUE.")
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
    }
    
    # Class preparation
    outputs <- list(name = "netclu_infomap")
    
    outputs$args <- list(weight = weight,
                         index = index,
                         nbmod = nbmod,
                         markovtime = markovtime,
                         seed = seed,
                         numtrials = numtrials,
                         twolevel = twolevel,
                         directed = directed,
                         bipartite_version = bipartite_version,
                         bipartite = bipartite,
                         site_col = site_col,
                         species_col = species_col,
                         return_node_type = return_node_type,
                         version = version,
                         delete_temp = delete_temp,
                         path_temp = path_temp,
                         binpath = binpath)
    
    outputs$inputs <- list(bipartite = isbip,
                           weight = weight,
                           pairwise = ifelse(isbip, FALSE, TRUE),
                           pairwise_metric = ifelse(isbip, NA, index),
                           dissimilarity = FALSE,
                           nb_sites = nbsites)
    
    outputs$algorithm <- list()
    
    # Export input in INFOMAP folder
    if (bipartite_version) { # Add tag if bipartite
      cat(paste0("*Bipartite ", N), "\n", file = paste0(path_temp, "/net.txt"))
      utils::write.table(netemp, paste0(path_temp, "/net.txt"),
                         append = TRUE,
                         row.names = FALSE, col.names = FALSE, sep = " ")
    } else {
      utils::write.table(netemp, paste0(path_temp, "/net.txt"),
                         row.names = FALSE,
                         col.names = FALSE, sep = " ")
    }
    
    # Prepare command to run INFOMAP
    cmd <- paste0("--silent --seed ", seed,
                  " --num-trials ", numtrials,
                  " --markov-time ", markovtime)
    if (nbmod > 0) {
      cmd <- paste0(cmd, " --preferred-number-of-modules ", nbmod)
    }
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
      cmd <- paste0(binpath, "/bin/INFOMAP/", version, "/infomap_lin ", cmd)
    } else if (os == "Windows") {
      cmd <- paste0(binpath, "/bin/INFOMAP/", version, "/infomap_win.exe ", cmd)
    } else if (os == "Darwin") {
      cmd <- paste0(binpath, "/bin/INFOMAP/", version, "/infomap_mac ", cmd)
    } else {
      stop("Linux, Windows or Mac distributions only.", call. = FALSE)
    }
    
    # Run INFOMAP
    system(command = cmd)
    
    # Control: if the command line did not work
    if (!("net.tree" %in% list.files(paste0(path_temp)))) {
      stop("Command line was wrongly implemented. Infomap did not run.",
           call. = FALSE)
    }
    
    # Retrieve output from net.tree
    tree <- utils::read.table(paste0(path_temp, "/net.tree"))
    
    # Import tree
    idinf <- as.numeric(tree[, 4]) # INFOMAP node ids
    tree <- as.character(tree[, 1]) # INFOMAP tree column
    
    # Extract community from tree
    cominf <- reformat_hierarchy(tree, algo = "infomap",
                                 integerize = !show_hierarchy)
    
    cominf <- cominf[, -1]
    nblev <- dim(cominf)[2]
    cominf <- cominf[, nblev:1] # Reverse column order
    
    com <- data.frame(ID = idnode[, 2], dum = 0) # Dummy level
    com[match(idinf, idnode[, 1]), 2] <- cominf[, 1]
    
    for (knblev in 2:nblev) {
      com$temp <- 0
      com[match(idinf, idnode[, 1]), (knblev + 1)] <- cominf[, knblev]
      colnames(com)[(knblev + 1)] <- paste0("V", (knblev + 1))
    }
    
    com <- com[, -2] # Remove dummy level
    
    # Remove temp folder
    if (delete_temp) {
      unlink(paste0(path_temp), recursive = TRUE)
    }
    
    # Rename and reorder columns
    com <- knbclu(com)
    
    # Add attributes and return_node_type
    if (isbip) {
      attr(com, "node_type") <- rep("site", dim(com)[1])
      attributes(com)$node_type[!is.na(match(com[, 1],
                                             idfeat$ID_NODE))] <- "species"
      if (return_node_type == "sites") {
        com <- com[attributes(com)$node_type == "site", ]
      }
      if (return_node_type == "species") {
        com <- com[attributes(com)$node_type == "species", ]
      }
    }
    
    # Put the warning back
    options(warn = defaultW)
    
    # Set algorithm in outputs
    outputs$algorithm$cmd <- cmd
    outputs$algorithm$version <- version
    outputs$algorithm$web <- "https://github.com/mapequation/infomap"
    
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
    if (!twolevel & show_hierarchy) {
      outputs$cluster_info$hierarchical_level <- 1:nrow(outputs$cluster_info)
    }
    
    # Return outputs
    class(outputs) <- append("bioregion.clusters", class(outputs))
    return(outputs)
  }
}
