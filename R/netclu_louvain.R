#' Louvain community finding
#'
#' This function finds communities in a (un)weighted undirected network based
#' on the Louvain algorithm.
#'
#' @param net the output object from [similarity()] or
#' [dissimilarity_to_similarity()].
#' If a `data.frame` is used, the first two columns represent pairs of sites
#' (or any pair of nodes), and the next column(s) are the similarity indices.
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
#' @param lang a string indicating what version of Louvain should be used
#' (`igraph` or `cpp`, see Details).
#' 
#' @param resolution a resolution parameter to adjust the modularity 
#' (1 is chosen by default, see Details).
#' 
#' @param seed for the random number generator (only when `lang = "igraph"`, 
#' NULL for random by default).
#'
#' @param q the quality function used to compute partition of the graph
#' (modularity is chosen by default, see Details).
#'
#' @param c the parameter for the Owsinski-Zadrozny quality function
#' (between 0 and 1, 0.5 is chosen by default).
#'
#' @param k the kappa_min value for the Shi-Malik quality function
#' (it must be > 0, 1 is chosen by default).
#'
#' @param bipartite a boolean indicating if the network is bipartite
#' (see Details).
#'
#' @param site_col name or number for the column of site nodes
#' (i.e. primary nodes).
#'
#' @param species_col name or number for the column of species nodes
#' (i.e. feature nodes).
#'
#' @param return_node_type a `character` indicating what types of nodes
#' (`site`, `species` or `both`) should be returned in the output
#' (`return_node_type = "both"` by default).
#'
#' @param binpath a `character` indicating the path to the bin folder
#' (see [install_binaries] and Details).
#' 
#' @param check_install a `boolean` indicating if the function should check that
#' the Louvain as been properly installed (see [install_binaries] and Details).
#'
#' @param path_temp a `character` indicating the path to the temporary folder
#' (see Details).
#'
#' @param delete_temp a `boolean` indicating if the temporary folder should
#' be removed (see Details).
#'
#' @param algorithm_in_output a `boolean` indicating if the original output
#' of [cluster_louvain][igraph::cluster_louvain] should be returned in the 
#' output (`TRUE` by default, see Value). 
#'
#' @details
#' Louvain is a network community detection algorithm proposed in
#' \insertCite{Blondel2008}{bioregion}. This function proposed two
#' implementations of the function (parameter `lang`): the
#' [igraph](https://cran.r-project.org/package=igraph)
#' implementation ([cluster_louvain][igraph::cluster_louvain]) and the C++
#' implementation (<https://sourceforge.net/projects/louvain/>, version 0.3).
#' 
#' The [igraph](https://cran.r-project.org/package=igraph)
#' implementation offers the possibility to adjust the resolution parameter of 
#' the modularity function (`resolution` argument) that the algorithm uses 
#' internally. Lower values typically yield fewer, larger clusters. The original
#' definition of modularity is recovered when the resolution parameter 
#' is set to 1 (by default).
#' 
#' The C++ implementation offers the possibility to choose among several 
#' quality functions,
#' `q = 0` for the classical Newman-Girvan criterion (also called
#' "Modularity"), 1 for the Zahn-Condorcet criterion, 2 for the
#' Owsinski-Zadrozny criterion (you should specify the value of the parameter
#' with the `c` argument), 3 for the Goldberg Density criterion, 4 for the
#' A-weighted Condorcet criterion, 5 for the Deviation to Indetermination
#' criterion, 6 for the Deviation to Uniformity criterion, 7 for the Profile
#' Difference criterion, 8 for the Shi-Malik criterion (you should specify the
#' value of kappa_min with `k` argument) and 9 for the Balanced Modularity
#' criterion.
#'
#' The C++ version of Louvain is based on the version 0.3
#' (<https://sourceforge.net/projects/louvain/>). This function needs
#'  binary files to run. They can be installed with
#' [install_binaries]. 
#' 
#' **If you changed the default path to the `bin` folder
#' while running [install_binaries] PLEASE MAKE SURE to set `binpath` 
#' accordingly.**
#' 
#' **If you did not used [install_binaries] to change the permissions and test 
#' the binary files PLEASE MAKE SURE to set `check_install` accordingly.**
#' 
#' The C++ version of Louvain generates temporary folders and/or files that are
#' stored in the `path_temp` folder ("louvain_temp" with an unique timestamp
#' located in the bin folder in `binpath` by default). This temporary folder
#' is removed by default (`delete_temp = TRUE`).
#'
#' @note
#' Although this algorithm was not primarily designed to deal with bipartite
#' network, it is possible to consider the bipartite network as unipartite
#' network (`bipartite = TRUE`).
#'
#' Do not forget to indicate which of the first two columns is dedicated to the
#' site nodes (i.e. primary nodes) and species nodes (i.e. feature nodes) using
#' the arguments `site_col` and `species_col`. The type of nodes returned in
#' the output can be chosen with the argument `return_node_type` equal to
#' `both` to keep both types of nodes, `sites` to preserve only the sites
#' nodes and `species` to preserve only the species nodes.
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
#' In the `algorithm` slot, if `algorithm_in_output = TRUE`, users can find an
#' the output of [cluster_louvain][igraph::cluster_louvain]
#' if `lang = "igraph"` and the following element if `lang = "cpp"`:
#'
#' \itemize{
#' \item{`cmd`: the command line use to run Louvain}
#' \item{`version`: the Louvain version}
#' \item{`web`: Louvain's website}
#' }.
#'
#' @author
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}),
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) and
#' Boris Leroy (\email{leroy.boris@gmail.com})
#'
#' @seealso [install_binaries()], [netclu_infomap()], [netclu_oslom()]
#'
#' @examples
#' comat <- matrix(sample(1000, 50), 5, 10)
#' rownames(comat) <- paste0("Site", 1:5)
#' colnames(comat) <- paste0("Species", 1:10)
#'
#' net <- similarity(comat, metric = "Simpson")
#' com <- netclu_louvain(net, lang = "igraph")
#'
#' @references
#' \insertRef{Blondel2008}{bioregion}
#'
#' @importFrom igraph graph_from_data_frame cluster_louvain
#'
#' @export

netclu_louvain <- function(net,
                           weight = TRUE,
                           cut_weight = 0,
                           index = names(net)[3],
                           lang = "igraph",
                           resolution = 1,
                           seed = NULL,
                           q = 0,
                           c = 0.5,
                           k = 1,
                           bipartite = FALSE,
                           site_col = 1,
                           species_col = 2,
                           return_node_type = "both",
                           binpath = "tempdir",
                           check_install = TRUE,
                           path_temp = "louvain_temp",
                           delete_temp = TRUE,
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

  # Control parameters LOUVAIN
  controls(args = lang, data = NULL, type = "character")
  if (!(lang %in% c("cpp", "igraph"))) {
    stop("Please choose lang among the following values:
cpp or igraph", call. = FALSE)
  }
  controls(args = resolution, data = NULL, type = "strict_positive_numeric")
  if(!is.null(seed)){
    controls(args = seed, data = NULL, type = "strict_positive_integer")
  }
  controls(args = q, data = NULL, type = "positive_integer")
  controls(args = c, data = NULL, type = "strict_positive_numeric")
  if (c >= 1) {
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
  outputs <- list(name = "netclu_louvain")

  outputs$args <- list(
    weight = weight,
    cut_weight = cut_weight,
    index = index,
    lang = lang,
    resolution = resolution,
    seed = seed,
    q = q,
    c = c,
    k = k,
    bipartite = bipartite,
    site_col = site_col,
    species_col = species_col,
    return_node_type = return_node_type,
    binpath = binpath,
    check_install = check_install,
    delete_temp = delete_temp,
    path_temp = path_temp,
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

  # igraph
  if (lang == "igraph") {
    
    # Run algo (with seed)
    net <- igraph::graph_from_data_frame(netemp, directed = FALSE)
    if(is.null(seed)){
      outalg <- igraph::cluster_louvain(net, resolution = resolution)
    }else{
      set.seed(seed)
      outalg <- igraph::cluster_louvain(net, resolution = resolution)
      rm(.Random.seed, envir=globalenv())
    }
    comtemp <- cbind(as.numeric(outalg$names), as.numeric(outalg$membership))

    com <- data.frame(ID = idnode[, 2], Com = NA)
    com[match(comtemp[, 1], idnode[, 1]), 2] <- comtemp[, 2]

    # Set algorithm in outputs
    if (!algorithm_in_output) {
      outalg <- NA
    }
    outputs$algorithm <- outalg
  }

  # cpp
  if (lang == "cpp") {
    
    # Control empty network
    if(dim(netemp)[1]==0){
      stop("The network is empty. 
         Please check your data or choose an appropriate cut_weight value.")
    }
    
    # Control and set binpath
    controls(args = binpath, data = NULL, type = "character")
    controls(args = check_install, data = NULL, type = "boolean")
    controls(args = path_temp, data = NULL, type = "character")
    controls(args = delete_temp, data = NULL, type = "boolean")
    if (binpath == "tempdir") {
      binpath <- tempdir()
    } else if (binpath == "pkgfolder") {
      binpath <- paste0(.libPaths()[1], "/bioregion")
    } else {
      if (!dir.exists(binpath)) {
        stop(paste0("Impossible to access ", binpath), call. = FALSE)
      }
    }
    binpath <- normalizePath(binpath)

    # Check OS
    os <- Sys.info()[["sysname"]]

    # Check if LOUVAIN has successfully been installed
    if (check_install &
        !file.exists(paste0(binpath, "/bin/LOUVAIN/check.txt"))) {
      message(
        "Louvain is not installed... Please have a look at
              https://bioRgeo.github.io/bioregion/articles/a1_install_binary_files.html
              for more details.\n",
        "It should be located in ",
        paste0(binpath, "/bin/LOUVAIN/")
      )
    } else {
      # Control temp folder + create temp folder
      if (path_temp == "louvain_temp") {
        path_temp <- paste0(
          binpath,
          "/bin/",
          path_temp,
          "_",
          round(as.numeric(as.POSIXct(Sys.time())))
        )
      } else {
        if (dir.exists(path_temp)) {
          stop(paste0(path_temp, " already exists. Please rename it or remove
                      it."),
            call. = FALSE
          )
        }
      }
      path_temp <- normalizePath(path_temp, mustWork = FALSE)
      dir.create(path_temp, showWarnings = FALSE, recursive = TRUE)
      if (!dir.exists(path_temp)) {
        stop(paste0("Impossible to create directory ", path_temp),
          call. = FALSE
        )
      }
      
      # Reclassify nodes 
      idnode1b <- as.character(netemp[, 1])
      idnode2b <- as.character(netemp[, 2])
      idnodeb <- c(idnode1b, idnode2b)
      idnodeb <- idnodeb[!duplicated(idnodeb)]
      idnodeb <- data.frame(IDb = 1:length(idnodeb), ID_NODEb = idnodeb)
      netemp[,1] <- idnodeb[match(netemp[,1],idnodeb[,2]),1]
      netemp[,2] <- idnodeb[match(netemp[,2],idnodeb[,2]),1]

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
        stop("Command line was wrongly implemented. Louvain did not run.", 
             call. = FALSE)
      }

      # Retrieve output from net.tree
      tree <- utils::read.table(paste0(path_temp, "/net.tree"))
      
      # Retrieve hierarchy
      tree <- reformat_hierarchy(tree, 
                                 algo = "louvain")
      
      tree[,1] <- idnodeb[match(tree[,1],idnodeb[,1]),2]
    
      com <- data.frame(ID = idnode[, 2], Com = NA)
      com[match(tree[, 1], idnode[, 1]), 2] <- tree[, 2]
      if(dim(tree)[2]>2){
        for (k in 3:dim(tree)[2]) {
          com$temp <- NA
          com[match(tree[,1], idnode[, 1]), k] <- tree[, k]
          colnames(com)[k] <- paste0("V", k)
        }
      }

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
    if (return_node_type == "site") {
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
      2, function(x) length(unique(x[!is.na(x)]))
    )
  )
  
  if (nrow(outputs$cluster_info)>1) {
    outputs$cluster_info$hierarchical_level <- 1:nrow(outputs$cluster_info)
    outputs$inputs$hierarchical <- TRUE
  }

  # Return outputs
  class(outputs) <- append("bioregion.clusters", class(outputs))
  return(outputs)
}
