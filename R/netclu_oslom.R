#' OsloM Community Finding
#'
#' This function finds communities in a (un)weighted (un)directed network based
#' on the OsloM algorithm (<http://oslom.org/>, version 2.4).
#'
#' @param net The output object from [similarity()] or
#' [dissimilarity_to_similarity()].
#' If a `data.frame` is used, the first two columns represent pairs of
#' sites (or any pair of nodes), and the next column(s) are the similarity
#' indices.
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
#' @param seed For the random number generator (NULL for random by default).
#'
#' @param reassign A `character` indicating if the nodes belonging to several
#' community should be reassigned and what method should be used (see Note).
#'
#' @param r The number of runs for the first hierarchical level
#' (10 by default).
#'
#' @param hr The number of runs for the higher hierarchical level (50 by
#' default, 0 if you are not interested in hierarchies).
#'
#' @param t The p-value, the default value is 0.10. Increase this value if you want
#' more modules.
#'
#' @param cp Kind of resolution parameter used to decide between taking some
#' modules or their union (default value is 0.5; a bigger value leads to bigger
#' clusters).
#'
#' @param directed A `boolean` indicating if the network is directed (from
#' column 1 to column 2).
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
#' @param binpath A `character` indicating the path to the bin folder
#' (see [install_binaries] and Details).
#'
#' @param check_install A `boolean` indicating if the function should check that
#' the OSLOM has been properly installed (see [install_binaries] and Details).
#'
#' @param path_temp A `character` indicating the path to the temporary folder
#' (see Details).
#'
#' @param delete_temp A `boolean` indicating if the temporary folder should
#' be removed (see Details).
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
#' In the `algorithm` slot, users can find the following elements:
#'
#' \itemize{
#' \item{`cmd`: The command line used to run OSLOM.}
#' \item{`version`: The OSLOM version.}
#' \item{`web`: The OSLOM's web site.}
#' }
#'
#' @details
#' OSLOM is a network community detection algorithm proposed in
#' Lancichinetti et al. (2011) that finds statistically significant
#' (overlapping) communities in (un)weighted and (un)directed networks.
#'
#' This function is based on the 2.4 C++ version of OSLOM
#' (<http://www.oslom.org/software.htm>). This function needs files
#' to run. They can be installed with [install_binaries]. 
#' 
#' **If you changed the default path to the `bin` folder
#' while running [install_binaries], PLEASE MAKE SURE to set `binpath` 
#' accordingly.**
#' 
#' **If you did not use [install_binaries] to change the permissions and test 
#' the binary files, PLEASE MAKE SURE to set `check_install` accordingly.**
#'
#' The C++ version of OSLOM generates temporary folders and/or files that are
#' stored in the `path_temp` folder (folder "oslom_temp" with a unique timestamp
#' located in the bin folder in `binpath` by default). This temporary folder is 
#' removed by default (`delete_temp = TRUE`).
#'
#' @note
#' Although this algorithm was not primarily designed to deal with bipartite
#' networks, it is possible to consider the bipartite network as unipartite
#' network (`bipartite = TRUE`). Do not forget to indicate which of the
#' first two columns is dedicated to the site nodes (i.e. primary nodes) and
#' species nodes (i.e. feature nodes) using the arguments `site_col` and
#' `species_col`. The type of nodes returned in the output can be chosen
#' with the argument `return_node_type` equal to `both` to keep both
#' types of nodes, `sites` to preserve only the sites nodes, and
#' `species` to preserve only the species nodes.
#'
#' Since OSLOM potentially returns overlapping communities, we propose two
#' methods to reassign the 'overlapping' nodes: randomly (`reassign = "random"`)
#' or based on the closest candidate community (`reassign = "simil"`) (only for
#' weighted networks, in this case the closest candidate community is
#' determined with the average similarity). By default, `reassign = "no"` and
#' all the information will be provided. The number of partitions will depend
#' on the number of overlapping modules (up to three). The suffix `_semel`,
#' `_bis`, and `_ter` are added to the column names. The first partition
#' (`_semel`) assigns a module to each node. A value of `NA` in the second
#' (`_bis`) and third (`_ter`) columns indicates that no overlapping module
#' was found for this node (i.e. non-overlapping nodes).
#' 
#' @references 
#' Lancichinetti A, Radicchi F, Ramasco JJ & Fortunato S (2011) Finding 
#' statistically significant communities in networks. \emph{PLOS ONE} 6, 
#' e18961.
#' 
#' @seealso 
#' For more details illustrated with a practical example, 
#' see the vignette: 
#' \url{https://biorgeo.github.io/bioregion/articles/a4_3_network_clustering.html}.
#' 
#' Associated functions: 
#' [netclu_greedy] [netclu_infomap] [netclu_louvain]
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
#' com <- netclu_oslom(net)
#'
#' @export
netclu_oslom <- function(net,
                         weight = TRUE,
                         cut_weight = 0,
                         index = names(net)[3],
                         seed = NULL,
                         reassign = "no",
                         r = 10,
                         hr = 50,
                         t = 0.1,
                         cp = 0.5,
                         directed = FALSE,
                         bipartite = FALSE,
                         site_col = 1,
                         species_col = 2,
                         return_node_type = "both",
                         binpath = "tempdir",
                         check_install = TRUE,
                         path_temp = "oslom_temp",
                         delete_temp = TRUE) {
  
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

  # Check if OSLOM has successfully been installed
  check <- FALSE
  controls(args = directed, data = NULL, type = "boolean")
  if (!directed) {
    if (check_install &
        !file.exists(paste0(binpath, "/bin/OSLOM/check.txt"))) {
      message(paste0("OSLOM is not installed... Please have a look at ",
                     "https://bioRgeo.github.io/bioregion/articles/a1_install_binary_files.html ",
                     "for more details.\n",
                     "It should be located in ",
                     binpath, 
                     "/bin/OSLOM/"))
    } else {
      check <- TRUE
    }
  } else {
    if (check_install &
        !file.exists(paste0(binpath, "/bin/OSLOM/check.txt"))) {
      message(paste0("OSLOM is not installed... Please have a look at ",
                     "https://bioRgeo.github.io/bioregion/articles/a3_1_install_binary_files.html ",
                     "for more details."))
    } else {
      if (!file.exists(paste0(binpath, "/bin/OSLOM/checkdir.txt"))) {
        message(paste0("The directed version of OSLOM is not installed... ", 
                       " Please have a look at ",
                       "https://bioRgeo.github.io/bioregion/articles/a3_1_install_binary_files.html ",
                       "for more details"))
      } else {
        check <- TRUE
      }
    }
  }

  if (check) {
    
    # Control parameters OSLOM
    controls(args = reassign, data = NULL, type = "character")
    if (!(reassign %in% c("no", "random", "simil"))) {
      stop(paste0("Please choose reassign from the following:\n",
                  "no, random or simil."), 
           call. = FALSE)   
    }
    controls(args = r, data = NULL, type = "strict_positive_integer")
    controls(args = hr, data = NULL, type = "positive_integer")
    if(!is.null(seed)){
      controls(args = seed, data = NULL, type = "strict_positive_integer")
    }
    controls(args = t, data = NULL, type = "strict_positive_numeric")
    if (t >= 1) {
      stop("t must be in the interval (0,1)!", call. = FALSE)
    }
    controls(args = cp, data = NULL, type = "strict_positive_numeric")
    if (cp >= 1) {
      stop("cp must be in the interval (0,1)!", call. = FALSE)
    }
    
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
    if (reassign == "simil" & !weight) {
      stop(paste0("A reassignement based on similarity should not be ",
                  "use when weight equal FALSE"))
    }
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

    # Control input directed
    if (!isbip) {
      controls(args = NULL, data = net, type = "input_net_isloop")
      controls(args = directed, data = net, type = "input_net_directed")
    } else {
      if (directed) {
        stop("directed cannot be set to TRUE if the network is bipartite!",
          call. = FALSE
        )
      }
    }
    
    # Control temp folder
    old_path_temp <- path_temp
    if (path_temp == "oslom_temp") {
      fold_temp <- paste0(path_temp,
                          "_",
                          round(as.numeric(as.POSIXct(Sys.time()))))
      path_temp <- paste0(
        binpath,
        "/bin/",
        fold_temp
      )
    } else {
      if (dir.exists(path_temp)) {
        stop(paste0(path_temp, 
                    " already exists. Please rename ", 
                    "it or remove it."),
          call. = FALSE
        )
      }
    }
    path_temp <- normalizePath(path_temp, mustWork = FALSE)
    dir.create(path_temp, showWarnings = FALSE, recursive = TRUE)
    if (!dir.exists(path_temp)) {
      stop(paste0("Impossible to create directory ", path_temp), call. = FALSE)
    }

    # Prepare input for OSLOM
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
    }

    # Class preparation
    outputs <- list(name = "netclu_oslom")

    outputs$args <- list(
      weight = weight,
      cut_weight = cut_weight,
      index = index,
      seed = seed,
      reassign = reassign,
      r = r,
      hr = hr,
      t = t,
      cp = cp,
      directed = directed,
      bipartite = bipartite,
      site_col = site_col,
      species_col = species_col,
      return_node_type = return_node_type,
      binpath = binpath,
      check_install = check_install,
      delete_temp = delete_temp,
      path_temp = path_temp
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

    # Export input for OSLOM
    utils::write.table(netemp, paste0(path_temp, "/net.txt"),
      row.names = FALSE,
      col.names = FALSE, sep = " "
    )

    # Prepare command to run OSLOM
    if(is.null(seed)){
      cmd <- paste0(
        "-r ", r, 
        " -hr ", hr, 
        " -t ", t, 
        " -cp ",
        cp
      )
    }else{
      cmd <- paste0(
        "-r ", r, 
        " -hr ", hr, 
        " -seed ", seed, 
        " -t ", t, 
        " -cp ",
        cp
      )
    }

    # Run OSLOM
    if (os == "Linux") {
      if (weight) {
        cmd <- paste0("-f ", path_temp, "/net.txt -w ", cmd)
      } else {
        cmd <- paste0("-f ", path_temp, "/net.txt -uw ", cmd)
      }
      if (directed) {
        cmd <- paste0(
          binpath, "/bin/OSLOM/oslom_dir_lin ", cmd,
          " > /dev/null 2>&1"
        )
      } else {
        cmd <- paste0(
          binpath, "/bin/OSLOM/oslom_undir_lin ", cmd,
          " > /dev/null 2>&1"
        )
      }
      system(command = cmd)
    } else if (os == "Windows") {
      if (weight) {
        cmd <- paste0("-f ", path_temp, "/net.txt -w ", cmd)
      } else {
        cmd <- paste0("-f ", path_temp, "/net.txt -uw ", cmd)
      }
      if (directed) {
        cmd <- paste0(binpath, "/bin/OSLOM/oslom_dir_win.exe ", cmd)
      } else {
        cmd <- paste0(binpath, "/bin/OSLOM/oslom_undir_win.exe ", cmd)
      }
      dir.create(paste0(path_temp, "/net.txt_oslo_files"),
        showWarnings = FALSE, recursive = TRUE
      )
      system(command = cmd, show.output.on.console = FALSE)
    } else if (os == "Darwin") {
      if(old_path_temp == "oslom_temp"){
        if (weight) {
          cmd <- paste0("-f ", fold_temp, "/net.txt -w ", cmd)
        } else {
          cmd <- paste0("-f ", fold_temp, "/net.txt -uw ", cmd)
        }
        if (directed) {
          cmd1 <- paste0("cd ", binpath, "/bin >/dev/null 2>&1")
          cmd2 <- paste0("OSLOM/oslom_dir_mac ", cmd, " > /dev/null 2>&1")
          cmd <- paste0(cmd1, " && ", cmd2)
        } else {
          cmd1 <- paste0("cd ", binpath, "/bin >/dev/null 2>&1")
          cmd2 <- paste0("OSLOM/oslom_undir_mac ", cmd, " > /dev/null 2>&1")
          cmd <- paste0(cmd1, " && ", cmd2)
        }
      }else{
        if (weight) {
          cmd <- paste0("-f ", path_temp, "/net.txt -w ", cmd)
        } else {
          cmd <- paste0("-f ", path_temp, "/net.txt -uw ", cmd)
        }
        if (directed) {
          cmd <- paste0(
            binpath, "/bin/OSLOM/oslom_dir_mac ", cmd,
            " > /dev/null 2>&1"
          )
        } else {
          cmd <- paste0(
            binpath, "/bin/OSLOM/oslom_undir_mac ", cmd,
            " > /dev/null 2>&1"
          )
        }
      }
      system(command = cmd)
    } else {
      stop("Linux, Windows or Mac distributions only.")
    }

    # Control: if the command line did not work
    if (!("tp" %in% list.files(paste0(path_temp, "/net.txt_oslo_files")))) {
      stop("Command line was wrongly implemented. OSLOM did not run.", 
             call. = FALSE)
    }

    # Number of levels
    nblev <- 1

    # Retrieve output from tp [TO COMMENT]
    com <- readLines(paste0(path_temp, "/net.txt_oslo_files/tp"))
    cl <- list()
    length(cl) <- length(com) / 2
    for (k in 1:length(com)) {
      if ((k / 2 - trunc(k / 2)) == 0) {
        cl[[(k / 2)]] <- as.numeric(as.matrix(strsplit(com[k],
          split = " "
        )[[1]]))
      }
    }

    tab <- unlist(cl)
    tab <- sort(tab[!duplicated(tab)])
    tab <- cbind(tab, 0, 0, 0)
    n <- nrow(tab)

    dupl <- rep(0, n)
    for (i in 1:length(cl)) {
      temp <- rep(0, n)
      temp[match(cl[[i]], tab[, 1])] <- i
      dupl <- dupl + 1 * (temp > 0)
      tab[match(cl[[i]], tab[, 1]), 2] <- i
    }

    for (i in 1:n) {
      if (dupl[i] > 1) {
        overcom <- NULL
        for (j in 1:length(cl)) {
          if (sum(cl[[j]] == tab[i, 1]) == 1) {
            overcom <- c(overcom, j)
          }
        }
        if (reassign == "random") {
          overcom <- overcom[sample(length(overcom), length(overcom))]
        }
        tab[i, 2] <- overcom[1]
        tab[i, 3] <- overcom[2]
        if (dupl[i] == 3) {
          tab[i, 4] <- overcom[3]
        }
      }
    }

    # Reassign tp [TO COMMENT]
    if (reassign == "simil") {
      dat <- netemp
      for (i in 1:n) {
        if (tab[i, 3] > 0) {
          test1 <- sum(dat[, 1] == tab[i, 1])
          if (test1 == 0) {
            dati1 <- matrix(0, nrow = 2, ncol = 2)
          } else if (test1 == 1) {
            dati1 <- matrix(0, nrow = 2, ncol = 2)
            dati1[1, ] <- as.numeric(dat[dat[, 1] == tab[i, 1], c(2, 3)])
          } else {
            dati1 <- dat[dat[, 1] == tab[i, 1], c(2, 3)]
          }
          colnames(dati1) <- c("ID", "SIM")

          test2 <- sum(dat[, 2] == tab[i, 1])
          if (test2 == 0) {
            dati2 <- matrix(0, nrow = 2, ncol = 2)
          } else if (test2 == 1) {
            dati2 <- matrix(0, nrow = 2, ncol = 2)
            dati2[1, ] <- as.numeric(dat[dat[, 2] == tab[i, 1], c(1, 3)])
          } else {
            dati2 <- dat[dat[, 2] == tab[i, 1], c(1, 3)]
          }
          colnames(dati2) <- c("ID", "SIM")

          dati <- rbind(dati1, dati2)

          check <- match(cl[[tab[i, 2]]], dati[, 1])
          check <- check[!is.na(check)]
          sim1 <- mean(dati[check, 2])

          check <- match(cl[[tab[i, 3]]], dati[, 1])
          check <- check[!is.na(check)]
          sim2 <- mean(dati[check, 2])

          if (sim2 > sim1) {
            tab[i, 2] <- tab[i, 3]
            sim1 <- sim2
          }

          if (tab[i, 4] > 0) {
            check <- match(cl[[tab[i, 4]]], dati[, 1])
            check <- check[!is.na(check)]
            sim2 <- mean(dati[check, 2])
            if (sim2 > sim1) {
              tab[i, 2] <- tab[i, 4]
            }
          }
        }
      }
      tabtp <- tab[, 1:2]
    } else if (reassign == "random") {
      tabtp <- tab[, 1:2]
    } else {
      tabtp <- tab
    }

    # Reshape tabtp
    comtp <- data.frame(ID = idnode[, 2], Com1 = 0)
    comtp[match(tabtp[, 1], idnode[, 1]), 2] <- tabtp[, 2]
    if (dim(tabtp)[2] > 2) {
      if (sum(tabtp[, 3]) > 0) {
        comtp$Com2 <- 0
        comtp[match(tabtp[, 1], idnode[, 1]), 3] <- tabtp[, 3]
      }
      if (sum(tabtp[, 4]) > 0) {
        comtp$Com3 <- 0
        comtp[match(tabtp[, 1], idnode[, 1]), 4] <- tabtp[, 4]
      }
    }
    com <- comtp

    # If tp1 exists (i.e. hierarchical level)
    if ("tp1" %in% list.files(paste0(path_temp, "/net.txt_oslo_files"))) {
      # Number of levels
      nblev <- 2

      # Retrieve output from tp1 [TO COMMENT]
      com <- readLines(paste0(path_temp, "/net.txt_oslo_files/tp1"))
      cl <- list()
      length(cl) <- length(com) / 2
      for (k in 1:length(com)) {
        if ((k / 2 - trunc(k / 2)) == 0) {
          cl[[(k / 2)]] <- as.numeric(as.matrix(strsplit(com[k],
            split = " "
          )[[1]]))
        }
      }

      tab <- unlist(cl)
      tab <- sort(tab[!duplicated(tab)])
      tab <- cbind(tab, 0, 0, 0)
      n <- nrow(tab)

      dupl <- rep(0, n)
      for (i in 1:length(cl)) {
        temp <- rep(0, n)
        temp[match(cl[[i]], tab[, 1])] <- i
        dupl <- dupl + 1 * (temp > 0)
        tab[match(cl[[i]], tab[, 1]), 2] <- i
      }

      for (i in 1:n) {
        if (dupl[i] > 1) {
          overcom <- NULL
          for (j in 1:length(cl)) {
            if (sum(cl[[j]] == tab[i, 1]) == 1) {
              overcom <- c(overcom, j)
            }
          }
          if (reassign == "random") {
            overcom <- overcom[sample(length(overcom), length(overcom))]
          }
          tab[i, 2] <- overcom[1]
          tab[i, 3] <- overcom[2]
          if (dupl[i] == 3) {
            tab[i, 4] <- overcom[3]
          }
        }
      }

      # Reassign tp1 [TO COMMENT]
      if (reassign == "simil") {
        dat <- netemp
        for (i in 1:n) {
          if (tab[i, 3] > 0) {
            test1 <- sum(dat[, 1] == tab[i, 1])
            if (test1 == 0) {
              dati1 <- matrix(0, nrow = 2, ncol = 2)
            } else if (test1 == 1) {
              dati1 <- matrix(0, nrow = 2, ncol = 2)
              dati1[1, ] <- as.numeric(dat[dat[, 1] == tab[i, 1], c(2, 3)])
            } else {
              dati1 <- dat[dat[, 1] == tab[i, 1], c(2, 3)]
            }
            colnames(dati1) <- c("ID", "SIM")

            test2 <- sum(dat[, 2] == tab[i, 1])
            if (test2 == 0) {
              dati2 <- matrix(0, nrow = 2, ncol = 2)
            } else if (test2 == 1) {
              dati2 <- matrix(0, nrow = 2, ncol = 2)
              dati2[1, ] <- as.numeric(dat[dat[, 2] == tab[i, 1], c(1, 3)])
            } else {
              dati2 <- dat[dat[, 2] == tab[i, 1], c(1, 3)]
            }
            colnames(dati2) <- c("ID", "SIM")

            dati <- rbind(dati1, dati2)

            check <- match(cl[[tab[i, 2]]], dati[, 1])
            check <- check[!is.na(check)]
            sim1 <- mean(dati[check, 2])

            check <- match(cl[[tab[i, 3]]], dati[, 1])
            check <- check[!is.na(check)]
            sim2 <- mean(dati[check, 2])

            if (sim2 > sim1) {
              tab[i, 2] <- tab[i, 3]
              sim1 <- sim2
            }

            if (tab[i, 4] > 0) {
              check <- match(cl[[tab[i, 4]]], dati[, 1])
              check <- check[!is.na(check)]
              sim2 <- mean(dati[check, 2])
              if (sim2 > sim1) {
                tab[i, 2] <- tab[i, 4]
              }
            }
          }
        }
        tabtph <- tab[, 1:2]
      } else if (reassign == "random") {
        tabtph <- tab[, 1:2]
      } else {
        tabtph <- tab
      }

      # Reshape tabtp1
      comtph <- data.frame(ID = idnode[, 2], HCom1 = 0)
      comtph[match(tabtph[, 1], idnode[, 1]), 2] <- tabtph[, 2]
      if (dim(tabtph)[2] > 2) {
        if (sum(tabtph[, 3]) > 0) {
          comtph$HCom2 <- 0
          comtph[match(tabtph[, 1], idnode[, 1]), 3] <- tabtph[, 3]
        }
        if (sum(tabtph[, 4]) > 0) {
          comtph$HCom3 <- 0
          comtph[match(tabtph[, 1], idnode[, 1]), 4] <- tabtph[, 4]
        }
      }

      com <- cbind(comtp, comtph)
      com <- com[, -(dim(comtp)[2] + 1)]
    }

    # If tp2 exists (i.e. hierarchical level)
    if ("tp2" %in% list.files(paste0(path_temp, "/oslomnet.txt_oslo_files"))) {
      # Number of levels
      nblev <- 3

      # Retrieve output from tp2 [TO COMMENT]
      com <- readLines(paste0(path_temp, "/net.txt_oslo_files/tp2"))
      cl <- list()
      length(cl) <- length(com) / 2
      for (k in 1:length(com)) {
        if ((k / 2 - trunc(k / 2)) == 0) {
          cl[[(k / 2)]] <- as.numeric(as.matrix(strsplit(com[k],
            split = " "
          )[[1]]))
        }
      }

      tab <- unlist(cl)
      tab <- sort(tab[!duplicated(tab)])
      tab <- cbind(tab, 0, 0, 0)
      n <- nrow(tab)

      dupl <- rep(0, n)
      for (i in 1:length(cl)) {
        temp <- rep(0, n)
        temp[match(cl[[i]], tab[, 1])] <- i
        dupl <- dupl + 1 * (temp > 0)
        tab[match(cl[[i]], tab[, 1]), 2] <- i
      }

      for (i in 1:n) {
        if (dupl[i] > 1) {
          overcom <- NULL
          for (j in 1:length(cl)) {
            if (sum(cl[[j]] == tab[i, 1]) == 1) {
              overcom <- c(overcom, j)
            }
          }
          if (reassign == "random") {
            overcom <- overcom[sample(length(overcom), length(overcom))]
          }
          tab[i, 2] <- overcom[1]
          tab[i, 3] <- overcom[2]
          if (dupl[i] == 3) {
            tab[i, 4] <- overcom[3]
          }
        }
      }

      # Reassign tp2 [TO COMMENT]
      if (reassign == "simil") {
        dat <- netemp
        for (i in 1:n) {
          if (tab[i, 3] > 0) {
            test1 <- sum(dat[, 1] == tab[i, 1])
            if (test1 == 0) {
              dati1 <- matrix(0, nrow = 2, ncol = 2)
            } else if (test1 == 1) {
              dati1 <- matrix(0, nrow = 2, ncol = 2)
              dati1[1, ] <- as.numeric(dat[dat[, 1] == tab[i, 1], c(2, 3)])
            } else {
              dati1 <- dat[dat[, 1] == tab[i, 1], c(2, 3)]
            }
            colnames(dati1) <- c("ID", "SIM")

            test2 <- sum(dat[, 2] == tab[i, 1])
            if (test2 == 0) {
              dati2 <- matrix(0, nrow = 2, ncol = 2)
            } else if (test2 == 1) {
              dati2 <- matrix(0, nrow = 2, ncol = 2)
              dati2[1, ] <- as.numeric(dat[dat[, 2] == tab[i, 1], c(1, 3)])
            } else {
              dati2 <- dat[dat[, 2] == tab[i, 1], c(1, 3)]
            }
            colnames(dati2) <- c("ID", "SIM")

            dati <- rbind(dati1, dati2)

            check <- match(cl[[tab[i, 2]]], dati[, 1])
            check <- check[!is.na(check)]
            sim1 <- mean(dati[check, 2])

            check <- match(cl[[tab[i, 3]]], dati[, 1])
            check <- check[!is.na(check)]
            sim2 <- mean(dati[check, 2])

            if (sim2 > sim1) {
              tab[i, 2] <- tab[i, 3]
              sim1 <- sim2
            }

            if (tab[i, 4] > 0) {
              check <- match(cl[[tab[i, 4]]], dati[, 1])
              check <- check[!is.na(check)]
              sim2 <- mean(dati[check, 2])
              if (sim2 > sim1) {
                tab[i, 2] <- tab[i, 4]
              }
            }
          }
        }
        tabtphh <- tab[, 1:2]
      } else if (reassign == "random") {
        tabtphh <- tab[, 1:2]
      } else {
        tabtphh <- tab
      }

      # Reshape tabtp2
      comtphh <- data.frame(ID = idnode[, 2], HHCom1 = 0)
      comtphh[match(tabtphh[, 1], idnode[, 1]), 2] <- tabtphh[, 2]
      if (dim(tabtphh)[2] > 2) {
        if (sum(tabtphh[, 3]) > 0) {
          comtphh$HHCom2 <- 0
          comtphh[match(tabtphh[, 1], idnode[, 1]), 3] <- tabtphh[, 3]
        }
        if (sum(tabtphh[, 4]) > 0) {
          comtphh$HHCom3 <- 0
          comtphh[match(tabtphh[, 1], idnode[, 1]), 4] <- tabtphh[, 4]
        }
      }

      com <- cbind(comtp, comtph, comtphh)
      com <- com[, -c(
        (dim(comtp)[2] + 1),
        (dim(comtp)[2] + dim(comtph)[2] + 1)
      )]
    }

    # If tp3 exists (i.e. hierarchical level)
    if ("tp3" %in% list.files(paste0(path_temp, "/net.txt_oslo_files"))) {
      # Number of levels
      nblev <- 4

      # Retrieve output from tp3 [TO COMMENT]
      com <- readLines(paste0(path_temp, "/net.txt_oslo_files/tp3"))
      cl <- list()
      length(cl) <- length(com) / 2
      for (k in 1:length(com)) {
        if ((k / 2 - trunc(k / 2)) == 0) {
          cl[[(k / 2)]] <- as.numeric(as.matrix(strsplit(com[k],
            split = " "
          )[[1]]))
        }
      }

      tab <- unlist(cl)
      tab <- sort(tab[!duplicated(tab)])
      tab <- cbind(tab, 0, 0, 0)
      n <- nrow(tab)

      dupl <- rep(0, n)
      for (i in 1:length(cl)) {
        temp <- rep(0, n)
        temp[match(cl[[i]], tab[, 1])] <- i
        dupl <- dupl + 1 * (temp > 0)
        tab[match(cl[[i]], tab[, 1]), 2] <- i
      }

      for (i in 1:n) {
        if (dupl[i] > 1) {
          overcom <- NULL
          for (j in 1:length(cl)) {
            if (sum(cl[[j]] == tab[i, 1]) == 1) {
              overcom <- c(overcom, j)
            }
          }
          if (reassign == "random") {
            overcom <- overcom[sample(length(overcom), length(overcom))]
          }
          tab[i, 2] <- overcom[1]
          tab[i, 3] <- overcom[2]
          if (dupl[i] == 3) {
            tab[i, 4] <- overcom[3]
          }
        }
      }

      # Reassign tp3 [TO COMMENT]
      if (reassign == "simil") {
        dat <- netemp
        for (i in 1:n) {
          if (tab[i, 3] > 0) {
            test1 <- sum(dat[, 1] == tab[i, 1])
            if (test1 == 0) {
              dati1 <- matrix(0, nrow = 2, ncol = 2)
            } else if (test1 == 1) {
              dati1 <- matrix(0, nrow = 2, ncol = 2)
              dati1[1, ] <- as.numeric(dat[dat[, 1] == tab[i, 1], c(2, 3)])
            } else {
              dati1 <- dat[dat[, 1] == tab[i, 1], c(2, 3)]
            }
            colnames(dati1) <- c("ID", "SIM")

            test2 <- sum(dat[, 2] == tab[i, 1])
            if (test2 == 0) {
              dati2 <- matrix(0, nrow = 2, ncol = 2)
            } else if (test2 == 1) {
              dati2 <- matrix(0, nrow = 2, ncol = 2)
              dati2[1, ] <- as.numeric(dat[dat[, 2] == tab[i, 1], c(1, 3)])
            } else {
              dati2 <- dat[dat[, 2] == tab[i, 1], c(1, 3)]
            }
            colnames(dati2) <- c("ID", "SIM")

            dati <- rbind(dati1, dati2)

            check <- match(cl[[tab[i, 2]]], dati[, 1])
            check <- check[!is.na(check)]
            sim1 <- mean(dati[check, 2])

            check <- match(cl[[tab[i, 3]]], dati[, 1])
            check <- check[!is.na(check)]
            sim2 <- mean(dati[check, 2])

            if (sim2 > sim1) {
              tab[i, 2] <- tab[i, 3]
              sim1 <- sim2
            }

            if (tab[i, 4] > 0) {
              check <- match(cl[[tab[i, 4]]], dati[, 1])
              check <- check[!is.na(check)]
              sim2 <- mean(dati[check, 2])
              if (sim2 > sim1) {
                tab[i, 2] <- tab[i, 4]
              }
            }
          }
        }
        tabtphhh <- tab[, 1:2]
      } else if (reassign == "random") {
        tabtphhh <- tab[, 1:2]
      } else {
        tabtphhh <- tab
      }

      # Reshape tabtp3
      comtphhh <- data.frame(ID = idnode[, 2], HHHCom1 = 0)
      comtphhh[match(tabtphhh[, 1], idnode[, 1]), 2] <- tabtphhh[, 2]
      if (dim(tabtphhh)[2] > 2) {
        if (sum(tabtphhh[, 3]) > 0) {
          comtphhh$HHHCom2 <- 0
          comtphhh[match(tabtphhh[, 1], idnode[, 1]), 3] <- tabtphhh[, 3]
        }
        if (sum(tabtphhh[, 4]) > 0) {
          comtphhh$HHHCom3 <- 0
          comtphhh[match(tabtphhh[, 1], idnode[, 1]), 4] <- tabtphhh[, 4]
        }
      }

      com <- cbind(comtp, comtph, comtphh, comtphhh)
      com <- com[, -c(
        (dim(comtp)[2] + 1),
        (dim(comtp)[2] + dim(comtph)[2] + 1),
        (dim(comtp)[2] + dim(comtph)[2] + dim(comtphh) + 1)
      )]
    }

    # Remove temporary files
    if (delete_temp) {
      unlink(paste0(path_temp), recursive = TRUE)
    }
    unlink("tp")
    unlink("time_seed.dat")

    # Rename and reorder columns
    com <- as.character(comtp[, 1])
    tempnblev <- 1
    if (nblev >= 4) {
      nov <- dim(comtphhh)[2] - 1
      if (nov == 1) {
        colnames(comtphhh)[2] <- paste0("K_", max(comtphhh[, 2]))
        com <- cbind(com, comtphhh[, 2])
        colnames(com)[dim(com)[2]] <- paste0("K_", max(comtphhh))
      } else {
        colnames(comtphhh)[2] <- paste0("K_", max(comtphhh[, 2]), "_semel")
        if (nov >= 2) {
          colnames(comtphhh)[3] <- paste0("K_", max(comtphhh[, 2]), "_bis")
        }
        if (nov == 3) {
          colnames(comtphhh)[4] <- paste0("K_", max(comtphhh[, 2]), "_ter")
        }
        com <- cbind(com, comtphhh[, -1])
      }
      tempnblev <- tempnblev + 1
    }
    if (nblev >= 3) {
      nov <- dim(comtphh)[2] - 1
      if (nov == 1) {
        colnames(comtphh)[2] <- paste0("K_", max(comtphh[, 2]))
        com <- cbind(com, comtphh[, 2])
        colnames(com)[dim(com)[2]] <- paste0("K_", max(comtphh[, 2]))
      } else {
        colnames(comtphh)[2] <- paste0("K_", max(comtphh[, 2]), "_semel")
        if (nov >= 2) {
          colnames(comtphh)[3] <- paste0("K_", max(comtphh[, 2]), "_bis")
        }
        if (nov == 3) {
          colnames(comtphh)[4] <- paste0("K_", max(comtphh[, 2]), "_ter")
        }
        com <- cbind(com, comtphh[, -1])
      }
      tempnblev <- tempnblev + 1
    }
    if (nblev >= 2) {
      nov <- dim(comtph)[2] - 1
      if (nov == 1) {
        colnames(comtph)[2] <- paste0("K_", max(comtph[, 2]))
        com <- cbind(com, comtph[, 2])
        colnames(com)[dim(com)[2]] <- paste0("K_", max(comtph[, 2]))
      } else {
        colnames(comtph)[2] <- paste0("K_", max(comtph[, 2]), "_semel")
        if (nov >= 2) {
          colnames(comtph)[3] <- paste0("K_", max(comtph[, 2]), "_bis")
        }
        if (nov == 3) {
          colnames(comtph)[4] <- paste0("K_", max(comtph[, 2]), "_ter")
        }
        com <- cbind(com, comtph[, -1])
      }
      tempnblev <- tempnblev + 1
    }
    if (nblev >= 1) {
      nov <- dim(comtp)[2] - 1
      if (nov == 1) {
        colnames(comtp)[2] <- paste0("K_", max(comtp[, 2]))
        com <- cbind(com, comtp[, 2])
        colnames(com)[dim(com)[2]] <- paste0("K_", max(comtp[, 2]))
      } else {
        colnames(comtp)[2] <- paste0("K_", max(comtp[, 2]), "_semel")
        if (nov >= 2) {
          colnames(comtp)[3] <- paste0("K_", max(comtp[, 2]), "_bis")
        }
        if (nov == 3) {
          colnames(comtp)[4] <- paste0("K_", max(comtp[, 2]), "_ter")
        }
        com <- cbind(com, comtp[, -1])
      }
      tempnblev <- tempnblev + 1
    }
    colnames(com)[1] <- "ID"
    com <- data.frame(com)
    for (k in 2:dim(com)[2]) {
      com[, k] <- as.numeric(as.character(com[, k]))
    }

    com[, 1] <- as.character(com[, 1])
    com[,-1][com[,-1]==0]=NA

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
    outputs$algorithm$cmd <- cmd
    outputs$algorithm$version <- "2.4"
    outputs$algorithm$web <- "http://oslom.org/"

    # Set clusters and cluster_info in output
    outputs$clusters <- com
    if(reassign == "no"){
      outputs$cluster_info <- data.frame(
        partition_name = names(outputs$clusters)[2:length(outputs$clusters),
                                                 drop = FALSE])
      outputs$cluster_info$n_clust <- as.numeric(do.call(rbind, 
                         strsplit(outputs$cluster_info$partition_name, 
                                  split="_"))[,2])
    }else{
      outputs$cluster_info <- data.frame(
        partition_name = names(outputs$clusters)[2:length(outputs$clusters),
                                                 drop = FALSE
        ],
        n_clust = apply(
          outputs$clusters[, 2:length(outputs$clusters), drop = FALSE],
          2, function(x) length(unique(x[!is.na(x)]))
        )
      )
    }  


    if (nblev>1) {
      outputs$inputs$hierarchical <- TRUE
      if(reassign == "no"){
        num1=outputs$cluster_info$n_clust
        num2=num1[!duplicated(num1)]
        outputs$cluster_info$hierarchical_level <- match(num1,num2)
      }else{
        outputs$cluster_info$hierarchical_level <- 1:nrow(outputs$cluster_info)
      }
    }
    
    # Return outputs
    class(outputs) <- append("bioregion.clusters", class(outputs))
    return(outputs)
  }
}
