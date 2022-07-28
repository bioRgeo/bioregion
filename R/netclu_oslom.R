#' OSLOM community finding
#'
#' This function finds communities in a (un)weighted (un)directed network based on the OSLOM algorithm
#' (\url{http://oslom.org/}, version 2.4).
#'
#' @param net a two- or three-column \code{data.frame} representing a network with the two first columns
#' as (un)directed links between pair of nodes and an optional third column indicating the weight of the link
#' @param weight a boolean indicating if the weights should be considered if there is a third column
#' @param reassign a string indicating if the nodes belonging to several community should be reassign
#' and what method should be used (see Details)
#' @param r the number of runs for the first hierarchical level (10 by default)
#' @param hr the number of runs for the higher hierarchical level (50 by default, 0 if you are not interested
#' in hierarchies)
#' @param seed for the random number generator
#' @param t the p−value, the default value is 0.10, increase this value you to get more modules
#' @param cp kind of resolution parameter used to decide between taking some modules or their union (default value is
#' 0.5, bigger value leads to bigger clusters)
#' @param directed a boolean indicating if the network is directed (from column 1 to column 2)
#' @param bipartite a boolean indicating if the network is bipartite (see Details)
#' @param primary_col name or number for the column of primary nodes (i.e. site)
#' @param feature_col name or number for the column of feature nodes (i.e. species)
#' @param remove_feature a boolean indicating if the feature nodes should be removed from the outputs (TRUE by default)
#' @param delete_temp a boolean indicating if the temporary folder should be removed (see Details)
#' @param path_temp a string indicating the path to the temporary folder (see Details)
#' @param binpath a string indicating the path to the bin folder (see \link{install_binaries} and Details)
#' @export
#' @details
#' OSLOM is a network community detection algorithm proposed in \insertCite{Lancichinetti2011}{bioRgeo} that
#' finds statistically significant (overlapping) communities in (un)weighted and (un)directed networks. Since a node
#' may belong to several communities we propose two methods to reassign the 'overlapping' nodes randomly
#' \code{reassign = 'random'} or based on the closest candidate community \code{reassign = 'simil'}
#' (only for weighted networks, in this case the closest candidate community is determined with the average similarity).
#' By default \code{reassign = 'no'} and all the information is provided.
#'
#' Although this algorithm was not primarily designed to deal with bipartite network, it is possible to consider
#' the bipartite network as unipartite network by using the arguments \code{bipartite}, \code{primary_col},
#' \code{feature_col} and \code{remove_feature}.
#'
#' This function is based on the 2.4 C++ version of OSLOM (\url{http://www.oslom.org/software.htm}).
#' This function needs executable files to run. They can be installed with \link{install_binaries}. If you set the path to
#' the folder that will host the bin folder  manually while running \link{install_binaries} please make sure to set \code{binpath}
#' accordingly.
#'
#' The C++ version of OSLOM generates temporary folders and/or files that are stored in the \code{path_temp} folder
#' (folder "oslom_temp" in the working directory by default). This temporary folder is removed by default
#' (\code{delete_temp = TRUE}).
#'
#' @return A \code{data.frame} providing one or several modules (according to the chosen option(s)) for each node. If
#' \code{reassign = simil} or \code{reassign = random} only one column by hierarchical level will be provided. If
#' \code{reassign = no} the number of columns will depend on the number of overlapping modules (up to three). A value of 0
#' indicates that no module were found (i.e. non-overlapping nodes).
#'
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}),
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) and
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' @seealso \link{netclu_infomap}, \link{netclu_louvain}
#' @examples
#' comat <- matrix(sample(1000, 50), 5, 10)
#' rownames(comat) <- paste0("Site", 1:5)
#' colnames(comat) <- paste0("Species", 1:10)
#'
#' net <- similarity(comat, metric = "Simpson")
#' # com=netclu_oslom(net) # run install_binaries() to use this function
#' @references
#' \insertRef{Lancichinetti2011}{bioRgeo}
#' @export
netclu_oslom <- function(net, weight = TRUE, reassign = "no", r = 10, hr = 50,
                         seed = 1, t = 0.1, cp = 0.5, directed = FALSE,
                         bipartite = FALSE, primary_col = 1, feature_col = 2, remove_feature = TRUE,
                         delete_temp = TRUE, path_temp = "oslom_temp", binpath = NULL) {

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

  # Check if OSLOM has successfully been installed
  if (!directed) {
    if (!file.exists(paste0(binpath, "/bin/OSLOM/check.txt"))) {
      stop("OSLOM is not installed... Please have a look at https//biorgeo.github.io/bioRgeo/articles/bin.html for more details.")
    }
  } else {
    if (!file.exists(paste0(binpath, "/bin/OSLOM/check.txt"))) {
      stop("OSLOM is not installed... Please have a look at https//biorgeo.github.io/bioRgeo/articles/bin.html for more details.")
    } else {
      if (!file.exists(paste0(binpath, "/bin/OSLOM/checkdir"))) {
        stop("The directed version of OSLOM is not installed... Please have a look at https//biorgeo.github.io/bioRgeo/articles/bin.html for more details")
      }
    }
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

  if (!(reassign %in% c("no", "random", "simil"))) {
    stop("The reassign method is not available.
     Please chose among the followings:
         no, random, simil")
  }

  if (reassign == "simil" & !weight) {
    stop("A reassignement based on similarity should not be use when weight equal FALSE")
  }

  if (!is.numeric(r)) {
    stop("r must be numeric")
  } else {
    if (r <= 0) {
      stop("r must be strictly higher than 0")
    }
    if ((r - floor(r)) > 0) {
      stop("r must be an integer strictly higher than 0")
    }
  }

  if (!is.numeric(hr)) {
    stop("hr must be numeric")
  } else {
    if (hr < 0) {
      stop("hr must be positive")
    }
    if ((hr - floor(hr)) > 0) {
      stop("hr must be an integer higher or equal to 0")
    }
  }

  if (!is.numeric(seed)) {
    stop("seed must be numeric")
  } else {
    if (seed <= 0) {
      stop("seed must be strictly higher than 0")
    }
    if ((seed - floor(seed)) > 0) {
      stop("seed must be an integer higher or equal to 0")
    }
  }

  if (!is.numeric(t)) {
    stop("t must be numeric")
  } else {
    if (t < 0 | t > 1) {
      stop("t must be in the interval (0,1)")
    }
  }

  if (!is.numeric(cp)) {
    stop("cp must be numeric")
  } else {
    if (cp < 0 | cp > 1) {
      stop("cp must be in the interval (0,1)")
    }
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

  # Prepare input for OSLOM
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
  }

  # Export input in OSLOM folder in WD
  utils::write.table(netemp, paste0(path_temp, "/net.txt"), row.names = FALSE, col.names = FALSE, sep = " ")

  # Prepare command to run OSLOM
  cmd <- paste0("-r ", r, " -hr ", hr, " -seed ", seed, " -t ", t, " -cp ", cp)
  if (weight) {
    cmd <- paste0("-f ", path_temp, "/net.txt -w ", cmd)
  } else {
    cmd <- paste0("-f ", path_temp, "/net.txt -uw ", cmd)
  }

  # Run OSLOM
  if (os == "Linux") {
    if (directed) {
      cmd <- paste0(binpath, "/bin/OSLOM/oslom_dir_lin ", cmd, " > /dev/null 2>&1")
    } else {
      cmd <- paste0(binpath, "/bin/OSLOM/oslom_undir_lin ", cmd, " > /dev/null 2>&1")
    }
    system(command = cmd)
  } else if (os == "Windows") {
    if (directed) {
      cmd <- paste0(binpath, "/bin/OSLOM/oslom_dir_win.exe ", cmd)
    } else {
      cmd <- paste0(binpath, "/bin/OSLOM/oslom_undir_win.exe ", cmd)
    }
    dir.create(paste0(path_temp, "/net.txt_oslo_files"), showWarnings = FALSE, recursive = TRUE)
    system(command = cmd, show.output.on.console = FALSE)
  } else if (os == "Darwin") {
    if (directed) {
      cmd <- paste0(binpath, "/bin/OSLOM/oslom_dir_mac ", cmd, " > /dev/null 2>&1")
    } else {
      cmd <- paste0(binpath, "/bin/OSLOM/oslom_undir_mac ", cmd, " > /dev/null 2>&1")
    }
    system(command = cmd)
  } else {
    stop("Linux, Windows or Mac distributions only.")
  }

  # Control: if the command line did not work
  if (!("tp" %in% list.files(paste0(path_temp, "/net.txt_oslo_files")))) {
    stop("Command line was wrongly implemented. OSLOM did not run.")
  }

  # Number of levels
  nblev <- 1

  # Retrieve output from tp [TO COMMENT]
  com <- readLines(paste0(path_temp, "/net.txt_oslo_files/tp"))
  cl <- list()
  length(cl) <- length(com) / 2
  for (k in 1:length(com)) {
    if ((k / 2 - trunc(k / 2)) == 0) {
      cl[[(k / 2)]] <- as.numeric(as.matrix(strsplit(com[k], split = " ")[[1]]))
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
        cl[[(k / 2)]] <- as.numeric(as.matrix(strsplit(com[k], split = " ")[[1]]))
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
        cl[[(k / 2)]] <- as.numeric(as.matrix(strsplit(com[k], split = " ")[[1]]))
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
    com <- com[, -c((dim(comtp)[2] + 1), (dim(comtp)[2] + dim(comtph)[2] + 1))]
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
        cl[[(k / 2)]] <- as.numeric(as.matrix(strsplit(com[k], split = " ")[[1]]))
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
    com <- com[, -c((dim(comtp)[2] + 1), (dim(comtp)[2] + dim(comtph)[2] + 1), (dim(comtp)[2] + dim(comtph)[2] + dim(comtphh) + 1))]
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

  # Remove feature nodes
  if (bipartite & remove_feature) {
    com <- com[match(idprim, com[, 1]), ]
  }

  # Return output
  com[, 1] <- as.character(com[, 1])
  return(com)
}
