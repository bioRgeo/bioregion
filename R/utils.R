controls <- function(args = NULL, data = NULL, type = "input_net") {
  
  # Input similarity ##########################################################
  if (type == "input_similarity") {
    if (inherits(data, "bioRgeo.pairwise.metric")) {
      if (attr(data, "type") == "dissimilarity") {
        stop(paste0(deparse(substitute(data)),
                    " seems to be a dissimilarity object. 
This function should be applied on similarities, not dissimilarities. 
Use dissimilarity_to_similarity() before using this function."), call. = FALSE)
      }
    } else {
      message(paste0(deparse(substitute(data)),
                     " is not a bioRgeo.pairwise.metric object. 
Note that some functions required dissimilarity metrics (hclu_ & nhclu) and
others similarity metrics (netclu_). 
Please carefully check your data before using the clustering functions."),
              call. = FALSE)
    }
  }
  
  # Input dissimilarity #######################################################
  if (type == "input_dissimilarity") {
    if (inherits(data, "bioRgeo.pairwise.metric")) {
      if (attr(data, "type") == "similarity") {
        stop(paste0(deparse(substitute(data)),
                    " seems to be a similarity object.
This function should be applied on dissimilarities, not similarities.
Use similarity_to_dissimilarity() before using this function."),
             call. = FALSE
        )
      }
    } else {
      message(paste0(deparse(substitute(data)),
                     " is not a bioRgeo.pairwise.metric object. 
Note that some functions required dissimilarity metrics (hclu_ & nhclu) and
others similarity metrics (netclu_). 
Please carefully check your data before using the clustering functions."))
    }
  }
  
  # Input network #############################################################
  if (type == "input_net") {
    if (!is.data.frame(data)) {
      stop(paste0(deparse(substitute(data)), " must be a data.frame."),
           call. = FALSE)
    }
    if (dim(data)[2] < 2) {
      stop(paste0(deparse(substitute(data)),
                  " must be a data.frame with at least two columns."),
           call. = FALSE)
    }
    pairs1 <- paste0(data[, 1], "_", data[, 2])
    if (sum(duplicated(pairs1)) > 0) {
      stop(paste0(
        "The first two columns of ", deparse(substitute(data)),
        " contain duplicated pairs of nodes!"
      ), call. = FALSE)
    }
    nbna <- sum(is.na(data))
    if (nbna > 0) {
      stop("NA(s) detected in the data.frame!", call. = FALSE)
    }
  }
  
  # Input network directed ####################################################
  if (type == "input_net_directed") {
    if (!is.logical(args)) {
      stop(paste0(deparse(substitute(args)), " must be a boolean."),
           call. = FALSE)
    }
    pairs1 <- paste0(data[, 1], "_", data[, 2])
    pairs2 <- paste0(data[, 2], "_", data[, 1])
    if (!args) {
      if (length(intersect(pairs1, pairs2)) > 0) {
        stop(paste0(deparse(substitute(data)),
                    " should not be directed if directed = FALSE."),
             call. = FALSE)
      }
    }
  }
  if (type == "input_net_isdirected") {
    pairs1 <- paste0(data[, 1], "_", data[, 2])
    pairs2 <- paste0(data[, 2], "_", data[, 1])
    if (length(intersect(pairs1, pairs2)) > 0) {
      message(paste0("It seems that the network is directed!
                        This function is designed for undirected networks!"))
    }
  }
  
  # Input network weight ######################################################
  if (type == "input_net_weight") {
    if (!is.logical(args)) {
      stop(paste0(deparse(substitute(args)), " must be a boolean."),
           call. = FALSE)
    }
    if (args & dim(data)[2] == 2) {
      stop(paste0(
        deparse(substitute(args)),
        " must be a data.frame with at least three columns if weight equal 
        TRUE."), call. = FALSE)
    }
  }
  
  # Input network index #######################################################
  if (type == "input_net_index") {
    if (is.character(args)) {
      if (!(args %in% colnames(data)[-(1:2)])) {
        stop(paste0(deparse(substitute(args)),
                    " is a character, it should be a column name (and not the
                    first or second column)."), call. = FALSE)
      }
    } else if (is.factor(args)) {
      args <- as.character(args)
      if (!(args %in% colnames(data)[-(1:2)])) {
        stop(paste0("If ", deparse(substitute(args)),
                    " is a character, it should be a column name (and not the
                    first or second column)."),
             call. = FALSE
        )
      }
    } else if (is.numeric(args)) {
      if (args %% 1 != 0) {
        stop(paste0("If ", deparse(substitute(args)),
                    " is numeric, it should be an integer."), call. = FALSE)
      } else {
        if ((args <= 2)) {
          stop(paste0(deparse(substitute(args)),
                      " should be stricltly higher than 2."), call. = FALSE)
        }
        if ((args > dim(data)[2])) {
          stop(paste0(deparse(substitute(args)),
                      " should be lower or equal to ", dim(data)[2], "."),
               call. = FALSE)
        }
      }
    } else {
      stop(paste0(deparse(substitute(args)),
                  " should be numeric or character."), call. = FALSE)
    }
  }
  
  # Input network index value #################################################
  if (type == "input_net_index_value") {
    if (!is.numeric(data[, 3])) {
      stop("The (dis)similarity metric must be numeric.", call. = FALSE)
    } else {
      minet <- min(data[, 3])
      if (minet < 0) {
        stop(
          "The (dis)similarity metric should contain only positive reals:
          negative value(s) detected!", call. = FALSE)
      }
    }
  }
  
  # Input network bip #########################################################
  if (type == "input_net_bip") {
    if (length(intersect(data[, 1], data[, 2])) > 0) {
      stop("The network is not bipartite!", call. = FALSE)
    }
  }
  
  # Input network bip_col #####################################################
  if (type == "input_net_bip_col") {
    if (is.character(args)) {
      if (!(args %in% colnames(data)[1:2])) {
        stop(paste0(
          "If ", deparse(substitute(args)),
          " is a character, it should be the first or second column name."),
          call. = FALSE)
      }
    } else if (is.factor(args)) {
      args <- as.character(args)
      if (!(args %in% colnames(data)[1:2])) {
        stop(paste0(
          "If ", deparse(substitute(args)),
          " is a character, it should be the first or second column name."),
          call. = FALSE)
      }
    } else if (is.numeric(args)) {
      if ((args != 1) & (args != 2)) {
        stop(paste0("If ", deparse(substitute(args)),
                    " is numeric, it should be equal to 1 or 2."),
             call. = FALSE)
      }
    } else {
      stop(paste0(deparse(substitute(args)),
                  " should be numeric or character."), call. = FALSE)
    }
  }
  
  # Input matrix ##############################################################
  if (type == "input_matrix") {
    if (!is.matrix(data)) {
      stop(paste0(deparse(substitute(data)), " must be a matrix."),
           call. = FALSE)
    }
    rowmat <- rownames(data)
    colmat <- colnames(data)
    if (sum(duplicated(rowmat)) > 0) {
      message("Duplicated rownames detected!")
    }
    if (sum(duplicated(colmat)) > 0) {
      message("Duplicated colnames detected!")
    }
    nbna <- sum(is.na(data))
    if (nbna > 0) {
      stop("NA(s) detected in the matrix!", call. = FALSE)
    }
  }
  
  # Character #################################################################
  if (type == "character") {
    if (!is.character(args)) {
      stop(paste0(deparse(substitute(args)), " must be a character."),
           call. = FALSE
      )
    }
    if (is.factor(args)) {
      args <- as.character(args)
    }
    return(args)
  }
  
  # Boolean ###################################################################
  if (type == "boolean") {
    if (!is.logical(args)) {
      stop(paste0(deparse(substitute(args)), " must be a boolean"),
           call. = FALSE
      )
    }
  }
  
  # Numeric ###################################################################
  if (type == "numeric") {
    if (!is.numeric(args)) {
      stop(paste0(deparse(substitute(args)), " must be numeric."),
           call. = FALSE
      )
    }
  }
  
  # Positive numeric ##########################################################
  if (type == "positive_numeric") {
    if (!is.numeric(args)) {
      stop(paste0(deparse(substitute(args)), " must be numeric."),
           call. = FALSE)
    } else {
      if (args < 0) {
        stop(paste0(deparse(substitute(args)), " must be higher than 0."),
             call. = FALSE
        )
      }
    }
  }
  
  # Strict positive numeric ###################################################
  if (type == "strict_positive_numeric") {
    if (!is.numeric(args)) {
      stop(paste0(deparse(substitute(args)), " must be numeric."),
           call. = FALSE)
    } else {
      if (args <= 0) {
        stop(paste0(deparse(substitute(args)),
                    " must be strictly higher than 0."), call. = FALSE)
      }
    }
  }
  
  # Integer ###################################################################
  if (type == "integer") {
    if (!is.numeric(args)) {
      stop(paste0(deparse(substitute(args)), " must be numeric."),
           call. = FALSE)
    } else {
      if (args %% 1 != 0) {
        stop(paste0(deparse(substitute(args)), " must be an integer."),
             call. = FALSE
        )
      }
    }
  }
  
  # Positive integer ##########################################################
  if (type == "positive_integer") {
    if (!is.numeric(args)) {
      stop(paste0(deparse(substitute(args)), " must be numeric."),
           call. = FALSE)
    } else {
      if (args %% 1 != 0) {
        stop(paste0(deparse(substitute(args)), " must be an integer."),
             call. = FALSE
        )
      } else {
        if (args < 0) {
          stop(paste0(deparse(substitute(args)), " must be higher than 0."),
               call. = FALSE
          )
        }
      }
    }
  }
  
  # Strict positive integer ###################################################
  if (type == "strict_positive_integer") {
    if (!is.numeric(args)) {
      stop(paste0(deparse(substitute(args)), " must be numeric."),
           call. = FALSE)
    } else {
      if (args %% 1 != 0) {
        stop(paste0(deparse(substitute(args)), " must be an integer."),
             call. = FALSE
        )
      } else {
        if (args <= 0) {
          stop(paste0(deparse(substitute(args)),
                      " must be strictly higher than 0."), call. = FALSE)
        }
      }
    }
  }
}

reformat_hierarchy <- function(input, algo = "infomap", integerize = FALSE) {
  
  # Input
  input <- as.character(as.vector(as.matrix(input)))
  
  # Algo
  if (algo == "infomap") {
    sep <- ":"
  }
  if (algo == "optics") {
    sep <- "."
  }
  
  # Nb levels
  nblev <- max(lengths(regmatches(input, gregexpr(paste0("\\", sep),
                                                  input)))) + 1
  
  # From input to table
  table <- tidyr::separate(
    data = data.frame(input),
    col = input,
    remove = FALSE,
    into = paste0("lvl", 1:nblev),
    sep = paste0("\\", sep),
    fill = "right"
  )
  
  # Replace NA by O
  table[which(is.na(table), arr.ind = TRUE)] <- 0
  
  # Replace last number by 0 for infomap
  if (algo == "infomap") {
    for (k in 3:(nblev + 1)) {
      table[table[, k] == 0, (k - 1)] <- 0
    }
    table[, (nblev + 1)] <- 0
  }
  
  # Output
  output <- table
  for (lvl in grep("lvl", colnames(output))[2:nblev]) {
    output[, lvl] <- paste(output[, lvl - 1], output[, lvl], sep = ".")
  }
  output[grep("lvl", colnames(output))] <- lapply(
    output[grep("lvl", colnames(output))],
    function(x) gsub("\\.0", "", x)
  )
  # Integerize
  if (integerize) {
    for (k in 2:(nblev + 1)) {
      output[, k] <- as.numeric(as.factor(output[, k]))
    }
  }
  
  return(output)
}

knbclu <- function(partitions, method = "length",
                   reorder = TRUE, rename_duplicates = TRUE) {
  
  # Identify the number of clusters per partition
  nb <- dim(partitions)[2] - 1
  
  
  if (method == "max") {
    nbclus <- as.numeric(apply(
      partitions[, 2:(nb + 1), drop = FALSE],
      2,
      function(x) max(x)
    ))
  } else if (method == "length") {
    nbclus <- apply(
      partitions[, 2:(nb + 1), drop = FALSE],
      2,
      function(x) length(unique(x))
    )
  }
  
  # Rename and reorder
  if (reorder) {
    ord <- cbind(2:(nb + 1), nbclus)
    ord <- ord[order(ord[, 2]), , drop = FALSE]
    partitions <- partitions[, c(1, ord[, 1])]
    colnames(partitions)[2:(nb + 1)] <- paste0("K_", ord[, 2])
  } else {
    colnames(partitions)[2:(nb + 1)] <- paste0("K_", nbclus)
  }
  
  # Rename duplicates
  if (rename_duplicates) {
    colnames(partitions)[2:(nb + 1)] <- make.unique.2(
      colnames(partitions)[2:(nb + 1)],
      sep = "_"
    )
  }
  
  # Convert in character
  for (k in 1:(nb + 1)) {
    partitions[, k] <- as.character(partitions[, k])
  }
  
  # Change colnames 1 en ID
  colnames(partitions)[1] <- "ID"
  
  partitions
}

# From https://stackoverflow.com/questions/7659891/r-make-unique-starting-in-1
make.unique.2 <- function(x, sep = ".") {
  stats::ave(x, x, FUN = function(a) {
    if (length(a) > 1) {
      paste(a, 1:length(a), sep = sep)
    } else {
      a
    }
  })
}
