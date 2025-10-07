# Controls #####################################################################
controls <- function(args = NULL, data = NULL, type = "input_net") {
  
  lstype <- c("input_nhandhclu",
              "input_similarity",
              "input_dissimilarity",
              "input_pairwise",
              "input_conversion_similarity",
              "input_conversion_dissimilarity",
              "input_net",
              "input_net_directed",
              "input_net_isdirected",
              "input_net_isloop",
              "input_net_weight",
              "input_net_index",
              "input_net_index_value",
              "input_net_index_positive_value",
              "input_net_bip",
              "input_net_bip_col", 
              "input_matrix",
              "input_dist",
              "input_data_frame_nhandhclu",
              "input_data_frame",
              "character",
              "character_vector",
              "boolean",
              "boolean_vector",
              "numeric",
              "numeric_vector",
              "positive_numeric",
              "positive_numeric_vector",
              "strict_positive_numeric",
              "strict_positive_numeric_vector",
              "integer",
              "integer_vector",
              "positive_integer",
              "positive_integer_vector",
              "strict_positive_integer",
              "strict_positive_integer_vector",
              "character_or_positive_integer")
  
  if(!(type %in% lstype)){
    stop("Control type not defined!", call.=FALSE)
  }
  
  # TODO: reformat all error messages to single lines, using the following
  # format:
  # paste0("This is a multiline ",
  #        "error sentence ",
  #        "with no problematic line ",
  #        "breaks")
  
  # Input nhandhclu ############################################################
  if (type == "input_nhandhclu") {
    if (!inherits(data, "bioregion.pairwise") &
        !inherits(data, "dist") &
        !is.data.frame(data)) {
      stop(paste0(deparse(substitute(data)), 
                  " is not a bioregion.pairwise object, ", 
                  "a dissimilarity matrix (class dist) or ",
                  "a data.frame with at least 3 columns ", 
                  "(site1, site2 and your dissimilarity index)."),
           call. = FALSE)
    }
  }
  
  # Input similarity ###########################################################
  if (type == "input_similarity") {
    if(!inherits(data, "bioregion.pairwise")) {
      # message(paste0(deparse(substitute(data)),
      #                " is not a bioregion.pairwise object.\n", 
      #                "Note that some functions required dissimilarity metrics ", 
      #                "(hclu_ & nhclu_) and others similarity metrics ",
      #                "(netclu_). Please carefully check your data before ", 
      #                "using the clustering functions."))
    }else{
      if(is.null(attr(data, "type"))){
        message(paste0(deparse(substitute(data)),
                       " is a bioregion.pairwise object but it has not ",
                       "been possible to identify the object's type ",
                       "(similarity or dissimilarity) probably because the ",
                       "bioregion.pairwise object has been altered.\n",
                       "Note that some functions required dissimilarity ",
                       "metrics (hclu_ & nhclu_) and others similarity ",
                       "metrics (netclu_). Please carefully check your data ",
                       "before using the clustering functions."))
      }else{
        if (attr(data, "type") == "dissimilarity") {
          stop(paste0(deparse(substitute(data)),
                      " seems to be a dissimilarity object. This function ",
                      "should be applied on similarities, not ",
                      "dissimilarities. Use dissimilarity_to_similarity() ",
                      "before using this function."), call. = FALSE)
        }
      }
    }
  }
  
  # Input dissimilarity ########################################################
  if (type == "input_dissimilarity") {
    if(!inherits(data, "bioregion.pairwise")) {
      # message(paste0(deparse(substitute(data)),
      #                " is not a bioregion.pairwise object.\n", 
      #                "Note that some functions required dissimilarity metrics ", 
      #                "(hclu_ & nhclu_) and others similarity metrics ",
      #                "(netclu_). Please carefully check your data before ", 
      #                "using the clustering functions."))
    }else{
      if(is.null(attr(data, "type"))){
        message(paste0(deparse(substitute(data)),
                       " is a bioregion.pairwise object but it has not ",
                       "been possible to identify the object's type ",
                       "(similarity or dissimilarity) probably because the ",
                       "bioregion.pairwise object has been altered.\n",
                       "Note that some functions required dissimilarity ",
                       "metrics (hclu_ & nhclu_) and others similarity ",
                       "metrics (netclu_). Please carefully check your data ",
                       "before using the clustering functions."))
      }else{
        if (attr(data, "type") == "similarity") {
          stop(paste0(deparse(substitute(data)),
                      " seems to be a similarity object. This function ",
                      "should be applied on dissimilarities, not ",
                      "similarities. Use similarity_to_dissimilarity() ",
                      "before using this function."), call. = FALSE)
        }
      }
    }
  }
  
  # Input pairwise #############################################################
  if (type == "input_pairwise") {
    
    if (!inherits(data, "bioregion.pairwise")) {
      stop(paste0(deparse(substitute(data)), 
                  " should be a bioregion.pairwise object created by ",
                  "similarity() or dissimilarity_to_similarity()."),
           call. = FALSE)
    }
    if(is.null(attr(data, "type"))){
      stop(paste0(deparse(substitute(data)),
                  " is a bioregion.pairwise object but it has not ",
                  "been possible to identify the object's type (similarity or ",
                  " dissimilarity) probably because the ",
                  "bioregion.pairwise object has been altered."),
           call. = FALSE)
    }
  }
  
  # Input conversion similarity ################################################
  if (type == "input_conversion_similarity") {
    
    if (!inherits(data, "bioregion.pairwise")) {
      stop(paste0(deparse(substitute(data)), 
                  " should be a bioregion.pairwise object created by ",
                  "similarity() or dissimilarity_to_similarity()."),
                  call. = FALSE)
    }
    if(is.null(attr(data, "type"))){
      stop(paste0(deparse(substitute(data)),
                  " is a bioregion.pairwise object but it has not ",
                  "been possible to identify the object's type (similarity or ",
                  " dissimilarity) probably because the ",
                  "bioregion.pairwise object has been altered."),
           call. = FALSE)
    }
    if (attr(data, "type") == "dissimilarity") {
      stop(paste0(deparse(substitute(data)), 
                  " is already composed of dissimilarity metrics. If you want ",
                  "to convert it to similarity, use ", 
                  "dissimilarity_to_similarity()."),
           call. = FALSE)
    }

  }
  
  # Input conversion dissimilarity #############################################
  if (type == "input_conversion_dissimilarity") {
    
    if (!inherits(data, "bioregion.pairwise")) {
      stop(paste0(deparse(substitute(data)), 
                  " should be a bioregion.pairwise object created by ",
                  "dissimilarity() or similarity_to_dissimilarity()."),
           call. = FALSE)
    }
    if(is.null(attr(data, "type"))){
      stop(paste0(deparse(substitute(data)),
                  " is a bioregion.pairwise object but it has not ",
                  "been possible to identify the object's type (similarity or ",
                  "dissimilarity) probably because the ",
                  "bioregion.pairwise object has been altered."),
           call. = FALSE)
    }
    if (attr(data, "type") == "similarity") {
      stop(paste0(deparse(substitute(data)), 
                  " is already composed of similarity metrics. If you want to ",
                  "convert it to dissimilarity, use ", 
                  "similarity_to_dissimilarity()."),
           call. = FALSE)
    }
  }  
  
  # Input network ##############################################################
  if (type == "input_net") {
    if (!is.data.frame(data)) {
      stop(paste0(deparse(substitute(data)), 
                  " must be a data.frame."),
           call. = FALSE)
    }
    if (dim(data)[2] < 2) {
      stop(paste0(deparse(substitute(data)),
                  " must be a data.frame with at least two columns."),
           call. = FALSE)
    }
    pairs1 <- paste0(data[, 1], "_", data[, 2])
    if (sum(duplicated(pairs1)) > 0) {
      stop(paste0("The first two columns of ", 
                  deparse(substitute(data)),
                  " contain duplicated pairs of nodes!"), 
           call. = FALSE)
    }
    nbna <- sum(is.na(data[,1:2]))
    if (nbna > 0) {
      stop(paste0("NA(s) detected in ", 
                  deparse(substitute(data)),"."), 
           call. = FALSE)
    }
  }
  
  # Input network loop or directed #############################################
  if (type == "input_net_directed") {
    if (length(args) > 1) {
      stop(paste0(deparse(substitute(args)), 
                  " must be of length 1."),
           call. = FALSE
      )
    }
    if (!is.logical(args)) {
      stop(paste0(deparse(substitute(args)), 
                  " must be a boolean."),
           call. = FALSE)
    }
    #data = data[data[,1] != data[,2],]
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
    #data = data[data[,1] != data[,2],]
    pairs1 <- paste0(data[, 1], "_", data[, 2])
    pairs2 <- paste0(data[, 2], "_", data[, 1])
    if (length(intersect(pairs1, pairs2)) > 0) {
       stop(paste0("The network is directed, ",
                   "this function is designed for undirected networks!"),
            call. = FALSE)
    }
  }
  if (type == "input_net_isloop") {
    if (sum(data[,1] == data[,2]) > 0) {
      stop("The network contains self-loop(s)!",
           call.=FALSE)
    }
  }
  
  # Input network weight #######################################################
  if (type == "input_net_weight") {
    if (length(args) > 1) {
      stop(paste0(deparse(substitute(args)), 
                  " must be of length 1."),
           call. = FALSE
      )
    }
    if (!is.logical(args)) {
      stop(paste0(deparse(substitute(args)), 
                  " must be a boolean."),
           call. = FALSE)
    }
    if (args & dim(data)[2] == 2) {
      stop(paste0(deparse(substitute(data)),
                  " must be a data.frame with at least three columns ",
                  "if weight equal TRUE."), 
           call. = FALSE)
    }

  }
  
  # Input network index ########################################################
  if (type == "input_net_index") {
    if (length(args) > 1) {
      stop(paste0(deparse(substitute(args)), 
                  " must be of length 1."),
           call. = FALSE
      )
    }
    if (is.character(args)) {
      if (!(args %in% colnames(data)[-(1:2)])) {
        stop(paste0("If ", deparse(substitute(args)),
                    " is a character, it should be a ",
                    "column name (and not the ",
                    "first or second column)."), 
             call. = FALSE)
      }
    } else if (is.numeric(args)) {
      if (args %% 1 != 0) {
        stop(paste0("If ", 
                    deparse(substitute(args)),
                    " is numeric, it should be an integer."), 
             call. = FALSE)
      } else {
        if ((args <= 2)) {
          stop(paste0(deparse(substitute(args)),
                      " should be strictly higher than 2."), 
               call. = FALSE)
        }
        if ((args > dim(data)[2])) {
          stop(paste0(deparse(substitute(args)),
                      " should be lower or equal to ", 
                      dim(data)[2], "."),
               call. = FALSE)
        }
      }
    } else {
      stop(paste0(deparse(substitute(args)),
                  " should be numeric or character."), 
           call. = FALSE)
    }
  }
  
  # Input network index value ##################################################
  if (type == "input_net_index_value") {
    nbna <- sum(is.na(data[,3]))
    if (nbna > 0) {
      stop("NA(s) detected in the weight column.", 
           call. = FALSE)
    }
    if (!is.numeric(data[, 3])) {
      stop("The weight column must be numeric.", 
           call. = FALSE)
    } 
  }
  
  # Input network index positive value #########################################
  if (type == "input_net_index_positive_value") {
    nbna <- sum(is.na(data[,3]))
    if (nbna > 0) {
      stop("NA(s) detected in the weight column.", call. = FALSE)
    }
    if (!is.numeric(data[, 3])) {
      stop("The weight column must be numeric.", call. = FALSE)
    } else {
      minet <- min(data[, 3])
      if (minet < 0) {
        stop(paste0("The weight column should ",
                    "contain only positive values."), 
             call. = FALSE)
      }
    }
  }
  
  # Input network bip #########################################################
  if (type == "input_net_bip") {
    if (length(intersect(data[, 1], data[, 2])) > 0) {
      stop("The network is not bipartite!", 
           call. = FALSE)
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
    } else if (is.numeric(args)) {
      if ((args != 1) & (args != 2)) {
        stop(paste0("If ", deparse(substitute(args)),
                    " is numeric, it should be equal to 1 or 2."),
             call. = FALSE)
      }
    } else {
      stop(paste0(deparse(substitute(args)),
                  " should be numeric or character."), 
           call. = FALSE)
    }
  }

  # Input matrix ###############################################################
  if (type == "input_matrix") {
    if (!is.matrix(data)) {
      stop(paste0(deparse(substitute(data)), " must be a matrix."),
           call. = FALSE)
    }
    rowmat <- rownames(data)
    colmat <- colnames(data)
    if (sum(duplicated(rowmat)) > 0) {
      stop("Duplicated rownames detected!", 
           call. = FALSE)
    }
    if (sum(duplicated(colmat)) > 0) {
      stop("Duplicated colnames detected!", 
           call. = FALSE)
    }
    nbna <- sum(is.na(data))
    if (nbna > 0) {
      stop("NA(s) detected in the matrix!", 
           call. = FALSE)
    }
  }
  
  # Input dist ###############################################################
  if (type == "input_dist") {
    if (!inherits(data, "dist")) {
      stop(paste0(deparse(substitute(data)), " must be a dist object."),
           call. = FALSE)
    }
    if (!is.numeric(data)) {
      stop(paste0(deparse(substitute(data)), " must be numeric."),
           call. = FALSE)
    }
    nbna <- sum(is.na(data))
    if (nbna > 0) {
      stop(paste0("NA(s) detected in ", 
                  paste0(deparse(substitute(data))), "."),
           call. = FALSE)
    }
  }

  # Input data.frame nhandhclu #################################################
  if (type == "input_data_frame_nhandhclu") {
    if (!is.data.frame(data)) {
      stop(paste0(deparse(substitute(data)), 
                  " must be a data.frame."),
           call. = FALSE)
    }
    if (dim(data)[2] < 3) {
      stop(paste0(deparse(substitute(data)),
                  " must be a data.frame with at least three columns."),
           call. = FALSE)
    }
    nbna <- sum(is.na(data))
    if (nbna > 0) {
      stop(paste0("NA(s) detected in ", 
                  deparse(substitute(data)),"."), 
           call. = FALSE)
    }

    if (sum(data[,1] == data[,2]) > 0) {
      stop(paste0(
        deparse(substitute(data)),
        " contains rows with the same site on both columns!"
      ), call. = FALSE)
    }
        
    pairs1 <- paste0(data[, 1], "_", data[, 2])
    if (sum(duplicated(pairs1)) > 0) {
      stop(paste0(
        "The first two columns of ", 
        deparse(substitute(data)),
        " contain duplicated pairs of sites!"
      ), call. = FALSE)
    }
    pairs2 <- paste0(data[, 2], "_", data[, 1])
    if (length(intersect(pairs1, pairs2)) > 0) {
      stop(paste0(
        "The first two columns of ", 
        deparse(substitute(data)),
        " contain (unordered) duplicated pairs of sites!"
      ), call. = FALSE)
    }
  }
  
  # Input data.frame ###########################################################
  if (type == "input_data_frame") {
    if (!is.data.frame(data)) {
      stop(paste0(deparse(substitute(data)), " must be a data.frame."),
           call. = FALSE)
    }
    #rowmat <- rownames(data)
    #colmat <- colnames(data)
    #if (sum(duplicated(rowmat)) > 0) {
    #  message("Duplicated rownames detected!")
    #}
    #if (sum(duplicated(colmat)) > 0) {
    #  message("Duplicated colnames detected!")
    #}
    nbna <- sum(is.na(data))
    if (nbna > 0) {
      stop("NA(s) detected in the data.frame!", call. = FALSE)
    }
  }
  
  # Character ##################################################################
  if (type == "character") {
    if (length(args) > 1) {
      stop(paste0(deparse(substitute(args)), " must be of length 1."),
           call. = FALSE
      )
    }
    if (!is.character(args)) {
      stop(paste0(deparse(substitute(args)), " must be a character."),
           call. = FALSE
      )
    }
  }
  
  # Character vector ###########################################################
  if (type == "character_vector") {
    if (!is.character(args)) {
      stop(paste0(deparse(substitute(args)), " must be a character."),
           call. = FALSE
      )
    }
  }
  
  # Boolean ####################################################################
  if (type == "boolean") {
    if (length(args) > 1) {
      stop(paste0(deparse(substitute(args)), " must be of length 1."),
           call. = FALSE
      )
    }
    if (!is.logical(args) || is.na(args)) {
      stop(paste0(deparse(substitute(args)), " must be a boolean."),
           call. = FALSE
      )
    }
  }
  
  # Boolean vector #############################################################
  if (type == "boolean_vector") {
    if (!is.logical(args) || is.na(args)) {
      stop(paste0(deparse(substitute(args)), " must be a boolean."),
           call. = FALSE
      )
    }
  }
  
  # Numeric ###################################################################
  if (type == "numeric") {
    if (length(args) > 1) {
      stop(paste0(deparse(substitute(args)), " must be of length 1."),
           call. = FALSE
      )
    }
    if (!is.numeric(args)) {
      stop(paste0(deparse(substitute(args)), " must be numeric."),
           call. = FALSE
      )
    }
  }
  
  # Numeric vector #############################################################
  if (type == "numeric_vector") {
    if (!is.numeric(args)) {
      stop(paste0(deparse(substitute(args)), " must be numeric."),
           call. = FALSE
      )
    }
  }
  
  # Positive numeric ##########################################################
  if (type == "positive_numeric") {
    if (length(args) > 1) {
      stop(paste0(deparse(substitute(args)), " must be of length 1."),
           call. = FALSE
      )
    }
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
  
  # Positive numeric vector ####################################################
  if (type == "positive_numeric_vector") {
    if (!is.numeric(args)) {
      stop(paste0(deparse(substitute(args)), " must be numeric."),
           call. = FALSE)
    } else {
      if (sum(args < 0) > 0) {
        stop(paste0(deparse(substitute(args)), 
                    " must be composed of values higher than 0."),
             call. = FALSE
        )
      }
    }
  }
  
  # Strict positive numeric ###################################################
  if (type == "strict_positive_numeric") {
    if (length(args) > 1) {
      stop(paste0(deparse(substitute(args)), " must be of length 1."),
           call. = FALSE
      )
    }
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
  
  # Strict positive numeric vector #############################################
  if (type == "strict_positive_numeric_vector") {
    if (!is.numeric(args)) {
      stop(paste0(deparse(substitute(args)), " must be numeric."),
           call. = FALSE)
    } else {
      if (sum(args <= 0) > 0) {
        stop(paste0(deparse(substitute(args)),
                    " must be composed of values strictly higher than 0."),
             call. = FALSE)
      }
    }
  }
  
  # Integer ###################################################################
  if (type == "integer") {
    if (length(args) > 1) {
      stop(paste0(deparse(substitute(args)), " must be of length 1."),
           call. = FALSE
      )
    }
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
  
  # Integer vector #############################################################
  if (type == "integer_vector") {
    if (!is.numeric(args)) {
      stop(paste0(deparse(substitute(args)), " must be numeric."),
           call. = FALSE)
    } else {
      if (sum(args %% 1 != 0) > 0) {
        stop(paste0(deparse(substitute(args)), " must be composed of integers."),
             call. = FALSE
        )
      }
    }
  }
  
  # Positive integer ##########################################################
  if (type == "positive_integer") {
    if (length(args) > 1) {
      stop(paste0(deparse(substitute(args)), " must be of length 1."),
           call. = FALSE
      )
    }
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
  
  # Positive integer vector ####################################################
  if (type == "positive_integer_vector") {
    if (!is.numeric(args)) {
      stop(paste0(deparse(substitute(args)), " must be numeric."),
           call. = FALSE)
    } else {
      if (sum(args %% 1 != 0) > 0) {
        stop(paste0(deparse(substitute(args)), 
                    " must be composed of integers."),
             call. = FALSE
        )
      } else {
        if (sum(args < 0) > 0) {
          stop(paste0(deparse(substitute(args)), 
                      " must be composed of values higher than 0."),
               call. = FALSE
          )
        }
      }
    }
  }
  
  # Strict positive integer ###################################################
  if (type == "strict_positive_integer") {
    if (length(args) > 1) {
      stop(paste0(deparse(substitute(args)), " must be of length 1."),
           call. = FALSE
      )
    }
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
  
  # Strict positive integer vector #############################################
  if (type == "strict_positive_integer_vector") {
    if (!is.numeric(args)) {
      stop(paste0(deparse(substitute(args)), " must be numeric."),
           call. = FALSE)
    } else {
      if (sum(args %% 1 != 0) > 0) {
        stop(paste0(deparse(substitute(args)), 
                    " must be composed of integers."),
             call. = FALSE
        )
      } else {
        if (sum(args <= 0) > 0) {
          stop(paste0(deparse(substitute(args)),
                      " must be composed of values strictly higher than 0."), 
               call. = FALSE)
        }
      }
    }
  }
  
  # Character or positive integer ##############################################
  if (type == "character_or_positive_integer") {
    if (length(args) > 1) {
      stop(paste0(deparse(substitute(args)), " must be of length 1."),
           call. = FALSE
      )
    }
    if (!is.character(args) && !is.numeric(args)) {
      stop(paste0(deparse(substitute(args)), 
                  " must be a character string or a positive integer."),
           call. = FALSE
      )
    }
    if (is.numeric(args)) {
      if (args %% 1 != 0) {
        stop(paste0(deparse(substitute(args)), " must be an integer."),
             call. = FALSE
        )
      }
      if (args <= 0) {
        stop(paste0(deparse(substitute(args)),
                    " must be strictly higher than 0."), 
             call. = FALSE)
      }
    }
  }
}

# Additional functions #########################################################
# reformat_hierarchy
reformat_hierarchy <- function(input, algo = "infomap", integerize = FALSE) {
  
  # Infomap
  if(algo == "infomap"){
    
    sep <- ":"
    
    # Input
    input <- as.character(as.vector(as.matrix(input)))
    
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
    
    for (k in 3:(nblev + 1)) {
      table[table[, k] == 0, (k - 1)] <- 0
    }
    table[, (nblev + 1)] <- 0
    
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
  
  }
  
  # Algo
  if (algo == "louvain") {
    
    id0 <- which(input[, 1] == 0)
    output <- input[(id0[1] + 1):(id0[2] - 1), ]
    colnames(output) <- c("ID", paste0("V",1))
    
    if(length(id0)>2){
      for(k in 2:(length(id0)-1)){
        output$temp <- NA
        tabk <- input[(id0[k] + 1):(id0[k+1] - 1),]
        output[,(k+1)] <- tabk[match(output[,k], tabk[,1]),2]
        colnames(output)[(k+1)] <- paste0("V",k)
        
      }
    }
    
  }
  
  return(output)
}

# knbclu
knbclu <- function(partitions, 
                   reorder = TRUE, 
                   rename_duplicates = TRUE) {
  
  # Identify the number of clusters per partition
  nb <- dim(partitions)[2] - 1
  nbclus <- apply(partitions[, 2:(nb + 1), drop = FALSE],
                  2,
                  function(x) length(unique(x[!is.na(x)])))
  
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

# make.unique.2
# from https://stackoverflow.com/questions/7659891/r-make-unique-starting-in-1
make.unique.2 <- function(x, sep = ".") {
  stats::ave(x, x, FUN = function(a) {
    if (length(a) > 1) {
      paste(a, 1:length(a), sep = sep)
    } else {
      a
    }
  })
}

# seedrng
seedrng <- function() {
  # as.numeric(as.POSIXct(Sys.time())) + sample(-10:10, 1)
  #sample(1:(.Machine$integer.max), 1)
  sample(1:10000, 1)
}



