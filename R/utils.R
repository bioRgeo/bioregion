# Controls #####################################################################
controls <- function(args = NULL, data = NULL, type = "input_net") {
  
  lstype <- c("input_bioregionalization",
              "input_nhandhclu",
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
  
  # Input bioregionalization ###################################################
  if (type == "input_bioregionalization") {
    if (!inherits(data, "bioregion.clusters")) {
      stop(paste0(deparse(substitute(data)), 
                  " must be a bioregion.clusters object."),
           call. = FALSE)
    }else{
      if(is.null(data$inputs$data_type)){
        stop(paste0(deparse(substitute(data)),
                    " is a bioregion.cluster object but it has been altered ",
                    "and some informations regarding the name of the algorithm ",
                    " data type and node type are missing."),
              call. = FALSE)
      }
      if(!is.na(data$inputs$data_type)){
        if(!(data$inputs$data_type %in% c("occurrence","abundance"))){
          stop(paste0(deparse(substitute(data)),
                      " is a bioregion.cluster object but it has been altered ",
                      "and some informations regarding the name of the algorithm ",
                      " data type and node type are missing."),
               call. = FALSE)
        }
      }
      if(is.null(data$inputs$node_type)){
        stop(paste0(deparse(substitute(data)),
                    " is a bioregion.cluster object but it has been altered ",
                    "and some informations regarding the name of the algorithm ",
                    " data type and node type are missing."),
             call. = FALSE)
      }
      if(!(data$inputs$node_type %in% c("site","species","both"))){
        stop(paste0(deparse(substitute(data)),
                    " is a bioregion.cluster object but it has been altered ",
                    "and some informations regarding the name of the algorithm ",
                    " data type and node type are missing."),
             call. = FALSE)
      }

      
    }
  }
  
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

# Species-to-bioregions/bioregionalization metrics
# Site-to-chorotype/chorological metrics 
sbgc <- function(clusters, 
                 bioregion_metrics,
                 bioregionalization_metrics,
                 comat,
                 type,  # sb or gc
                 data){ # occurrence, abundance or both
  
  # Initialization output
  res1 <- NULL
  res12 <- NULL
  res11 <- NULL
  res2 <- NULL
  res21 <- NULL
  res22 <- NULL
  
  # sb
  col1 <- "Species"
  col2 <- "Bioregion"
  colcoren <- c("n_sb", "n_s", "n_b")
  colcorew <- c("w_sb", "w_s", "w_b")
  
  # gc
  if(type == "gc"){
    comat <- t(comat)  
    col1 <- "Site"
    col2 <- "Chorotypes"
    colcoren <- c("n_gc", "n_g", "n_c")
    colcorew <- c("w_gc", "w_g", "w_c")
  }
  
  # Occurrence
  if((data != "abundance") |
     (data == "abundance" & "Rho" %in% bioregion_metrics) |
     (data == "abundance" & "NSpecificity" %in% bioregion_metrics) |
     (data == "abundance" & "Indval" %in% bioregion_metrics) |
     (data == "abundance" & "NIndval" %in% bioregion_metrics)){
    
    # comat_bin
    comat_bin <- comat
    comat_bin[comat_bin > 0] <- 1
  
    # CoreTerms
    temp <- aggregate(comat_bin, list(clusters), sum)
    nij_mat <- t(as.matrix(temp[,-1]))
    rownames(nij_mat) <- colnames(temp)[-1]
    colnames(nij_mat) <- temp[,1]

    ni_mat <- matrix(apply(comat_bin, 2, sum), 
                     dim(comat_bin)[2], 
                     dim(nij_mat)[2])
    rownames(ni_mat) <- rownames(nij_mat)
    colnames(ni_mat) <- colnames(nij_mat)
    
    temp <- aggregate(comat_bin, list(clusters), length)
    nj_mat <- t(as.matrix(temp[,-1]))
    rownames(nj_mat) <- rownames(nij_mat)
    colnames(nj_mat) <- colnames(nij_mat)
    
    n <- sum(nj_mat[1,])
    
    # Normalized for NSpecificity & NIndVal
    if("NSpecificity" %in% bioregion_metrics |
       "NIndval" %in% bioregion_metrics){
      Nnij_mat <- nij_mat / nj_mat
      Nnij_mat <- Nnij_mat / apply(Nnij_mat, 1, sum)
      Nnij_mat[is.na(Nnij_mat)] <- 0
    }

    
    # Output bioregions
    if(!is.null(bioregion_metrics) & data != "abundance"){
      
      res11 <- cbind(mat_to_net(nij_mat, weight = TRUE, remove_zeroes = FALSE),
                     mat_to_net(ni_mat, weight = TRUE, remove_zeroes = FALSE)[,3],
                     mat_to_net(nj_mat, weight = TRUE, remove_zeroes = FALSE)[,3])
      colnames(res11) <- c(col1, col2, colcoren) 
      
      nij <- res11[,3]
      ni <- res11[,4]
      nj <- res11[,5]

      # Specificity 
      if("Specificity" %in% bioregion_metrics){
        res11$Specificity_occ <- nij / ni
      }
      
      # NSpecificity 
      if("NSpecificity" %in% bioregion_metrics){
        
        tempnspe <- mat_to_net(Nnij_mat, weight = TRUE, remove_zeroes = FALSE)
        
        res11$NSpecificity_occ <- tempnspe[,3]
      }
      
      # Fidelity 
      if("Fidelity" %in% bioregion_metrics){
        res11$Fidelity_occ <- nij / nj
      }
      
      # IndVal 
      if("IndVal" %in% bioregion_metrics){
        res11$IndVal_occ <- (nij / ni) * (nij / nj)
      }
      
      # NIndVal 
      if("NIndVal" %in% bioregion_metrics){
        
        tempniv <- mat_to_net(Nnij_mat, weight = TRUE, remove_zeroes = FALSE)
        
        res11$NIndVal_occ <- tempniv[,3] * (nij / nj)
      }
      
      # Rho
      if("Rho" %in% bioregion_metrics){
        
        num <- nij-((ni*nj)/n)
        den <- sqrt((nj*(n-nj)/(n-1))*(ni/n)*(1-(ni/n)))
        # den[is.na(den)] <- ??
        
        res11$Rho_occ <- num/den
      }
      
      # CoreTerms
      if(!("CoreTerms" %in% bioregion_metrics)){
        res11 <- res11[,-c(3,4,5)]
      }
    }
    
    # Output bioregionalizations
    if(!is.null(bioregionalization_metrics) & data != "abundance"){
      
      res21 <- data.frame(nij_mat[,1],nij_mat[,1])
      res21[,1] <- rownames(nij_mat)
      colnames(res21) <- c(col1, "Dummy")
      
      if("P" %in% bioregionalization_metrics){
        res21$P_occ <- 1 - apply((nij_mat / ni_mat)*(nij_mat / ni_mat), 1 , sum)
      }
      
      res21 <- res21[,-2]
      rownames(res21) <- 1:dim(res21)[1]
      
    }
  }
  
  # Abundance
  if(data != "occurrence"){
    
    # CoreTerms
    temp <- aggregate(comat, list(clusters), sum)
    wij_mat <- t(as.matrix(temp[,-1]))
    rownames(wij_mat) <- colnames(temp)[-1]
    colnames(wij_mat) <- temp[,1]
    
    wi_mat <- matrix(apply(comat, 2, sum), 
                     dim(comat)[2], 
                     dim(wij_mat)[2])
    rownames(wi_mat) <- rownames(wij_mat)
    colnames(wi_mat) <- colnames(wij_mat)
    
    w2i_mat <- matrix(apply(comat*comat, 2, sum), 
                     dim(comat)[2], 
                     dim(wij_mat)[2])
    rownames(w2i_mat) <- rownames(wij_mat)
    colnames(w2i_mat) <- colnames(wij_mat)
    
    wj_mat <- matrix(apply(wij_mat,2,sum), 
                     dim(comat)[2], 
                     dim(wij_mat)[2],
                     byrow=TRUE)
    rownames(wj_mat) <- rownames(wij_mat)
    colnames(wj_mat) <- colnames(wij_mat)
    
    # Normalized for NSpecificity & NIndVal
    if("NSpecificity" %in% bioregion_metrics |
       "NIndval" %in% bioregion_metrics){
      Nwij_mat <- wij_mat / nj_mat
      Nwij_mat <- Nwij_mat / apply(Nwij_mat, 1, sum)
      Nwij_mat[is.na(Nwij_mat)] <- 0
    }
    
    if("Rho" %in% bioregion_metrics){
      muij_mat <- wij_mat / nj_mat
      muij_mat[is.na(muij_mat)] <- 0
      
      mui_mat <- wi_mat / n
      vari_mat <- (1/(n-1)) * (w2i_mat-wi_mat*wi_mat/n)
    }
    
    # Output bioregions
    if(!is.null(bioregion_metrics)){
      
      res12 <- cbind(mat_to_net(wij_mat, weight = TRUE, remove_zeroes = FALSE),
                     mat_to_net(wi_mat, weight = TRUE, remove_zeroes = FALSE)[,3],
                     mat_to_net(wj_mat, weight = TRUE, remove_zeroes = FALSE)[,3])
      colnames(res12) <- c(col1, col2, colcorew)
      
      wij <- res12[,3]
      wi <- res12[,4]
      wj <- res12[,5]
      
      # Specificity
      if("Specificity" %in% bioregion_metrics){
        res12$Specificity_abund <- wij / wi
      }
      
      # NSpecificity
      if("NSpecificity" %in% bioregion_metrics){
        tempnspe <- mat_to_net(Nwij_mat, weight = TRUE, remove_zeroes = FALSE)
        
        res12$NSpecificity_abund <- tempnspe[,3]
      }
      
      # Fidelity 
      if("Fidelity" %in% bioregion_metrics){
        res12$Fidelity_abund <- wij / wj
      }
      
      # IndVal 
      if("IndVal" %in% bioregion_metrics){
        tempindval <- cbind(mat_to_net(nij_mat, weight = TRUE, 
                                       remove_zeroes = FALSE),
                            mat_to_net(nj_mat, weight = TRUE, 
                                       remove_zeroes = FALSE)[,3])
        
        nij <- tempindval[,3]
        nj <- tempindval[,4]
        
        res12$IndVal_abund <- (wij / wi) * (nij / nj)
      }
      
      # NIndVal 
      if("NIndVal" %in% bioregion_metrics){
        tempnindval <- cbind(mat_to_net(Nwij_mat, weight = TRUE, 
                                        remove_zeroes = FALSE),
                             mat_to_net(nij_mat, weight = TRUE, 
                                       remove_zeroes = FALSE)[,3],
                             mat_to_net(nj_mat, weight = TRUE, 
                                       remove_zeroes = FALSE)[,3])
        
        Nwij <- tempnindval[,3]
        nij <- tempnindval[,4]
        nj <- tempnindval[,5]
        
        res12$NIndVal_abund <- Nwij * (nij / nj)
      }
      
      # Rho
      if("Rho" %in% bioregion_metrics){
        temprho <- cbind(mat_to_net(muij_mat, weight = TRUE, remove_zeroes = FALSE),
                         mat_to_net(mui_mat, weight = TRUE, remove_zeroes = FALSE)[,3],
                         mat_to_net(vari_mat, weight = TRUE, remove_zeroes = FALSE)[,3],
                         mat_to_net(nj_mat, weight = TRUE, remove_zeroes = FALSE)[,3])
        
        muij <- temprho[,3]
        mui <- temprho[,4]
        vari <- temprho[,5]
        nj <- temprho[,6]
        
        num <- muij-mui
        den <- sqrt((n-nj)/(n-1)*(vari/nj))
        # den[is.na(den)] <- ??
        
        res12$Rho_abund <- num/den
      }
      
      # CoreTerms
      if(!("CoreTerms" %in% bioregion_metrics)){
        res12 <- res12[,-c(3,4,5)]
      }
    }
    
    # Output bioregionalizations
    if(!is.null(bioregionalization_metrics)){
      
      res22 <- data.frame(wij_mat[,1],wij_mat[,1])
      res22[,1] <- rownames(wij_mat)
      colnames(res22) <- c(col1, "Dummy")
      
      if("P" %in% bioregionalization_metrics){
        res22$P_abund <- 1 - apply((wij_mat / wi_mat)*(wij_mat / wi_mat),1,sum)
      }
      
      res22 <- res22[,-2]
      rownames(res22) <- 1:dim(res22)[1]
      
    }
  }  
  
  # Combine outputs
  if(!is.null(res11) & !is.null(res12)){
    res1 <- cbind(res11, res12[, -c(1,2), drop = FALSE])
  }else{
    if(!is.null(res11)){
      res1 <- res11
    }
    if(!is.null(res12)){
      res1 <- res12
    }
  }
  if(!is.null(res21) & !is.null(res22)){
    res2 <- cbind(res21, res22[, -1, drop = FALSE])
  }else{
    if(!is.null(res21)){
      res2 <- res21
    }
    if(!is.null(res22)){
      res2 <- res22
    }
  }
  
  # Return output
  res <- list()
  res$bioregion1 <- res1
  res$bioregion2 <- res2
  
  return(res)
  
}

# Site-to-bioregions/bioregionalization metrics
gb <- function(clusters,
               bioregion_metrics,
               bioregionalization_metrics,
               comat,
               similarity,
               #data,  # occurrence, abundance or both
               include_cluster){ 
  
  # Initialization output
  res1 <- NULL
  res2 <- NULL
  
  # Check needed inputs
  comat_needed <- (("Richness" %in% bioregion_metrics) |
                   ("Rich_Endemics" %in% bioregion_metrics) |
                   ("Prop_Endemics" %in% bioregion_metrics))
  
  sim_needed <- (("MeanSim" %in% bioregion_metrics) |
                 ("SdSim" %in% bioregion_metrics) |
                 ("Silhouette" %in% bioregionalization_metrics))
  
  # Precompute muij if sim_needed
  if(sim_needed){
    
    diag(similarity) <- NA
    
    temp <- aggregate(similarity, list(clusters), mean, na.rm=TRUE)
    muij_mat <- t(as.matrix(temp[,-1]))
    rownames(muij_mat) <- colnames(temp)[-1]
    colnames(muij_mat) <- temp[,1]
    muij_mat[is.na(muij_mat)] <- 0
    
  }
  
  # Output bioregions
  if(!is.null(bioregion_metrics)){
    
    # Create base matrix site x site (diag=1, 0 otherwise)
    if(comat_needed){
      comat_bin <- comat
      comat_bin[comat_bin > 0] <- 1
      base <- comat_bin %*% t(comat_bin)
      base[!diag(dim(base)[1])] <- 0
    }else{
      base <- similarity
      base[!diag(dim(base)[1])] <- 0
      diag(base) <- 1
    }
     
    # Initialize res1
    temp <- aggregate(base, list(clusters), sum)
    res1 <- t(as.matrix(temp[,-1]))
    res1[res1>0] <- 1
    rownames(res1) <- colnames(temp)[-1]
    colnames(res1) <- temp[,1]
    
    res1 <- mat_to_net(res1, weight=TRUE, remove_zeroes = FALSE)
    colnames(res1) <- c("Site", "Bioregion", "Assigned")
    
    # Richness
    if("Richness" %in% bioregion_metrics |
       "Prop_Endemics" %in% bioregion_metrics){
      
      temp <- comat_bin %*% t(comat_bin)
      temp[!diag(dim(temp)[1])] <- 0
      temp <- aggregate(temp, list(clusters), sum)
      
      ng_mat <- t(as.matrix(temp[,-1]))
      ng_mat[] <- apply(ng_mat, 1, sum)
      
      rownames(ng_mat) <- colnames(temp)[-1]
      colnames(ng_mat) <- temp[,1]

      res1$Richness <- mat_to_net(ng_mat, 
                                  weight = TRUE, 
                                  remove_zeroes = FALSE)[,3]
    }
      
    # Rich_Endemics
    if("Rich_Endemics" %in% bioregion_metrics|
       "Prop_Endemics" %in% bioregion_metrics){
      
      # Species x cluster (1 if species in cluster)
      temp <- aggregate(comat_bin, list(clusters), max)
      is_sb <- t(as.matrix(temp[, -1]))
      rownames(is_sb) <- colnames(temp)[-1]
      colnames(is_sb) <- temp[,1]
      
      # Set 0 for none endemic
      is_sb[apply(is_sb, 1, sum) > 1] = 0
      
      # Rich_Endemics
      nge_mat <- comat_bin %*% is_sb

      res1$Rich_Endemics <- mat_to_net(nge_mat, 
                                       weight = TRUE, 
                                       remove_zeroes = FALSE)[,3]

    }
    
    # Prop_Endemics
    if("Prop_Endemics" %in% bioregion_metrics){
      
      res1$Prop_Endemics <- res1$Rich_Endemics / res1$Richness
      
      if(!("Richness" %in% bioregion_metrics)){
        res1$Richness <- NULL
      }
      if(!("Rich_Endemics" %in% bioregion_metrics)){
        res1$Rich_Endemics <- NULL
      }
      
    }  
    
    if("MeanSim" %in% bioregion_metrics){
      res1$MeanSim <- mat_to_net(muij_mat, 
                                 weight = TRUE, 
                                 remove_zeroes = FALSE)[,3]
    }
    if("SdSim" %in% bioregion_metrics){
      
      temp <- aggregate(similarity, list(clusters), sd, na.rm=TRUE)
      sdij_mat <- t(as.matrix(temp[,-1]))
      rownames(sdij_mat) <- colnames(temp)[-1]
      colnames(sdij_mat) <- temp[,1]
      sdij_mat[is.na(sdij_mat)] <- 0
      
      res1$SdSim <- mat_to_net(sdij_mat, 
                               weight = TRUE, 
                               remove_zeroes = FALSE)[,3]
    }

    # include_cluster
    if(!include_cluster){
      res1 <- res1[,-3]
    }
    
  }
  
  # Output bioregionalizations
  if(!is.null(bioregionalization_metrics)){
    
    res2 <- data.frame(muij_mat[,1], muij_mat, clusters)
    res2[,1] <- rownames(muij_mat)
    colnames(res2) <- c("Site", colnames(muij_mat), "Assigned")
    
    if(dim(res2)[2] == 3){ # Only one site
      nob <- TRUE
      res2$a <- res2[,2]
      res2$b <- NA
    }else{
      nob <- FALSE
      res2$a <- apply(res2, 1, function(x) {
        # x[2:(ncol-1)] = MeanSim
        meansim_values <- as.numeric(x[2:(ncol(res2)-1)])
        # Assigned bioregion
        assigned <- x[ncol(res2)]
        # POTENTIAL PROBLEM WITH NUMERIC WHEN > 10 [5 become " 5"]

        # Extract meansim corresponding to the assigned bioregion
        a_val <- meansim_values[which(colnames(res2)[2:(ncol(res2)-1)] == assigned)]
        return(a_val)
      })
      res2$b <- apply(res2, 1, function(x) {
        meansim_values <- as.numeric(x[2:(ncol(res2)-2)])
        assigned <- x[ncol(res2)-1]  # colonne assigned
        # Put NA for the Assigned
        meansim_values[colnames(res2)[2:(ncol(res2)-2)] == assigned] <- NA
        # b = max among other bioregions
        b_val <- max(meansim_values, na.rm = TRUE)
        return(b_val)
      })
    }
    
    res2 <- res2[, c(1, (dim(res2)[2]-1), dim(res2)[2])]
    
    if("Silhouette" %in% bioregionalization_metrics){
      if(nob){
        res2$Silhouette <- NA
      }else{
        res2$Silhouette <- (res2$a - res2$b) / pmax(res2$a,res2$b)
      }
    }

    res2 <- res2[,-c(2,3)]
    rownames(res2) <- 1:dim(res2)[1]
    
  }
  
  # Return output
  res <- list()
  res$bioregion1 <- res1
  res$bioregion2 <- res2
  
  return(res)
  
}



################################################################################

#' Detect data type from metric name
#' 
#' Determines whether a similarity/dissimilarity metric is based on
#' occurrence (presence/absence) or abundance (quantitative) data.
#' 
#' @param metric Character string with metric name
#' 
#' @return Character: "occurrence", "abundance", or "unknown"
#' 
#' @details
#' Occurrence metrics (using abc): Jaccard, Jaccardturn, Sorensen, Simpson, abc
#' Betapart occurrence metrics (case insensitive): beta.sim, beta.sne,
#' beta.sor, beta.jtu, beta.jne, beta.jac
#' Abundance metrics (using ABC): Bray, Brayturn, ABC
#' Betapart abundance metrics (case insensitive): beta.bray.bal, beta.bray.gra,
#' beta.bray, beta.ruz.bal, beta.ruz.gra, beta.ruz
#' Unknown: Euclidean, custom formulas, or NA
#' 
#' @noRd
detect_data_type_from_metric <- function(metric) {
  if (is.na(metric) || is.null(metric)) {
    return(NA)
  }

  occurrence_metrics <- c("abc", "Jaccard", "Jaccardturn", "Sorensen", "Simpson")
  abundance_metrics <- c("ABC", "Bray", "Brayturn")

  betapart_occurrence <- c("beta.sim", "beta.sne", "beta.sor",
                           "beta.jtu", "beta.jne", "beta.jac")
  betapart_abundance <- c("beta.bray.bal", "beta.bray.gra", "beta.bray",
                          "beta.ruz.bal", "beta.ruz.gra", "beta.ruz")

  metric_lower <- tolower(metric)

  if (metric %in% occurrence_metrics) {
    return("occurrence")
  } else if (metric %in% abundance_metrics) {
    return("abundance")
  } else if (metric_lower %in% betapart_occurrence) {
    return("occurrence")
  } else if (metric_lower %in% betapart_abundance) {
    return("abundance")
  } else if (metric == "Euclidean") {
    return(NA)
  } else {
    # Custom formula or unknown metric
    return(NA)
  }
}
