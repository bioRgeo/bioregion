#' Convert a matrix or list of matrices to a bioregion (dis)similarity object
#'
#' Converts a (dis)similarity `matrix` or a `list` of such matrices into a 
#' `bioregion.pairwise.metric` object compatible with the `bioregion` package. 
#' The input can come from base R, `dist` objects, or outputs from other 
#' packages.
#'
#' @param mat A `matrix`, a `dist` object, or a `list` of these representing 
#'   pairwise similarity or dissimilarity values to convert into a 
#'   `bioregion.pairwise.metric` object. This function can also directly handle 
#'   outputs from other R packages (see the `pkg` argument).
#'   
#' @param metric_name Optional `character` vector or single `character` string 
#'   specifying the name of the (dis)similarity metric(s), which will appear as 
#'   column names in the output (see Note).
#'
#' @param pkg An optional `character` string indicating the name of the package 
#'   from which `mat` was generated (`NULL` by default, see Details). 
#'   Available options are `"adespatial"`, `"betapart"`, `"ecodist"`, or 
#'   `"vegan"`.
#'
#' @param is_similarity A `logical` value indicating whether the input data 
#'   represents similarity (`TRUE`) or dissimilarity (`FALSE`).
#' 
#' @details
#' This function can directly handle outputs from nine functions across four 
#' packages:
#'
#' - **adespatial**: [beta.div][adespatial::beta.div], 
#'   [beta.div.comp][adespatial::beta.div.comp]
#' - **betapart**: [beta.pair][betapart::beta.pair], 
#'   [beta.pair.abund][betapart::beta.pair.abund], 
#'   [betapart.core][betapart::betapart.core]
#' - **ecodist**: [distance][ecodist::distance], 
#'   [bcdist][ecodist::bcdist]
#' - **vegan**: [vegdist][vegan::vegdist], 
#'   [designdist][vegan::designdist]
#'
#' See the documentation of these packages for more information:
#' - https://cran.r-project.org/package=adespatial
#' - https://cran.r-project.org/package=betapart
#' - https://cran.r-project.org/package=ecodist
#' - https://cran.r-project.org/package=vegan
#'
#' @note
#' If no specific package is specified (i.e., `pkg = NULL`), site names will be 
#' based on the row names of the first matrix. If row names are `NULL`, they 
#' will be generated automatically. If `mat` is a named list, those names will 
#' be used as column names only if `metric_name = NULL`.
#' 
#' @return 
#' A dissimilarity or similarity object of class `bioregion.pairwise.metric`, 
#' compatible with the `bioregion` package.
#' 
#' @seealso 
#' For more details illustrated with a practical example, 
#' see the vignette: 
#' \url{https://biorgeo.github.io/bioregion/articles/a3_pairwise_metrics.html}.
#' 
#' Associated functions: 
#' [dissimilarity] [similarity] 
#' 
#' @author
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) \cr
#' Boris Leroy (\email{leroy.boris@gmail.com}) \cr
#' Pierre Denelle (\email{pierre.denelle@gmail.com})
#' 
#' @examples
#' 
#' mat <- matrix(runif(100), 10, 10)
#' rownames(mat) <- paste0("s",1:10)
#' 
#' pair <- as_bioregion_pairwise(list(mat,mat,mat), 
#'                               metric_name = NULL,
#'                               pkg = NULL,
#'                               is_similarity = FALSE)
#'                               
#' pair
#' 
#' @export
as_bioregion_pairwise <- function(mat,
                                  metric_name = NULL,
                                  pkg = NULL,
                                  is_similarity = FALSE) {
  # Control metric_name
  if(!is.null(metric_name)){
    controls(args = metric_name, data = NULL, type = "character_vector")
    if(!is.null(pkg)){
      message("metric_name will be ignored when pkg is not NULL.")
      metric_name = NULL
    }
  }
  
  # Control is_similarity
  controls(args = is_similarity, data = NULL, type = "boolean")
  if(is_similarity & !is.null(pkg)){
    message("is_similarity will be ignored when pkg is not NULL.")
  }
  netype <- "dissimilarity"
  
  # Control pkg
  if(!is.null(pkg)){
    controls(args = pkg, data = NULL, type = "character")
    if (!(pkg %in% c("adespatial", "betapart", "ecodist", "vegan"))) {
      stop(paste0("Please choose pkg from the following:\n",
                  "adespatial, betapart, ecodist or vegan."), 
           call. = FALSE)
    }
  }  
  
  # pkg 
  if(is.null(pkg)){
    if(is_similarity){
      netype <- "similarity"
    }
  }else{
    
    # adespatial
    if(pkg == "adespatial"){
      
      betadiv <- c("beta","SCBD","LCBD","p.LCBD","p.adj","method","note","D")
      betadivcomp <- c("repl","rich","D","part","Note") 
      betadivcompabc <- c("repl","rich","D","part","Note","a","b","c") 
      
      if(!any(identical(names(mat),betadiv),
          identical(names(mat),betadivcomp),
          identical(names(mat),betadivcompabc))){
        stop("mat does not seem to be an output from the adespatial package.", 
             call. = FALSE)
      }else{
        if(identical(names(mat),betadiv)){
          if(is.na(mat$D[1])){
            stop("D is NULL. Check that save.D=TRUE.", 
                 call. = FALSE)
          }else{
            metric_name <- mat$method[1]
            mat <- mat$D
          }
        }
        if(identical(names(mat),betadivcomp) | 
           identical(names(mat),betadivcompabc)){
           if(identical(names(mat),betadivcompabc)){
             metric_name <- c(mat$Note, "a", "b", "c")
             mat <- list(mat$D, mat$a, mat$b, mat$c)
           }else{
             metric_name <- mat$Note
             mat <- mat$D
           }
        }
      }        
    }
    
    # betapart
    if(pkg == "betapart"){
      
      betapairj <- c("beta.jtu", "beta.jne", "beta.jac")
      betapairs <- c("beta.sim", "beta.sne", "beta.sor")
      betapairb <- c("beta.bray.bal", "beta.bray.gra", "beta.bray")
      betapairr <- c("beta.ruz.bal", "beta.ruz.gra", "beta.ruz")
      betacore <- c("data", "sumSi", "St", "a", "shared", "not.shared", 
                   "sum.not.shared", "max.not.shared", "min.not.shared")
      
      if(!any(identical(names(mat),betapairj),
              identical(names(mat),betapairs),
              identical(names(mat),betapairb),
              identical(names(mat),betapairr),
              identical(names(mat),betacore))){
        stop("mat does not seem to be an output from the betapart package.", 
             call. = FALSE)
      }else{
        if(identical(names(mat),betacore)){
          a <- mat$shared
          b <- mat$not.shared
          c <- t(b)
          mat <- list(a=a, b=b, c=c)
        }
      }
    }
    
    # ecodist
    if(pkg == "ecodist"){
      if(is.null(attr(mat, "method"))){
        stop("mat does not seem to be an output from the ecodist package.", 
             call. = FALSE)
      }else{
        metric_name <- attr(mat, "method")
      }
    }
    
    # vegan
    if(pkg == "vegan"){
      if(is.null(attr(mat, "method"))){
        stop("mat does not seem to be an output from the vegan package.", 
             call. = FALSE)
      }else{
        metric_name <- attr(mat, "method")
      }
    }
  }

  
  # Convert
  if(!(inherits(mat, c("matrix","dist","list")))){
    stop("mat must be a matrix, a dist object, or a list of these.", 
         call. = FALSE)
  }
  
  # if not transform mat in list of length 1 
  if (!inherits(mat, "list")) {
    lmat <- list()
    lmat[[1]] <- mat
    if(!is.null(metric_name) & length(metric_name)!=1){
      stop("metric_name should have the same length as mat.", 
           call. = FALSE)
    }
  }else{
    lmat <- mat
    if(!is.null(metric_name) & length(metric_name)!=length(lmat)){
      stop("metric_name should have the same length as mat.", 
           call. = FALSE)
    }
  }
  
  # Control loop over lmat
  nlmat <- NULL
  for(k in 1:length(lmat)){
    
    if(!(inherits(lmat[[k]], c("matrix","dist")))){
      stop("mat must be a matrix, a dist object, or a list of these.", 
           call. = FALSE)
    }  
    lmat[[k]] <- as.matrix(lmat[[k]])
    
    nlmatk <- dim(lmat[[k]])[1]
    nlmat <- c(nlmat,nlmatk)
    mlmatk <- dim(lmat[[k]])[2]
    if(nlmatk < 2 | (nlmatk != mlmatk)){
      stop(paste0("mat must be or contain only numeric (without NAs), ", 
                  "square (dis)similarity matrices between at least ",
                  "two sites."),
           call. = FALSE)
    }
    if (!is.numeric(lmat[[k]])) {
      stop(paste0("mat must be or contain only numeric (without NAs), ", 
                  "square (dis)similarity matrices between at least ",
                  "two sites."),
           call. = FALSE)
    }
    if (sum(is.na(lmat[[k]])) > 0) {
      stop(paste0("mat must be or contain only numeric (without NAs), ", 
                  "square (dis)similarity matrices between at least ",
                  "two sites."),
           call. = FALSE)
    }
  }
  if(length(unique(nlmat)) > 1){
    stop(paste0("mat must contain only square matrices with the same number ",
                "sites."),
         call. = FALSE)
  }
  
  # From mat_to_net
  if(is.null(rownames(lmat[[1]]))){
    rownames(lmat[[1]]) <- 1:dim(lmat[[1]])[1]
    colnames(lmat[[1]]) <- rownames(lmat[[1]])
  }else(
    colnames(lmat[[1]]) <- rownames(lmat[[1]])
  )
  
  net <- mat_to_net(lmat[[1]], 
                    weight = TRUE,
                    remove_zeroes = FALSE,
                    include_diag = FALSE,
                    include_lower = FALSE)
  
  if(length(lmat) > 1){
    for(k in 2:length(lmat)){
      netk <- mat_to_net(lmat[[k]], 
                         weight = TRUE,
                         remove_zeroes = FALSE,
                         include_diag = FALSE,
                         include_lower = FALSE)
      net <- cbind(net, netk[,3])
    }
  }
  
  # Rename
  colnames(net)[1:2] <- c("Site1","Site2")
  if(length(lmat)==1){
    if(!is.null(metric_name)){
      colnames(net)[3] <- metric_name
    }else{
      colnames(net)[3] <- "Metric"
    }
  }else{
    if(!is.null(metric_name)){
      colnames(net)[-c(1,2)] <- metric_name
    }else{
      if(!is.null(names(lmat))){
        colnames(net)[-c(1,2)] <- names(lmat)
      }else{
        colnames(net)[-c(1,2)] <- paste0("Metric",1:length(lmat))
      }
    }
  }
  
  attr(net, "type") <- netype
  attr(net, "nb_sites") <- dim(lmat[[1]])[1]
  attr(net, "nb_species") <- NA
  
  class(net) <- append("bioregion.pairwise.metric", class(net))
  
  return(net)
  
}
