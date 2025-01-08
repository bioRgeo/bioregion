#' Non-hierarchical clustering: Affinity Propagation
#'
#' This function performs non-hierarchical clustering using the Affinity 
#' Propagation algorithm.
#'
#' @param similarity The output object from [similarity()] or 
#' [dissimilarity_to_similarity()], or a `dist` object. If a `data.frame` is 
#' used, the first two columns should represent pairs of sites (or any pair of 
#' nodes), and the subsequent column(s) should contain the similarity indices.
#'
#' @param index The name or number of the similarity column to use. By default, 
#' the third column name of `similarity` is used.
#'
#' @param p Input preference, which can be a vector specifying individual 
#' preferences for each data point. If scalar, the same value is used for all 
#' data points. If `NA`, exemplar preferences are initialized based on the 
#' distribution of non-Inf values in the similarity matrix, controlled by `q`.
#'
#' @param q If `p = NA`, exemplar preferences are initialized according to the 
#' distribution of non-Inf values in the similarity matrix. By default, the 
#' median is used. A value between 0 and 1 specifies the sample quantile, 
#' where `q = 0.5` results in the median.
#'
#' @param maxits The maximum number of iterations to execute.
#'
#' @param convits The algorithm terminates if the exemplars do not change for 
#' `convits` iterations.
#'
#' @param lam The damping factor, a value in the range [0.5, 1). Higher values 
#' correspond to heavier damping, which may help prevent oscillations.
#'
#' @param details If `TRUE`, detailed information about the algorithm's progress 
#' is stored in the output object.
#'
#' @param nonoise If `TRUE`, disables the addition of a small amount of noise to 
#' the similarity object, which prevents degenerate cases.
#'
#' @param seed The seed for the random number generator.
#'
#' @param K The desired number of clusters. If not `NULL`, the function 
#' [apcluster][apcluster::apclusterK] is called.
#'
#' @param prc A parameter needed when `K` is not `NULL`. The algorithm stops if 
#' the number of clusters deviates by less than `prc` percent from the desired 
#' value `K`. Set to 0 to enforce exactly `K` clusters.
#'
#' @param bimaxit A parameter needed when `K` is not `NULL`. Specifies the 
#' maximum number of bisection steps to perform. No warning is issued if the 
#' number of clusters remains outside the desired range.
#'
#' @param exact A flag indicating whether to compute the initial preference 
#' range exactly.
#'
#' @param algorithm_in_output A `boolean` indicating whether to include the 
#' original output of [apcluster][apcluster::apcluster] in the result. Defaults 
#' to `TRUE`.
#'
#' @return
#' A `list` of class `bioregion.clusters` with five slots:
#' \enumerate{
#' \item{**name**: A `character` string containing the name of the algorithm.}
#' \item{**args**: A `list` of input arguments as provided by the user.}
#' \item{**inputs**: A `list` describing the characteristics of the clustering 
#' process.}
#' \item{**algorithm**: A `list` of objects associated with the clustering 
#' procedure, such as original cluster objects (if `algorithm_in_output = TRUE`).}
#' \item{**clusters**: A `data.frame` containing the clustering results.}}
#'
#' If `algorithm_in_output = TRUE`, the `algorithm` slot includes the output of 
#' [apcluster][apcluster::apcluster].
#'
#' @details
#' This function is based on the [apcluster](https://cran.r-project.org/package=apcluster) 
#' package ([apcluster][apcluster::apcluster]).
#' 
#' @references
#' Frey B & Dueck D (2007) Clustering by Passing Messages Between Data Points. 
#' \emph{Science} 315, 972-976.
#' 
#' @seealso 
#' For more details illustrated with a practical example, 
#' see the vignette: 
#' \url{https://biorgeo.github.io/bioregion/articles/a4_2_non_hierarchical_clustering.html}.
#' 
#' Associated functions: 
#' [nhclu_clara] [nhclu_clarans] [nhclu_dbscan] [nhclu_kmeans] [nhclu_affprop] 
#'
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) \cr
#' Boris Leroy (\email{leroy.boris@gmail.com}) \cr
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) 
#' 
#' @examples
#' comat_1 <- matrix(sample(0:1000, size = 10*12, replace = TRUE,
#' prob = 1/1:1001), 10, 12)
#' rownames(comat_1) <- paste0("Site", 1:10)
#' colnames(comat_1) <- paste0("Species", 1:12)
#' comat_1 <- cbind(comat_1,
#'                  matrix(0, 10, 8,
#'                         dimnames = list(paste0("Site", 1:10),
#'                                         paste0("Species", 13:20))))
#' 
#' comat_2 <- matrix(sample(0:1000, size = 10*12, replace = TRUE,
#'                          prob = 1/1:1001), 10, 12)
#' rownames(comat_2) <- paste0("Site", 11:20)
#' colnames(comat_2) <- paste0("Species", 9:20)
#' comat_2 <- cbind(matrix(0, 10, 8,
#'                         dimnames = list(paste0("Site", 11:20),
#'                                         paste0("Species", 1:8))),
#'                  comat_2)
#' 
#' comat <- rbind(comat_1, comat_2)
#'
#' dissim <- dissimilarity(comat, metric = "Simpson")
#' sim <- dissimilarity_to_similarity(dissim)
#' 
#' clust1 <- nhclu_affprop(sim)
#' 
#' clust2 <- nhclu_affprop(sim, q = 1)
#'
#' # Fixed number of clusters
#' clust3 <- nhclu_affprop(sim, K = 2, prc = 10, bimaxit = 20, exact = FALSE)
#' 
#' @importFrom apcluster apcluster apclusterK
#'         
#' @export
nhclu_affprop <- function(similarity, 
                          index = names(similarity)[3],
                          seed = NULL,
                          p = NA, 
                          q = NA, 
                          maxits = 1000, 
                          convits = 100,
                          lam = 0.9, 
                          # includeSim = FALSE,
                          details = FALSE, 
                          nonoise = FALSE, 
                          K = NULL, 
                          prc = NULL, 
                          bimaxit = NULL, 
                          exact = NULL,
                          algorithm_in_output = TRUE){
  
  # 1. Controls ---------------------------------------------------------------
  controls(args = NULL, data = similarity, type = "input_nhandhclu")
  if(!inherits(similarity, "dist")){
    controls(args = NULL, data = similarity, type = "input_similarity")
    controls(args = NULL, data = similarity, 
             type = "input_data_frame_nhandhclu")
    controls(args = index, data = similarity, type = "input_net_index")
    net <- similarity
    # Convert tibble into dataframe
    if(inherits(net, "tbl_df")){
      net <- as.data.frame(net)
    }
    net[, 3] <- net[, index]
    net <- net[, 1:3]
    controls(args = NULL, data = net, type = "input_net_index_value")
    dist.obj <- stats::as.dist(
      net_to_mat(net,
                 weight = TRUE, squared = TRUE, symmetrical = TRUE))
  } else {
    controls(args = NULL, data = similarity, type = "input_dist")
    dist.obj <- similarity
    if(is.null(labels(dist.obj))){
      attr(dist.obj, "Labels") <- paste0(1:attr(dist.obj, "Size"))
      message("No labels detected, they have been assigned automatically.")
    }
  }
  
  #controls(args = NULL, data = similarity,
  #        type = "input_conversion_similarity")
	
  if(!is.null(seed)){
    controls(args = seed, data = NULL, type = "strict_positive_integer")
  }
  
  if(length(p) == 1){
    if(!is.na(p)){
      if(!is.null(p)){
        controls(args = p, data = NULL, type = "numeric")
      }
    }  
  }else{
    controls(args = p, data = NULL, type = "numeric_vector")
    if(length(p) != length(unique(c(similarity[, 1], similarity[, 2])))){
      stop(paste0("If p is a vector, its length should be equal to the ",
                  "number of sites."), 
           call. = FALSE)
    }
  }
  
  if(length(q) > 1){
    stop("q must be of length 1.", call. = FALSE)
  }
  if(!is.na(q)){
    controls(args = q, data = NULL, type = "positive_numeric")
    if(q < 0 | q > 1){
      stop("q should be in the interval [0, 1].", call. = FALSE)
    }
  }
  
  controls(args = maxits, data = NULL, type = "positive_integer")
  controls(args = convits, data = NULL, type = "positive_integer")
  
  controls(args = lam, data = NULL, type = "positive_numeric")
  if(lam < 0.5 | lam >= 1){
    stop("lam should be in the interval [0.5, 1).", call. = FALSE)
  }
  
  # controls(args = includeSim, data = NULL, type = "boolean")
  controls(args = details, data = NULL, type = "boolean")
  controls(args = nonoise, data = NULL, type = "boolean")
  
  # Argument for desired number of clusters: positive integer, if not null
  # (default value) then we call apcluter::apclusterK() (with argument K)
  if(!is.null(K)){
    controls(args = K, data = NULL, type = "strict_positive_integer")
    if(is.null(prc)){
      stop(paste0("When K is not NULL, you need to define prc. ",
                  "prc is a percentage value."), call. = FALSE)
    }
    if(is.null(bimaxit)){
      stop("When K is not NULL, you need to define bimaxit.", call. = FALSE)
    }
    if(is.null(exact)){
      stop("When K is not NULL, you need to define exact.", call. = FALSE)
    }
  }
  
  if(!is.null(prc)){
    if(is.null(K)){
      message(paste0("prc argument will be considered only if K is not ",
                     "set to NULL."))
    }
    controls(args = prc, data = NULL, type = "positive_numeric")
    if (prc < 0 | prc > 100) {
      stop("prc should be in the interval [0, 100].")
    }
  }
  
  if(!is.null(bimaxit)){
    if(is.null(K)){
      message(paste0("bimaxit argument will be considered only if K is not ",
                     "set to NULL."))
    }
    controls(args = bimaxit, data = NULL, type = "strict_positive_integer")
  }
  
  if(!is.null(exact)){
    if(is.null(K)){
      message(paste0("exact argument will be considered only if K is not ",
                     "set to NULL."))
    }
    controls(args = exact, data = NULL, type = "boolean")
  }
  
  controls(args = algorithm_in_output, data = NULL, type = "boolean")
  
  sim <- NULL
  
  # 2. Function ---------------------------------------------------------------
  outputs <- list(name = "nhclu_affprop")
  
  outputs$args <- list(index = index,
                       seed = seed,
                       p = p,
                       q = q,
                       maxits = maxits,
                       convits = convits,
                       lam = lam,
                       # includeSim = includeSim,
                       details = details,
                       nonoise = nonoise,
                       K = K,
                       prc = prc,
                       bimaxit = bimaxit,
                       exact = exact,
                       algorithm_in_output = algorithm_in_output)
  
  outputs$inputs <- list(bipartite = FALSE,
                         weight = TRUE,
                         pairwise = TRUE,
                         pairwise_metric = ifelse(!inherits(similarity, 
                                                            "dist"), 
                                                  ifelse(is.numeric(index), 
                                                         names(net)[3], index), 
                                                  NA),
                         dissimilarity = FALSE,
                         nb_sites = attr(dist.obj, "Size"),
                         hierarchical = FALSE)
  
  outputs$algorithm <- list()
  
  outputs$clusters <- data.frame(matrix(ncol = 1,
                                        nrow = length(labels(dist.obj)),
                                        dimnames = list(labels(dist.obj),
                                                        "name")))
  
  outputs$clusters$name <- labels(dist.obj)
  
  # Square similarity matrix
  sim_square <- net_to_mat(similarity, 
                           weight = TRUE, 
                           squared = TRUE,
                           symmetrical = TRUE)
  
  ## 2.1. apclusterK ----------------------------------------------------------
  if(!is.null(K)){
    if(is.null(seed)){
      outputs$algorithm <- apcluster::apclusterK(s = sim_square,
                                                 K = K,
                                                 prc = prc,
                                                 bimaxit = bimaxit,
                                                 exact = FALSE,
                                                 maxits = maxits,
                                                 convits = convits,
                                                 lam = lam,
                                                 includeSim = FALSE,
                                                 details = details,
                                                 nonoise = nonoise,
                                                 seed = NA)
    }else{
      set.seed(seed)
      outputs$algorithm <- apcluster::apclusterK(s = sim_square,
                                                 K = K,
                                                 prc = prc,
                                                 bimaxit = bimaxit,
                                                 exact = FALSE,
                                                 maxits = maxits,
                                                 convits = convits,
                                                 lam = lam,
                                                 includeSim = FALSE,
                                                 details = details,
                                                 nonoise = nonoise)
      rm(.Random.seed, envir=globalenv())
    }
  }
  ## 2.2. apcluster -----------------------------------------------------------
  else{
    if(is.null(seed)){
      outputs$algorithm <- apcluster::apcluster(s = sim_square,
                                                p = p,
                                                q = q,
                                                maxits = maxits,
                                                convits = convits,
                                                lam = lam,
                                                includeSim = FALSE,
                                                details = details,
                                                nonoise = nonoise,
                                                seed = NA)
    }else{
      set.seed(seed)
      outputs$algorithm <- apcluster::apcluster(s = sim_square,
                                                p = p,
                                                q = q,
                                                maxits = maxits,
                                                convits = convits,
                                                lam = lam,
                                                includeSim = FALSE,
                                                details = details,
                                                nonoise = nonoise)
      rm(.Random.seed, envir=globalenv())
    }
  }
  
  # names(outputs$algorithm) <- paste0("K_", n_clust)
  
  # Convert output of apcluster into a data.frame with bioregions per site
  names(outputs$algorithm@clusters) <-
    paste0("K_", 1:length(outputs$algorithm@clusters))
  outputs$algorithm@clusters <- lapply(outputs$algorithm@clusters, names)
  
  outputs_df <- mapply(cbind, outputs$algorithm@clusters,
                       "K_" = names(outputs$algorithm@clusters),
                       SIMPLIFY = FALSE)
  outputs_df <- do.call(rbind, outputs_df)
  colnames(outputs_df) <- c("Site",
                            paste0("K_", length(outputs$algorithm@clusters)))
  
  outputs$clusters <- as.data.frame(outputs_df) #data.frame(outputs$clusters,
  # outputs_df)
  
  outputs$clusters <- knbclu(outputs$clusters, reorder = TRUE)
  
  outputs$cluster_info <- data.frame(
    partition_name = colnames(outputs$clusters)[2],
    n_clust = apply(outputs$clusters[, 2, drop = FALSE],
                    2, function(x) length(unique(x[!is.na(x)]))))
  
  class(outputs) <- append("bioregion.clusters", class(outputs))
  
  # Set algorithm in output
  if (!algorithm_in_output) {
    outputs$algorithm <- NA
  }
  
  return(outputs)
}
