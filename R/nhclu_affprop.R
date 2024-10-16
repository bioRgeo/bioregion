#' Non hierarchical clustering: Affinity Propagation
#'
#' This function performs non hierarchical
#' clustering on the Affinity Propagation algorithm.
#'
#' @param similarity the output object from [similarity()] or
#' [dissimilarity_to_similarity()], or a `dist` object. If a `data.frame` is
#' used, the first two columns represent pairs of sites (or any pair of nodes),
#' and the next column(s) are the dissimilarity indices.
#'
#' @param index name or number of the similarity column to use. By default, 
#' the third column name of `similarity` is used.
#' 
#' @param p input preference; can be a vector that specifies individual
#' preferences for each data point. If scalar, the same value is used for all
#' data points. If NA, exemplar preferences are initialized according to the
#' distribution of non-Inf values in the similarity matrix. How this is done is
#' controlled by the parameter q.
#' 
#' @param q if p=NA, exemplar preferences are initialized according to the
#' distribution of non-Inf values in the similarity matrix. If q=NA, exemplar
#' preferences are set to the median of non-Inf values in the similarity
#' matrix. If q is a value between 0 and 1, the sample quantile with
#' threshold q is used, whereas q=0.5 again results in the median.
#' 
#' @param maxits maximal number of iterations that should be executed
#' 
#' @param convits the algorithm terminates if the examplars have not changed
#' for convits iterations
#' 
#' @param lam damping factor; should be a value in the range [0.5, 1);
#' higher values correspond to heavy damping which may be needed if
#' oscillations occur
#' 
#' @param details if TRUE, more detailed information about the algorithm's
#' progress is stored in the output object
#' 
#' @param nonoise small amount of noise added to the similarity object to
#' prevent degenerate cases; disabled when set to TRUE.
#' 
#' @param seed seed of the random number generator.
#' 
#' @param K desired number of clusters. If not null, then the function
#' [apcluster][apcluster::apclusterK] is called.
#' 
#' @param prc argument needed when K is not null. The algorithm stops if the
#' number of clusters does not deviate more than prc percent from desired value
#' K; set to 0 if you want to have exactly K clusters.
#' 
#' @param bimaxit argument needed when K is not null. maximum number of
#' bisection steps to perform; note that no warning is issued if the number of
#' clusters is still not in the desired range.
#' 
#' @param exact flag indicating whether or not to compute the initial
#' preference range exactly.
#' 
#' 
#' @param algorithm_in_output a `boolean` indicating if the original output
#' of [apcluster][apcluster::apcluster] should be returned in the output
#' (`TRUE` by default, see Value).
#' 
#' @details
#' Based on [apcluster](https://cran.r-project.org/package=apcluster)
#' package ([apcluster][apcluster::apcluster]).
#' 
#' @return
#' A `list` of class `bioregion.clusters` with five slots:
#' \enumerate{
#' \item{**name**: `character` containing the name of the algorithm}
#' \item{**args**: `list` of input arguments as provided by the user}
#' \item{**inputs**: `list` of characteristics of the clustering process}
#' \item{**algorithm**: `list` of all objects associated with the
#'  clustering procedure, such as original cluster objects}
#' \item{**clusters**: `data.frame` containing the clustering results}}
#' 
#' In the `algorithm` slot, if `algorithm_in_output = TRUE`, users can
#' find the output of [apcluster][apcluster::apcluster].
#'
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}),
#' Boris Leroy (\email{leroy.boris@gmail.com}), and
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) 
#' 
#' @seealso  [nhclu_pam]
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
#' @references
#' \insertRef{Frey2007}{bioregion}
#' 
#' @importFrom apcluster apcluster apclusterK
#'         
#' @export

nhclu_affprop <- function(similarity, index = names(similarity)[3],
                          p = NA, q = NA, maxits = 1000, convits = 100,
                          lam = 0.9, # includeSim = FALSE,
                          details = FALSE, nonoise = FALSE, seed = NULL,
                          K = NULL, prc = NULL, bimaxit = NULL, exact = NULL,
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
  
  # Control other arguments
  controls(args = NULL, data = similarity,
           type = "input_conversion_similarity")
  
  if(!is.na(q)){
    if (length(q) > 1) {
      stop("q should be a numeric value in the range [0, 1]")
    }
    if (!is.numeric(q)) {
      stop("q should be a numeric value in the range [0, 1]")
    } else if (q < 0.5) {
      stop("q should be a numeric value in the range [0, 1]")
    }
  }
  
  if(!is.na(p)){
    if(!is.null(p)){
      if (!is.numeric(p)) {
        stop("p should be a number or a vector")
      }
      if (length(p) > 1) {
        if(length(p) != length(unique(c(sim[, 1], sim[, 2])))){
          stop("vector 'p' is shorter than number of sites.")
        }
      }
    }
  }
  
  controls(args = maxits, data = NULL, type = "positive_integer")
  controls(args = convits, data = NULL, type = "positive_integer")
  if (length(lam) > 1) {
    stop("lam should be a numeric value in the range [0.5, 1)")
  }
  if (!is.numeric(lam)) {
    stop("lam should be a numeric value in the range [0.5, 1)")
  } else if (lam < 0.5) {
    stop("lam should be a numeric value in the range [0.5, 1)")
  }
  # controls(args = includeSim, data = NULL, type = "boolean")
  controls(args = details, data = NULL, type = "boolean")
  controls(args = nonoise, data = NULL, type = "boolean")
  if(!is.null(seed)){
    controls(args = seed, data = NULL, type = "strict_positive_integer")
  }
  
  # Argument for desired number of clusters: positive integer, if not null
  # (default value) then we call apcluter::apclusterK() (with argument K)
  if(!is.null(K)){
    controls(args = K, data = NULL, type = "strict_positive_integer")
    if(is.null(prc)){
      stop("When K is not NULL, you need to define prc. prc is a percentage
           value.")
    }
    if(is.null(bimaxit)){
      stop("When K is not NULL, you need to define bimaxit.")
    }
    if(is.null(exact)){
      stop("When K is not NULL, you need to define exact.")
    }
  }
  
  if(!is.null(prc)){
    if(is.null(K)){
      warning("The prc argument will be considered only if K is not set to
            NULL.")
    }
    if (length(prc) > 1) {
      stop("prc should be a numeric value in the range [0, 100]")
    }
    if (!is.numeric(prc)) {
      stop("prc should be a numeric value in the range [0, 100]")
    } else if (prc < 0 | prc > 100) {
      stop("frac should be a numeric value in the range [0, 100]")
    }
  }
  
  if(!is.null(bimaxit)){
    if(is.null(K)){
      warning("The bimaxit argument will be considered only if K is not set to
            NULL.")
    }
    controls(args = bimaxit, data = NULL, type = "strict_positive_integer")
  }
  
  if(!is.null(exact)){
    if(is.null(K)){
      warning("The exact argument will be considered only if K is not set to
            NULL.")
    }
    controls(args = exact, data = NULL, type = "boolean")
  }
  
  sim <- NULL
  
  # 2. Function ---------------------------------------------------------------
  outputs <- list(name = "nhclu_affprop")
  
  outputs$args <- list(index = index,
                       p = p,
                       q = q,
                       maxits = maxits,
                       convits = convits,
                       lam = lam,
                       # includeSim = includeSim,
                       details = details,
                       nonoise = nonoise,
                       seed = seed,
                       algorithm_in_output = algorithm_in_output)
  
  outputs$inputs <- list(bipartite = FALSE,
                         weight = TRUE,
                         pairwise = TRUE,
                         pairwise_metric = ifelse(!inherits(dissimilarity, 
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
  sim_square <- net_to_mat(similarity, weight = TRUE, squared = TRUE,
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
