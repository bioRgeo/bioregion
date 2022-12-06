#' Hierarchical clustering based on dissimilarity or beta-diversity
#'
#' This function generates a hierarchical tree from a dissimilarity
#' (beta-diversity) `data.frame`, calculates the cophenetic correlation
#' coefficient, and can get clusters from the tree if requested by the user.
#' The function implements randomization of the dissimilarity matrix to
#' generate the tree, with a selection method based on the optimal cophenetic
#' correlation coefficient. Typically, the dissimilarity `data.frame` is a
#' `bioRgeo.pairwise.metric` object obtained by running `similarity`
#' or `similarity` and then `similarity_to_dissimilarity`.
#'
#' @param dissimilarity the output object from [dissimilarity()] or
#'  [similarity_to_dissimilarity()], or a `dist` object. 
#'  If a `data.frame` is used, the first two 
#' columns represent pairs of sites (or any pair of nodes), and the next
#' column(s) are the dissimilarity indices.
#' 
#' @param index name or number of the dissimilarity column to use. By default, 
#' the third column name of `dissimilarity` is used.
#'  
#' @param method name of the hierarchical classification method, as in
#' [fastcluster::hclust()][fastcluster::hclust]. Should be one of `"ward.D"`,
#' `"ward.D2"`, `"single"`, `"complete"`, `"average"`
#' (= UPGMA), `"mcquitty"` (= WPGMA), `"median"` (= WPGMC) or
#' `"centroid"` (= UPGMC).
#' 
#' @param randomize a boolean indicating if the dissimilarity matrix should be
#' randomized, to account for the order of sites in the dissimilarity matrix.
#' @param n_runs number of trials to randomize the dissimilarity matrix.
#' 
#' @param keep_trials a boolean indicating if all random trial results.
#' should be stored in the output object (set to FALSE to save space if your
#' `dissimilarity` object is large).
#' 
#' @param optimal_tree_method a character vector indicating how the final tree
#' should be obtained from all trials. The only option currently is
#' `"best"`, which means the tree with the best cophenetic correlation
#' coefficient will be chosen.
#' 
#' @param n_clust an integer or a vector of integers indicating the number of
#' clusters to be obtained from the hierarchical tree, or the output from
#' [partition_metrics]. Should not be used at the same time as
#' `cut_height`.
#' 
#' @param cut_height a numeric vector indicating the height(s) at which the
#' tree should be cut. Should not be used at the same time as `n_clust`.
#' 
#' @param find_h a boolean indicating if the height of cut should be found for
#' the requested `n_clust`.
#' 
#' @param h_max a numeric indicating the maximum possible tree height for
#' the chosen `index`.
#' 
#' @param h_min a numeric indicating the minimum possible height in the tree
#' for the chosen `index`.
#' 
#' @details
#' The default method for the hierarchical tree is `"average"`, i.e.
#' UPGMA as it has been recommended as the best method to generate a tree
#' from beta diversity dissimilarity \insertCite{Kreft2010}{bioRgeo}
#'
#' Clusters can be obtained by two methods:
#' \itemize{
#' \item{Specifying a desired number of clusters in `n_clust`}
#' \item{Specifying one or several heights of cut in `cut_height`}}
#'
#' To find an optimal number of clusters, see [partition_metrics()]
#'
#' @return
#' A `list` of class `bioRgeo.clusters` with five slots:
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
#' \item{`trials`: a list containing all randomization trials. Each trial
#' contains the dissimilarity matrix, with site order randomized, the
#' associated tree and the cophenetic correlation coefficient (Spearman) for
#' that tree}
#' \item{`final.tree`: a `hclust` object containing the final
#' hierarchical tree to be used}
#' \item{`final.tree.coph.cor`: the cophenetic correlation coefficient
#' between the initial dissimilarity matrix and `final.tree`}
#' }
#'
#' @references
#' \insertRef{Kreft2010}{bioRgeo}
#' @author
#' Boris Leroy (\email{leroy.boris@gmail.com}),
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) and
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) 
#' 
#' @seealso [cut_tree] 
#' @examples
#' \dontrun{
#' comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
#' 20, 25)
#' rownames(comat) <- paste0("Site",1:20)
#' colnames(comat) <- paste0("Species",1:25)
#'
#' dissim <- dissimilarity(comat, metric = "all")
#'
#' # User-defined number of clusters
#' tree1 <- hclu_hierarclust(dissim, n_clust = 5)
#' tree1
#' plot(tree1)
#' str(tree1)
#' tree1$clusters
#' 
#' # User-defined height cut
#' # Only one height
#' tree2 <- hclu_hierarclust(dissim, cut_height = .05)
#' tree2
#' tree2$clusters
#' 
#' # Multiple heights
#' tree3 <- hclu_hierarclust(dissim, cut_height = c(.05, .15, .25))
#' 
#' tree3$clusters # Mind the order of height cuts: from deep to shallow cuts
#' # Info on each partition can be found in table cluster_info
#' tree3$cluster_info
#' plot(tree3)
#' 
#' # Recut the tree afterwards
#' tree3.1 <- cut_tree(tree3, n = 5)
#' 
#' tree4 <- hclu_hierarclust(dissim, n_clust = 1:19)
#' }
#' 
#' @importFrom stats as.dist cophenetic cor
#' @importFrom fastcluster hclust
#' 
#' @export

hclu_hierarclust <- function(dissimilarity,
                             index = names(dissimilarity)[3],
                             method = "average",
                             randomize = TRUE,
                             n_runs = 30,
                             keep_trials = FALSE,
                             optimal_tree_method = "best", # best or consensus
                             n_clust = NULL,
                             cut_height = NULL,
                             find_h = TRUE,
                             h_max = 1,
                             h_min = 0){
  
  # 1. Controls ---------------------------------------------------------------
  if(inherits(dissimilarity, "bioRgeo.pairwise.metric")){
    if(attr(dissimilarity, "type") == "similarity") {
      stop("dissimilarity seems to be a similarity object.
         hclu_hierarclust() should be applied on dissimilarity, not
         similarities.
         Use similarity_to_dissimilarity() before using hclu_hierarclust()")
    }
    if(is.numeric(index)){
      index <- names(dissimilarity)[index]
    }
    if(!(index %in% colnames(dissimilarity))) {
      stop("Argument index should be one of the column names of dissimilarity")
    }
    
  } else if(!any(inherits(dissimilarity, "bioRgeo.pairwise.metric"),
                 inherits(dissimilarity, "dist"))){
    if(is.numeric(index)) {
      index <- names(dissimilarity)[index]
    } 
    if(is.null(index) || !(index %in% colnames(dissimilarity))) {
      stop("dissimilarity is not a bioRgeo.pairwise.metric object, a
           dissimilarity matrix (class dist) or a data.frame with at least 3
           columns (site1, site2, and your dissimilarity index).")
    }
  }
  
  if(!is.null(n_clust)) {
    if(is.numeric(n_clust)) {
      if(any(!(n_clust %% 1 == 0))) {
        stop("n_clust must an integer or a vector of integers determining the
             number of clusters.")
      }
    } else if(inherits(n_clust, "bioRgeo.partition.metrics")){
      if(!is.null(n_clust$algorithm$optimal_nb_clusters)) {
        n_clust <- n_clust$algorithm$optimal_nb_clusters
      } else {
        stop("n_clust does not have an optimal number of clusters. Did you
        specify partition_optimisation = TRUE in partition_metrics()?")
      }
    } else{
      stop("n_clust must be one of those:
        * an integer determining the number of clusters
        * a vector of integers determining the numbers of clusters for each cut
        * the output from partition_metrics()")
    }
    if(!is.null(cut_height)){
      stop("Please provide either n_clust or cut_height, but not both at the
           same time.")
    }
  }
  
  if(!is.character(method) || length(method) != 1 ||
     !(all(method %in% c("ward.D", "ward.D2", "single", "complete", "average",
                         "mcquitty", "median", "centroid" )))){
    stop("method is a character string indicating what hierarchical
         classification method to use. See help for available options.")
  }
  
  if(!is.logical(randomize)){
    stop("randomize must be a Boolean.")
  }
  
  if(!is.numeric(n_runs) || n_runs < 0){
    stop("n_runs must be a positive integer.")
  }
  
  if(!is.logical(keep_trials)){
    stop("keep_trials must be a Boolean.")
  }
  
  if(optimal_tree_method != "best"){
    stop("optimal_tree_method must be a character string. Only available
         option at the moment is best.")
  }
  
  if(!is.null(cut_height)){
    if(!is.numeric(cut_height) || cut_height < 0){
      stop("cut_height must be a positive integer.")
    }
  }
  
  if(!is.logical(find_h)){
    stop("find_h must be a Boolean.")
  }
  
  if(!is.numeric(h_max) || h_max < 0){
    stop("h_max must be a positive integer.")
  }
  
  if(!is.numeric(h_min) || h_min < 0){
    stop("h_min must be a positive integer.")
  }
  
  if(h_min > h_max){
    stop("h_min must be inferior to h_max.")
  }
  
  # 2. Function ---------------------------------------------------------------
  outputs <- list(name = "hierarchical_clustering")
  
  if(!inherits(dissimilarity, "dist")){
    # dist.obj <- .dfToDist(dissimilarity, metric = index)
    dist.obj <- stats::as.dist(
      net_to_mat(dissimilarity[, c(colnames(dissimilarity)[1:2], index)],
                 weight = TRUE, squared = TRUE, symmetrical = TRUE))
    
  } else {
    dist.obj <- dissimilarity
  }
  
  # Adding dynamic_tree_cut = FALSE for compatibility with generic functions
  dynamic_tree_cut = FALSE
  outputs$args <- list(index = index,
                       method = method,
                       randomize = randomize,
                       n_runs = n_runs,
                       optimal_tree_method = optimal_tree_method,
                       n_clust = n_clust,
                       cut_height = cut_height,
                       find_h = find_h,
                       h_max = h_max,
                       h_min = h_min,
                       dynamic_tree_cut = dynamic_tree_cut)
  
  outputs$inputs <- list(bipartite = FALSE,
                         weight = TRUE,
                         pairwise_metric = TRUE,
                         dissimilarity = TRUE,
                         nb_sites = attr(dist.obj, "Size"))
  
  # outputs$dist.matrix <- dist.obj
  
  if(randomize) {
    message(paste0("Randomizing the dissimilarity matrix with ", n_runs,
                   " trials"))
    for(run in 1:n_runs){
      outputs$algorithm$trials[[run]] <- list()
      outputs$algorithm$trials[[run]]$dist.matrix <-
        .randomizeDistance(dist.obj)
      
      # Compute hierarchical tree
      outputs$algorithm$trials[[run]]$hierartree <-
        fastcluster::hclust(outputs$algorithm$trials[[run]]$dist.matrix,
                            method = method)
      
      # Calculate cophenetic correlation coefficient
      coph <- as.matrix(
        stats::cophenetic(outputs$algorithm$trials[[run]]$hierartree))
      # Correct site order to not mess the correlation
      coph <- coph[match(attr(outputs$algorithm$trials[[run]]$dist.matrix,
                              "Labels"),
                         rownames(coph)),
                   match(attr(outputs$algorithm$trials[[run]]$dist.matrix,
                              "Labels"),
                         colnames(coph))]
      dist.mat <- as.matrix(outputs$algorithm$trials[[run]]$dist.matrix)
      
      outputs$algorithm$trials[[run]]$cophcor <-
        stats::cor(dist.mat[lower.tri(dist.mat)], coph[lower.tri(coph)],
                   method = "spearman")
    }
    
    if(optimal_tree_method == "best")  {
      coph.coeffs <- sapply(1:n_runs, function(x)
      {
        outputs$algorithm$trials[[x]]$cophcor
      })
      
      message(paste0(" -- range of cophenetic correlation coefficients among
                     trials: ", round(min(coph.coeffs), 2),
                     " - ", round(max(coph.coeffs), 2)))
      
      # There might be multiple trees with the highest cophenetic correlation
      # coefficient, so we arbitrarily take the first one
      best.run <- which(coph.coeffs == max(coph.coeffs))[1]
      
      outputs$algorithm$final.tree <-
        outputs$algorithm$trials[[best.run]]$hierartree
      outputs$algorithm$final.tree.coph.cor <- max(coph.coeffs)
    }
    
    message(paste0(
      "Optimal tree has a ",
      round(outputs$algorithm$final.tree.coph.cor, 2),
      " cophenetic correlation coefficient with the initial dissimilarity
      matrix\n"))
    
  } else  {
    outputs$algorithm$final.tree <- fastcluster::hclust(dist.obj,
                                                        method = method)
    
    coph <- as.matrix(stats::cophenetic(outputs$algorithm$final.tree))
    coph <- coph[match(attr(dist.obj, "Labels"),
                       rownames(coph)),
                 match(attr(dist.obj, "Labels"),
                       colnames(coph))]
    dist.mat <- as.matrix(dist.obj)
    
    outputs$algorithm$final.tree.coph.cor <-
      stats::cor(dist.mat[lower.tri(dist.mat)], coph[lower.tri(coph)],
                 method = "spearman")
    
    message(paste0("Output tree has a ",
                   round(outputs$algorithm$final.tree.coph.cor, 2),
                   " cophenetic correlation coefficient with the initial
                   dissimilarity matrix\n"))
  }
  
  class(outputs) <- append("bioRgeo.clusters", class(outputs))
  
  if(any(!is.null(n_clust) | !is.null(cut_height))){
    outputs <- cut_tree(outputs,
                        n_clust = n_clust,
                        cut_height = cut_height,
                        find_h = find_h,
                        h_max = h_max,
                        h_min = h_min,
                        dynamic_tree_cut = dynamic_tree_cut)
  } else{
    outputs$clusters <- NA
    outputs$cluster_info <- NA
  }
  
  if(!keep_trials){
    outputs$algorithm$trials <- "Trials not stored in output"
  }
  
  return(outputs)
}

.randomizeDistance <- function(distmatrix){
  distmatrix <- as.matrix(distmatrix)
  ord <- sample(rownames(distmatrix))
  return(stats::as.dist(distmatrix[ord, ord]))
}
