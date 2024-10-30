#' Hierarchical clustering based on dissimilarity or beta-diversity
#'
#' This function generates a hierarchical tree from a dissimilarity
#' (beta-diversity) `data.frame`, calculates the cophenetic correlation
#' coefficient, and can get clusters from the tree if requested by the user.
#' The function implements randomization of the dissimilarity matrix to
#' generate the tree, with two different methods to generate the final tree.
#' Typically, the dissimilarity `data.frame` is a
#' `bioregion.pairwise.metric` object obtained by running `similarity`
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
#' [hclust][fastcluster::hclust]. Should be one of "ward.D",
#' "ward.D2", "single", "complete", "average"
#' (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or
#' "centroid" (= UPGMC).
#' 
#' @param randomize a `boolean` indicating if the dissimilarity matrix should be
#' randomized, to account for the order of sites in the dissimilarity matrix.
#' 
#' @param n_runs number of trials to randomize the dissimilarity matrix.
#' 
#' @param keep_trials a `boolean` indicating if all random trial results.
#' should be stored in the output object (set to FALSE to save space if your
#' `dissimilarity` object is large).
#' 
#' @param optimal_tree_method a `character` indicating how the final tree
#' should be obtained from all trials. The only option currently is
#' "best", which means the tree with the best cophenetic correlation
#' coefficient will be chosen.
#' 
#' @param n_clust an `integer` or an `integer` vector indicating the number of
#' clusters to be obtained from the hierarchical tree, or the output from
#' [partition_metrics]. Should not be used at the same time as
#' `cut_height`.
#' 
#' @param cut_height a `numeric` vector indicating the height(s) at which the
#' tree should be cut. Should not be used at the same time as `n_clust`.
#' 
#' @param find_h a `boolean` indicating if the height of cut should be found for
#' the requested `n_clust`.
#' 
#' @param h_max a `numeric` indicating the maximum possible tree height for
#' the chosen `index`.
#' 
#' @param h_min a `numeric` indicating the minimum possible height in the tree
#' for the chosen `index`.
#' 
#' @param consensus_p a `numeric`, (only if `optimal_tree_method = "consensus"`, 
#' indicating the threshold proportion of trees that must 
#' support a region/cluster for it to be included in the final consensus tree.
#'  
#' @details
#' The function is based on [hclust][fastcluster::hclust].
#' The default method for the hierarchical tree is `average`, i.e.
#' UPGMA as it has been recommended as the best method to generate a tree
#' from beta diversity dissimilarity \insertCite{Kreft2010}{bioregion}.
#'
#' Clusters can be obtained by two methods:
#' \itemize{
#' \item{Specifying a desired number of clusters in `n_clust`}
#' \item{Specifying one or several heights of cut in `cut_height`}}
#'
#' To find an optimal number of clusters, see [partition_metrics()]
#' 
#' It is important to pay attention to the fact that the order of rows
#' in the input distance matrix
#' influences the tree topology as explained in 
#' \insertRef{Dapporto2013}{bioregion}. To address this, the function generates
#'  multiple trees by randomizing the distance matrix. 
#' Two methods are available to obtain the final tree:
#' \itemize{
#' \item{`optimal_tree_method = "best"`: This method selects the tree with 
#' the highest cophenetic correlation coefficient, representing the best fit 
#' between the hierarchical structure and the original distance matrix. }
#' \item{`optimal_tree_method = "consensus"`: This method constructs a consensus 
#' tree using phylogenetic methods with the function 
#' [consensus][ape::consensus].
#' When using this option, you must set the `consensus_p` parameter, which 
#' indicates 
#' the proportion of trees that must contain a region/cluster for it to be 
#' included 
#' in the final consensus tree. 
#' Consensus trees lack an inherent height because they represent a majority 
#' structure rather than an actual hierarchical clustering. To assign heights, 
#' we use a non-negative least squares method ([nnls.tree][phangorn::nnls.tree]) 
#' based on the initial distance matrix, ensuring that the consensus 
#' tree preserves 
#' approximate distances among clusters.}
#' }
#' AIt is currently unresolved which method is best to use, and so we recommend
#' users to consult the literature.  
#' Consensus trees can offer more stable solutions 
#' in certain contexts \insertRef{Dapporto2013}{bioregion}, but at large 
#' spatial scales, they may yield inconsistent results 
#' (see appendix S2 in \insertRef{Leroy2019}{bioregion}). 
#' Additionally, consensus trees often have a lower cophenetic 
#' correlation coefficient 
#' than the best individual tree, so users should consider these factors when 
#' choosing the method.
#'
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
#' \insertRef{Kreft2010}{bioregion}
#' @author
#' Boris Leroy (\email{leroy.boris@gmail.com}),
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) and
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) 
#' 
#' @seealso [cut_tree]
#' @examples
#' comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
#' 20, 25)
#' rownames(comat) <- paste0("Site",1:20)
#' colnames(comat) <- paste0("Species",1:25)
#'
#' dissim <- dissimilarity(comat, metric = "Simpson")
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
#' # Make multiple cuts
#' tree4 <- hclu_hierarclust(dissim, n_clust = 1:19)
#' 
#' # Change the method to get the final tree (see details)
#' tree5 <- hclu_hierarclust(dissim,
#'                           optimal_tree_method = "consensus",
#'                           n_clust = 10,
#'                           consensus_p = 0.75)
#' 
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
                             h_min = 0,
                             consensus_p = 0.5){
  
  # 1. Controls ---------------------------------------------------------------
  controls(args = NULL, data = dissimilarity, type = "input_nhandhclu")
  if(!inherits(dissimilarity, "dist")){
    controls(args = NULL, data = dissimilarity, type = "input_dissimilarity")
    controls(args = NULL, data = dissimilarity, 
             type = "input_data_frame_nhandhclu")
    controls(args = index, data = dissimilarity, type = "input_net_index")
    net <- dissimilarity
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
    controls(args = NULL, data = dissimilarity, type = "input_dist")
    dist.obj <- dissimilarity
    if(is.null(labels(dist.obj))){
      attr(dist.obj, "Labels") <- paste0(1:attr(dist.obj, "Size"))
      message("No labels detected, they have been assigned automatically.")
    }
  }
  
  controls(args = method, data = NULL, type = "character")
  if(!(method %in% c("ward.D", "ward.D2", "single", "complete", "average",
                          "mcquitty", "median", "centroid" ))){
    stop("Please choose method among the followings values:
ward.D, ward.D2, single, complete, average, mcquitty, median or centroid", 
         call. = FALSE)
  }
  controls(args = randomize, data = NULL, type = "boolean")
  controls(args = n_runs, data = NULL, type = "strict_positive_integer")
  controls(args = keep_trials, data = NULL, type = "boolean")
  controls(args = optimal_tree_method, data = NULL, type = "character")
  if(!(optimal_tree_method %in% c("best", "consensus"))){
    stop("Please choose optimal_tree_method among the followings values:
best, consensus", 
    call. = FALSE)
  }
  
  if(!is.null(n_clust)) {
    if(is.numeric(n_clust)) {
        controls(args = n_clust, data = NULL, 
                 type = "strict_positive_integer_vector")
    } else if(inherits(n_clust, "bioregion.partition.metrics")){
      if(!is.null(n_clust$algorithm$optimal_nb_clusters)) {
        n_clust <- n_clust$algorithm$optimal_nb_clusters
      } else {
        stop("n_clust does not have an optimal number of clusters. Did you
        specify partition_optimisation = TRUE in partition_metrics()?", 
             call. = FALSE)
      }
    } else{
      stop("n_clust must be one of those:
        * an integer determining the number of clusters
        * a vector of integers determining the numbers of clusters for each cut
        * the output from partition_metrics()", 
           call. = FALSE)
    }
    if(!is.null(cut_height)){
      stop("Please provide either n_clust or cut_height, but not both at the
           same time.", 
           call. = FALSE)
    }
  }
  if(!is.null(cut_height)){
    controls(args = cut_height, data = NULL, type = "positive_numeric_vector")
  }
  controls(args = find_h, data = NULL, type = "boolean")
  controls(args = h_min, data = NULL, type = "positive_numeric")
  controls(args = h_max, data = NULL, type = "positive_numeric")
  if(h_min > h_max){
    stop("h_min must be inferior to h_max.")
  }
  controls(args = consensus_p, data = NULL, type = "positive_numeric")
  if(consensus_p < 0.5 | consensus_p > 1) {
    stop("consensus_p must be between 0.5 and 1")
  }
  
  # 2. Function ---------------------------------------------------------------
  outputs <- list(name = "hclu_hierarclust")
  
  # Adding dynamic_tree_cut = FALSE for compatibility with generic functions
  dynamic_tree_cut <- FALSE
  outputs$args <- list(index = index,
                       method = method,
                       randomize = randomize,
                       n_runs = n_runs,
                       optimal_tree_method = optimal_tree_method,
                       keep_trials = keep_trials,
                       n_clust = n_clust,
                       cut_height = cut_height,
                       find_h = find_h,
                       h_max = h_max,
                       h_min = h_min,
                       consensus_p = consensus_p,
                       dynamic_tree_cut = dynamic_tree_cut)
  
  outputs$inputs <- list(bipartite = FALSE,
                         weight = TRUE,
                         pairwise = TRUE,
                         pairwise_metric = ifelse(!inherits(dissimilarity, 
                                                            "dist"), 
                                                  ifelse(is.numeric(index), 
                                                         names(net)[3], index), 
                                                  NA),
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
    } else if (optimal_tree_method == "consensus") {
      # New method - calculate the consensus tree
      if (n_runs < 2) {
        stop("At least two trees are required to calculate a consensus.")
      }
      trees <- lapply(outputs$algorithm$trials,
                      function(trial) trial$hierartree)
      trees <- lapply(trees, 
                      function(tree) ape::as.phylo(tree))
      consensus_tree <- ape::consensus(trees, p = 0.5)
      consensus_tree <- phangorn::nnls.tree(dist.obj, consensus_tree, 
                                            method = "ultrametric",
                                            trace = 0)

      consensus_tree <- ape::multi2di(consensus_tree)
      # Calculate cophenetic correlation coefficient
      coph <- as.matrix(
        stats::cophenetic(consensus_tree))
      consensus_tree <- ape::as.hclust.phylo(consensus_tree)
      outputs$algorithm$final.tree <- consensus_tree
      
      # Correct site order to not mess the correlation
      coph <- coph[match(attr(dist.obj,
                              "Labels"),
                         rownames(coph)),
                   match(attr(dist.obj,
                              "Labels"),
                         colnames(coph))]
      dist.mat <- as.matrix(dist.obj)
      
      outputs$algorithm$final.tree.coph.cor <-
        stats::cor(dist.mat[lower.tri(dist.mat)], coph[lower.tri(coph)],
                   method = "spearman")
      
    }
    
    message(paste0(
      "Final tree has a ",
      round(outputs$algorithm$final.tree.coph.cor, 2),
      " cophenetic correlation coefficient with the initial dissimilarity
      matrix\n"))
    
  } else  {
    outputs$algorithm$final.tree <- fastcluster::hclust(dist.obj,
                                                        method = method)
    
    coph <- as.matrix(stats::cophenetic(outputs$algorithm$final.tree))
    coph <- coph[match(attr(dist.obj, "Labels"), rownames(coph)),
                 match(attr(dist.obj, "Labels"), colnames(coph))]
    dist.mat <- as.matrix(dist.obj)
    
    outputs$algorithm$final.tree.coph.cor <-
      stats::cor(dist.mat[lower.tri(dist.mat)], coph[lower.tri(coph)],
                 method = "spearman")
    
    message(paste0("Output tree has a ",
                   round(outputs$algorithm$final.tree.coph.cor, 2),
                   " cophenetic correlation coefficient with the initial
                   dissimilarity matrix\n"))
  }
  
  class(outputs) <- append("bioregion.clusters", class(outputs))
  
  if(any(!is.null(n_clust) | !is.null(cut_height))){
    outputs <- cut_tree(outputs,
                        n_clust = n_clust,
                        cut_height = cut_height,
                        find_h = find_h,
                        h_max = h_max,
                        h_min = h_min,
                        dynamic_tree_cut = dynamic_tree_cut)
    outputs$inputs$hierarchical <- ifelse(ncol(outputs$clusters) > 2,
                                          TRUE,
                                          FALSE)
  } else {
    outputs$clusters <- NA
    outputs$cluster_info <- NA
    outputs$inputs$hierarchical <- FALSE
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



