#' Hierarchical clustering based on dissimilarity or beta-diversity
#'
#' This function generates a hierarchical tree from a dissimilarity
#' (beta-diversity) `data.frame`, calculates the cophenetic correlation
#' coefficient, and optionally retrieves clusters from the tree upon user 
#' request. The function includes a randomization process for the dissimilarity 
#' matrix to generate the tree, with two methods available for constructing the 
#' final tree. Typically, the dissimilarity `data.frame` is a
#' `bioregion.pairwise` object obtained by running `similarity`,
#' or by running `similarity` followed by `similarity_to_dissimilarity`.
#'
#' @param dissimilarity The output object from [dissimilarity()] or
#'  [similarity_to_dissimilarity()], or a `dist` object. 
#'  If a `data.frame` is used, the first two columns represent pairs of sites 
#'  (or any pair of nodes), and the subsequent column(s) contain the 
#'  dissimilarity indices.
#' 
#' @param index The name or number of the dissimilarity column to use. By 
#' default, the third column name of `dissimilarity` is used.
#' 
#' @param method The name of the hierarchical classification method, as in
#' [hclust][fastcluster::hclust]. Should be one of `"ward.D"`,
#' `"ward.D2"`, `"single"`, `"complete"`, `"average"`
#' (= UPGMA), `"mcquitty"` (= WPGMA), `"median"` (= WPGMC), or
#' `"centroid"` (= UPGMC).
#' 
#' @param randomize A `boolean` indicating whether the dissimilarity matrix 
#' should be randomized to account for the order of sites in the dissimilarity
#'  matrix.
#' 
#' @param n_runs The number of trials for randomizing the dissimilarity matrix.
#' 
#' @param keep_trials A `boolean` indicating whether all random trial results
#' should be stored in the output object. Set to `FALSE` to save space if your
#' `dissimilarity` object is large. Note that this cannot be set to `TRUE` if
#' `optimal_tree_method = "iterative_consensus_tree"`.
#' 
#' @param optimal_tree_method A `character` string indicating how the final tree
#' should be obtained from all trials. Possible values are 
#' `"iterative_consensus_tree"` (default), `"best"`, and `"consensus"`. 
#' **We recommend `"iterative_consensus_tree"`. See Details.** 
#' 
#' @param n_clust An `integer` vector or a single `integer` indicating the 
#' number of clusters to be obtained from the hierarchical tree, or the output 
#' from [bioregionalization_metrics]. This parameter should not be used 
#' simultaneously with `cut_height`.
#' 
#' @param cut_height A `numeric` vector indicating the height(s) at which the
#' tree should be cut. This parameter should not be used simultaneously with 
#' `n_clust`.
#' 
#' @param find_h A `boolean` indicating whether the height of the cut should be 
#' found for the requested `n_clust`.
#' 
#' @param h_max A `numeric` value indicating the maximum possible tree height 
#' for the chosen `index`.
#' 
#' @param h_min A `numeric` value indicating the minimum possible height in the
#'  tree for the chosen `index`.
#' 
#' @param consensus_p A `numeric` value (applicable only if 
#' `optimal_tree_method = "consensus"`) indicating the threshold proportion of 
#' trees that must support a region/cluster for it to be included in the final 
#' consensus tree.
#'  
#' @param verbose A `boolean` (applicable only if 
#' `optimal_tree_method = "iterative_consensus_tree"`) indicating whether to 
#' display progress messages. Set to `FALSE` to suppress these messages.
#' 
#' @return
#' A `list` of class `bioregion.clusters` with five slots:
#' \enumerate{
#' \item{**name**: A `character` string containing the name of the algorithm.}
#' \item{**args**: A `list` of input arguments as provided by the user.}
#' \item{**inputs**: A `list` describing the characteristics of the clustering process.}
#' \item{**algorithm**: A `list` containing all objects associated with the
#'  clustering procedure, such as the original cluster objects.}
#' \item{**clusters**: A `data.frame` containing the clustering results.}}
#'
#' In the `algorithm` slot, users can find the following elements:
#'
#' \itemize{
#' \item{`trials`: A list containing all randomization trials. Each trial
#' includes the dissimilarity matrix with randomized site order, the
#' associated tree, and the cophenetic correlation coefficient (Spearman) for
#' that tree.}
#' \item{`final.tree`: An `hclust` object representing the final
#' hierarchical tree to be used.}
#' \item{`final.tree.coph.cor`: The cophenetic correlation coefficient
#' between the initial dissimilarity matrix and the `final.tree`.}
#' }
#'  
#' @details
#' The function is based on [hclust][fastcluster::hclust].
#' The default method for the hierarchical tree is `average`, i.e.
#' UPGMA as it has been recommended as the best method to generate a tree
#' from beta diversity dissimilarity (Kreft & Jetz, 2010).
#'
#' Clusters can be obtained by two methods:
#' \itemize{
#' \item{Specifying a desired number of clusters in `n_clust`}
#' \item{Specifying one or several heights of cut in `cut_height`}}
#'
#' To find an optimal number of clusters, see [bioregionalization_metrics()]
#' 
#' It is important to pay attention to the fact that the order of rows
#' in the input distance matrix influences the tree topology as explained in 
#' Dapporto (2013). To address this, the function generates multiple trees by 
#' randomizing the distance matrix. 
#' 
#' Two methods are available to obtain the final tree:
#' \itemize{
#' 
#' \item{`optimal_tree_method = "iterative_consensus_tree"`: The Iterative 
#' Hierarchical Consensus Tree (IHCT) method reconstructs a consensus tree by 
#' iteratively splitting the dataset into two subclusters based on the pairwise 
#' dissimilarity of sites across `n_runs` trees based on `n_runs` randomizations
#' of the distance matrix. At each iteration, it 
#' identifies the majority membership of sites into two stable groups across
#' all trees,
#' calculates the height based on the selected linkage method (`method`),
#' and enforces monotonic constraints on 
#' node heights to produce a coherent tree structure. 
#' This approach provides a robust, hierarchical representation of site 
#' relationships, balancing 
#' cluster stability and hierarchical constraints.}
#' 
#' \item{`optimal_tree_method = "best"`: This method selects one tree among with 
#' the highest cophenetic correlation coefficient, representing the best fit 
#' between the hierarchical structure and the original distance matrix. }
#' 
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
#' 
#' We recommend using the `"iterative_consensus_tree"` as all the branches of
#' this tree will always reflect the majority decision among many randomized 
#' versions of the distance matrix. This method is inspired by 
#' Dapporto et al. (2015), which also used the majority decision
#' among many randomized versions of the distance matrix, but it expands it 
#' to reconstruct the entire topology of the tree iteratively. 
#' 
#' We do not recommend using the basic `consensus` method because in many 
#' contexts it provides inconsistent results, with a meaningless tree topology
#' and a very low cophenetic correlation coefficient. 
#' 
#' For a fast exploration of the tree, we recommend using the `best` method
#' which will only select the tree with the highest cophenetic correlation
#' coefficient among all randomized versions of the distance matrix. 
#'
#' @references
#' Kreft H & Jetz W (2010) A framework for delineating biogeographical regions
#' based on species distributions. \emph{Journal of Biogeography} 37, 2029-2053.
#' 
#' Dapporto L, Ramazzotti M, Fattorini S, Talavera G, Vila R & Dennis, RLH 
#' (2013) Recluster: an unbiased clustering procedure for beta-diversity 
#' turnover. \emph{Ecography} 36, 1070--1075.
#' 
#' Dapporto L, Ciolli G, Dennis RLH, Fox R & Shreeve TG (2015) A new procedure 
#' for extrapolating turnover regionalization at mid-small spatial scales, 
#' tested on British butterflies. \emph{Methods in Ecology and Evolution} 6
#' , 1287--1297. 
#' 
#' @seealso
#' For more details illustrated with a practical example, 
#' see the vignette: 
#' \url{https://biorgeo.github.io/bioregion/articles/a4_1_hierarchical_clustering.html}.
#' 
#' Associated functions: 
#' [cut_tree]
#' 
#' @author
#' Boris Leroy (\email{leroy.boris@gmail.com}) \cr
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) \cr
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) 
#' 
#' @examples
#' comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
#' 20, 25)
#' rownames(comat) <- paste0("Site",1:20)
#' colnames(comat) <- paste0("Species",1:25)
#'
#' dissim <- dissimilarity(comat, metric = "Simpson")
#'
#' # User-defined number of clusters
#' tree1 <- hclu_hierarclust(dissim, 
#'                           n_clust = 5)
#' tree1
#' plot(tree1)
#' str(tree1)
#' tree1$clusters
#' 
#' # User-defined height cut
#' # Only one height
#' tree2 <- hclu_hierarclust(dissim, 
#'                           cut_height = .05)
#' tree2
#' tree2$clusters
#' 
#' # Multiple heights
#' tree3 <- hclu_hierarclust(dissim, 
#'                           cut_height = c(.05, .15, .25))
#' 
#' tree3$clusters # Mind the order of height cuts: from deep to shallow cuts
#' # Info on each partition can be found in table cluster_info
#' tree3$cluster_info
#' plot(tree3)
#' 
#' @importFrom stats as.dist cophenetic cor
#' @importFrom fastcluster hclust
#' 
#' @export
hclu_hierarclust <- function(dissimilarity,
                             index = names(dissimilarity)[3],
                             method = "average",
                             randomize = TRUE,
                             n_runs = 100,
                             keep_trials = FALSE,
                             optimal_tree_method = "iterative_consensus_tree", 
                             n_clust = NULL,
                             cut_height = NULL,
                             find_h = TRUE,
                             h_max = 1,
                             h_min = 0,
                             consensus_p = 0.5,
                             verbose = TRUE){
  
# TODO: Add show_hierarchy to hclu_hierarclust AND cut_tree
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
    dissimilarity <- mat_to_net(as.matrix(dissimilarity), weight = TRUE)
    if(is.null(index)) {
      colnames(dissimilarity)[3] <- index <- "Dissimilarity"
    } else {
      colnames(dissimilarity)[3] <- index
    }
  }

  
  controls(args = method, data = NULL, type = "character")
  if(!(method %in% c("ward.D", "ward.D2", "single", "complete", "average",
                          "mcquitty", "median", "centroid" ))){
    stop(paste0("Please choose method from the following:\n",
                "ward.D, ward.D2, single, complete, average, mcquitty, median ", 
                "or centroid"), 
         call. = FALSE)
  }
  controls(args = randomize, data = NULL, type = "boolean")
  controls(args = n_runs, data = NULL, type = "strict_positive_integer")
  controls(args = keep_trials, data = NULL, type = "boolean")
  controls(args = optimal_tree_method, data = NULL, type = "character")
  if(!(optimal_tree_method %in% c("iterative_consensus_tree",
                                  "best", "consensus"))){
    stop(paste0("Please choose optimal_tree_method from the following:\n",
                "iterative_consensus_tree, best or consensus"), 
    call. = FALSE)
  }

  if(!is.null(n_clust)) {
    if(is.numeric(n_clust)) {
        controls(args = n_clust, 
                 data = NULL, 
                 type = "strict_positive_integer_vector")
    } else if(inherits(n_clust, "bioregion.bioregionalization.metrics")){
      if(!is.null(n_clust$algorithm$optimal_nb_clusters)) {
        n_clust <- n_clust$algorithm$optimal_nb_clusters
      } else {
        stop(paste0("n_clust does not have an optimal number of clusters. ",
                    "Did you specify partition_optimisation = TRUE in ",
                    "bioregionalization_metrics()?"), 
             call. = FALSE)
      }
    } else{
      stop("n_clust must be one of those:
        * an integer determining the number of clusters
        * a vector of integers determining the numbers of clusters for each cut
        * the output from bioregionalization_metrics()", 
           call. = FALSE)
    }
    if(!is.null(cut_height)){
      stop(paste0("Please provide either n_clust or cut_height, ",
                  "but not both at the same time."), 
           call. = FALSE)
    }
  }
  if(!is.null(cut_height)){
    controls(args = cut_height, data = NULL, type = "positive_numeric_vector")
  }
  controls(args = find_h, data = NULL, type = "boolean")
  if(find_h){
    controls(args = h_min, data = NULL, type = "positive_numeric")
    controls(args = h_max, data = NULL, type = "positive_numeric")
    if(h_min > h_max){
      stop("h_min must be inferior to h_max.",
           call. = FALSE)
    }
  }
  controls(args = consensus_p, data = NULL, type = "positive_numeric")
  if(consensus_p < 0.5 | consensus_p > 1) {
    stop("consensus_p must be between 0.5 and 1.",
         call. = FALSE)
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
                       verbose = verbose,
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
    if(optimal_tree_method == "iterative_consensus_tree") {
      message(paste0("Building the iterative hierarchical consensus tree...",
                     " Note that this",
                     " process can take time especially if you have a lot of",
                     " sites."))

      if(method == "mcquitty") {
        warning("mcquitty (WPGMA) method may not be properly implemented",
                " in Iterative Hierarchical Tree Construction (IHCT), because of the ",
                "hybrid divise-agglomerative nature of IHCT. ",
                "In WPGMA, heights are updated iteratively from bottom to top as ",
                "the tree is constructed. ",
                "In IHCT, divisions are created from top to bottom, based on a ",
                "majority decision among many trees. ", 
                "Hence, it is not possible to exactly compute WPGMA calculations ",
                "with IHCT - final height calculations are approximated into UPGMA.")
      }
      
      consensus_tree <- iterative_consensus_tree(dissimilarity, 
                                                 sites = unique(c(dissimilarity[, 1],
                                                                  dissimilarity[, 2])), 
                                                 index = index,
                                                 method = method,
                                                 depth = 1, 
                                                 # tree_structure = list(), 
                                                 previous_height = Inf, 
                                                 verbose = verbose,
                                                 n_runs = n_runs,
                                                 monotonicity_direction = "bottom-up")
      consensus_tree <- reconstruct_hclust(consensus_tree)
      
      # Compute hierarchical tree
      outputs$algorithm$final.tree <- consensus_tree
      
      evals <- tree_eval(consensus_tree,
                         dist.obj)
      
      outputs$algorithm$final.tree.coph.cor <- evals$cophcor
      # outputs$algorithm$final.tree.2norm <- evals$norm2
      outputs$algorithm$final.tree.msd <- evals$msd
      

    } else {
      message(paste0("Randomizing the dissimilarity matrix with ", 
                     n_runs,
                     " trials"))
      
      results <- vector("list", n_runs)
      
      for (run in 1:n_runs) {
        trial <- list()

        trial$dist.matrix <- .randomizeDistance(dist.obj)

        trial$hierartree <- fastcluster::hclust(trial$dist.matrix, method = method)

        evals <- tree_eval(trial$hierartree, trial$dist.matrix)
        trial$cophcor <- evals$cophcor
        # trial$`2norm` <- evals$norm2
        trial$msd <- evals$msd

        results[[run]] <- trial
      }
      
      if (optimal_tree_method == "best") {
        coph.coeffs <- sapply(results, function(x) x$cophcor)
        
        message(paste0(" -- range of cophenetic correlation coefficients ",
                       "among trials: ",
                       round(min(coph.coeffs), 4), 
                       " - ", 
                       round(max(coph.coeffs), 4)))

        best.run <- which.max(coph.coeffs)
        final.tree <- results[[best.run]]$hierartree
        final.tree.metrics <- results[[best.run]]
        
      } else if (optimal_tree_method == "consensus") {
        if (n_runs < 2) {
          stop("At least two trees are required to calculate a consensus.",
               call. = FALSE)
        }

        trees <- lapply(results, function(trial) ape::as.phylo(trial$hierartree))

        consensus_tree <- ape::consensus(trees, p = 0.5)
        consensus_tree <- phangorn::nnls.tree(dist.obj, consensus_tree, method = "ultrametric", trace = 0)
        consensus_tree <- ape::multi2di(consensus_tree)

        tree_ape_for_coph <- consensus_tree
        consensus_tree <- ape::as.hclust.phylo(consensus_tree)
        
        final.tree <- consensus_tree
        evals <- tree_eval(tree_ape_for_coph, dist.obj)
        final.tree.metrics <- list(cophcor = evals$cophcor, 
                                   # `2norm` = evals$norm2, 
                                   msd = evals$msd)
      }

      outputs$algorithm$final.tree <- final.tree
      outputs$algorithm$final.tree.coph.cor <- final.tree.metrics$cophcor
      # outputs$algorithm$final.tree.2norm <- final.tree.metrics$`2norm`
      outputs$algorithm$final.tree.msd <- final.tree.metrics$msd

    }

    
    message(paste0("\nFinal tree has a ",
                   round(outputs$algorithm$final.tree.coph.cor, 4),
                   " cophenetic correlation coefficient with the initial ",
                   "dissimilarity matrix\n"))
    
  } else  {
    outputs$algorithm$final.tree <- fastcluster::hclust(dist.obj,
                                                        method = method)
    
    
    evals <- tree_eval(outputs$algorithm$final.tree,
                       dist.obj)
    
    outputs$algorithm$final.tree.coph.cor <- evals$cophcor
    # outputs$algorithm$final.tree.2norm <- evals$norm2
    outputs$algorithm$final.tree.msd <- evals$msd
    
    message(paste0("Output tree has a ",
                   round(outputs$algorithm$final.tree.coph.cor, 4),
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

# squared_2norm <- function(D, U) {
#   diff_matrix <- D - U
#   singular_values <- svd(diff_matrix)$d
#   return(max(singular_values)^2)
# }
# 
# tree_eval <- function(tree, dist.obj, method = "pearson") {
#   coph <- as.matrix(
#     stats::cophenetic(tree))
#   
#   # Correct site order to not mess the correlation
#   coph <- coph[match(attr(dist.obj,
#                           "Labels"),
#                      rownames(coph)),
#                match(attr(dist.obj,
#                           "Labels"),
#                      colnames(coph))]
#   dist.mat <- as.matrix(dist.obj)
#   
#   return(list(cophcor = stats::cor(dist.mat[lower.tri(dist.mat)],
#                                    coph[lower.tri(coph)],
#                                    method = method),
#               norm2 = squared_2norm (dist.mat, coph),
#               msd = mean((dist.mat - coph)^2)
#   ))
# 
#   # cophcor: Sokal & Rohlf 1962 Taxon
#   # norm2: Mérigot et al. 2010 Ecology
#   # msd: Maire et al. 2015 GEB
# }


# Try to optimize tree_eval for faster calculations
tree_eval <- function(tree, dist.obj, method = "pearson") {
  coph <- as.matrix(stats::cophenetic(tree))

  labels <- attr(dist.obj, "Labels")
  coph <- coph[match(labels, rownames(coph)), 
               match(labels, colnames(coph))]
  
  dist.mat <- as.matrix(dist.obj)
  
  lower_tri_idx <- lower.tri(dist.mat)

  cophcor <- stats::cor(dist.mat[lower_tri_idx],
                        coph[lower_tri_idx], method = method)
  diff_matrix <- dist.mat - coph
  
  # singular_values <- svd(diff_matrix, nu = 0, nv = 0)$d
  # norm2 <- max(singular_values)^2
  
  msd <- mean(diff_matrix[lower_tri_idx]^2)
  #   # cophcor: Sokal & Rohlf 1962 Taxon
  #   # norm2: Mérigot et al. 2010 Ecology
  #   # msd: Maire et al. 2015 GEB
  return(list(cophcor = cophcor, 
              # norm2 = norm2,
              msd = msd))
}
