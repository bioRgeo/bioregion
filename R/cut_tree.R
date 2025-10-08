#' Cut a hierarchical tree
#'
#' This function is designed to work on a hierarchical tree and cut it 
#' at user-selected heights. It works with outputs from either
#' `hclu_hierarclust` or `hclust` objects. The function allows for cutting 
#' the tree based on the chosen number(s) of clusters or specified height(s). 
#' Additionally, it includes a procedure to automatically determine the cutting 
#' height for the requested number(s) of clusters.
#'
#' @param tree A `bioregion.hierar.tree` or an `hclust` object.
#' 
#' @param n_clust An `integer` vector or a single `integer` indicating the 
#' number of clusters to be obtained from the hierarchical tree, or the output 
#' from [bioregionalization_metrics()]. This should not be used concurrently 
#' with `cut_height`.
#' 
#' @param cut_height A `numeric` vector specifying the height(s) at which the
#' tree should be cut. This should not be used concurrently with `n_clust` or
#' `optim_method`.
#' 
#' @param find_h A `boolean` indicating whether the cutting height should be 
#' determined for the requested `n_clust`.
#' 
#' @param h_max A `numeric` value indicating the maximum possible tree height 
#' for determining the cutting height when `find_h = TRUE`.
#' 
#' @param h_min A `numeric` value specifying the minimum possible height in the 
#' tree for determining the cutting height when `find_h = TRUE`.
#' 
#' @param dynamic_tree_cut A `boolean` indicating whether the dynamic tree cut 
#' method should be used. If `TRUE`, `n_clust` and `cut_height` are ignored.
#' 
#' @param dynamic_method A `character` string specifying the method to be used
#' for dynamically cutting the tree: either `"tree"` (clusters searched only
#' within the tree) or `"hybrid"` (clusters searched in both the tree and the 
#' dissimilarity matrix).
#' 
#' @param dynamic_minClusterSize An `integer` indicating the minimum cluster 
#' size for the dynamic tree cut method (see
#' [dynamicTreeCut::cutreeDynamic()][dynamicTreeCut::cutreeDynamic]).
#' 
#' @param dissimilarity Relevant only if `dynamic_method = "hybrid"`. Provide 
#' the dissimilarity `data.frame` used to build the `tree`.
#' 
#' @param show_hierarchy A `boolean` specifying if the hierarchy of clusters
#' should be identifiable in the outputs (`FALSE` by default).
#' 
#' @param ... Additional arguments passed to
#' [dynamicTreeCut::cutreeDynamic()][dynamicTreeCut::cutreeDynamic] to
#' customize the dynamic tree cut method.
#' 
#' @return 
#' If `tree` is an output from [hclu_hierarclust()], the same
#' object is returned with updated content (i.e., `args` and `clusters`). If
#' `tree` is an `hclust` object, a `data.frame` containing the clusters is
#' returned.
#'
#' @details
#' The function supports two main methods for cutting the tree. First, the tree 
#' can be cut at a uniform height (specified by `cut_height` or determined 
#' automatically for the requested `n_clust`). Second, the dynamic tree cut 
#' method (Langfelder et al., 2008) can be applied, which adapts to the shape 
#' of branches in the tree, cutting at varying heights based on cluster 
#' positions.
#'
#' The dynamic tree cut method has two variants:
#' \itemize{
#' \item{The tree-based variant (`dynamic_method = "tree"`) uses a top-down 
#' approach, relying solely on the tree and the order of clustered objects.}
#' \item{The hybrid variant (`dynamic_method = "hybrid"`) employs a bottom-up 
#' approach, leveraging both the tree and the dissimilarity matrix to identify 
#' clusters based on dissimilarity among sites. This approach is useful for 
#' detecting outliers within clusters.}
#' }
#'
#' @note
#' The `find_h` argument is ignored if `dynamic_tree_cut = TRUE`, 
#' as cutting heights cannot be determined in this case.

#' 
#' @references 
#' Langfelder P, Zhang B & Horvath S (2008) Defining clusters from a
#' hierarchical cluster tree: the Dynamic Tree Cut package for R.
#' \emph{BIOINFORMATICS} 24, 719-720.
#' 
#' @seealso 
#' For more details illustrated with a practical example, 
#' see the vignette: 
#' \url{https://biorgeo.github.io/bioregion/articles/a4_1_hierarchical_clustering.html}.
#' 
#' Associated functions: 
#' [hclu_hierarclust]
#'
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) \cr
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) \cr
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' 
#' @examples
#' comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
#' 20, 25)
#' rownames(comat) <- paste0("Site", 1:20)
#' colnames(comat) <- paste0("Species", 1:25)
#'
#' simil <- similarity(comat, metric = "all")
#' dissimilarity <- similarity_to_dissimilarity(simil)
#'
#' # User-defined number of clusters
#' tree1 <- hclu_hierarclust(dissimilarity,
#'                           n_clust = 5)
#' tree2 <- cut_tree(tree1, cut_height = .05)
#' tree3 <- cut_tree(tree1, n_clust = c(3, 5, 10))
#' tree4 <- cut_tree(tree1, cut_height = c(.05, .1, .15, .2, .25))
#' tree5 <- cut_tree(tree1, n_clust = c(3, 5, 10), find_h = FALSE)
#'
#' hclust_tree <- tree2$algorithm$final.tree
#' clusters_2 <- cut_tree(hclust_tree, n_clust = 10)
#'
#' cluster_dynamic <- cut_tree(tree1, dynamic_tree_cut = TRUE,
#'                             dissimilarity = dissimilarity)
#' 
#' @importFrom stats as.dist cutree na.omit
#' @importFrom dynamicTreeCut cutreeDynamic
#' 
#'@export
cut_tree <- function(tree,
                     n_clust = NULL,
                     cut_height = NULL,
                     find_h = TRUE,
                     h_max = 1,
                     h_min = 0,
                     dynamic_tree_cut = FALSE,
                     dynamic_method = "tree",
                     dynamic_minClusterSize = 5,
                     dissimilarity = NULL,
                     show_hierarchy = FALSE,
                     ...){
  
  # Control n_clust
  if(!is.null(n_clust)){
    if(is.numeric(n_clust)){
      if(any(!(n_clust %% 1 == 0))) {
        stop(paste0("n_clust must an integer or a vector of ",
                    "integers determining the number of clusters."), 
             call. = FALSE)
      }
    } else if(inherits(n_clust, "bioregion.partition.metrics")){
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
  
  # Control cut_height
  if(!is.null(cut_height)){
    controls(args = cut_height, data = NULL, type = "positive_numeric_vector")
  }
  
  # Control find_h, h_min and h_max
  controls(args = find_h, data = NULL, type = "boolean")
  if(find_h){
    controls(args = h_min, data = NULL, type = "positive_numeric")
    controls(args = h_max, data = NULL, type = "positive_numeric")
    if(h_min > h_max){
      stop("h_min must be inferior to h_max.",
           call. = FALSE)
    }
  }
  
  # Control show_hierarchy
  controls(args = show_hierarchy, data = NULL, type = "boolean")
  
  # Control dynamic_tree_cut, dynamic_method, dynamic_minClusterSize 
  # and dissimilarity
  controls(args = dynamic_tree_cut, data = NULL, type = "boolean")
  if(dynamic_tree_cut){
    if(!is.null(n_clust)){
      message(paste0("The dynamic tree cut method was requested, ",
                     "argument n_clust will be ignored."))
      n_clust <- NULL
    }
    if(!is.null(cut_height)){
      message(paste0("The dynamic tree cut method was requested, ",
                     "argument cut_height will be ignored."))
      cut_height <- NULL
    }
    
    controls(args = dynamic_method, data = NULL, type = "character")
    if(!(dynamic_method %in% c("tree", "hybrid"))){
      stop(paste0("Please choose dynamic_method from the following:\n",
                  "tree or hybrid"), 
           call. = FALSE)
    }
    if(dynamic_method == "hybrid"){
      if(inherits(dissimilarity, "bioregion.pairwise"))
      {
        if(attr(dissimilarity, "type") == "similarity")
        {
          stop(paste0("dissimilarity seems to be a similarity object. ",
                      "hclu_hierarclust() should be applied on dissimilarity, ",
                      "not similarity. Use similarity_to_dissimilarity() ",
                      "before using hclu_hierarclust()"), 
               call. = FALSE)
        }
        index <- colnames(dissimilarity)[3]
        
        dist_matrix <- stats::as.dist(
          net_to_mat(
            dissimilarity[, c(1, 2,
                              which(colnames(dissimilarity) == index))],
            weight = TRUE, squared = TRUE, symmetrical = TRUE))
        
        
      } else if(!any(inherits(dissimilarity, "bioregion.pairwise"),
                     inherits(dissimilarity, "dist"))){
        if(!is.data.frame(dissimilarity)){
          stop(paste0("dissimilarity is not a bioregion.pairwise ",
                      "object, a dissimilarity matrix (class dist) or a ",
                      "data.frame with at least 3 columns (site1, site2, and ",
                      "your dissimilarity index)"), 
               call. = FALSE)
        }
        if(ncol(dissimilarity) != 3){
          stop(paste0("dissimilarity is not a bioregion.pairwise ",
                      "object, a dissimilarity matrix (class dist) or a ",
                      "data.frame with at least 3 columns (site1, site2, and ",
                      "your dissimilarity index)"), 
               call. = FALSE)
        }
        dist_matrix <- stats::as.dist(
          net_to_mat(
            dissimilarity[, c(1, 2,
                              which(colnames(dissimilarity) == index))],
            weight = TRUE, squared = TRUE, symmetrical = TRUE))
        
      } else{
        dist_matrix <- dissimilarity
      }
    }
  }
  
  arg_added <- list(...)
  if(inherits(tree, "bioregion.clusters")){
    if (tree$name == "hclu_hierarclust") {
      cur.tree <- tree$algorithm$final.tree
      # Update args
      tree$args[c("n_clust", "cut_height", "find_h", "h_max", "h_min",
                  "dynamic_tree_cut", "show_hierarchy")] <-
        list(n_clust, cut_height, find_h, h_max, h_min, dynamic_tree_cut, 
             show_hierarchy)
      
      if(dynamic_tree_cut){
        tree$args[c("dynamic_method",
                    "dynamic_minClusterSize")] <-
          list(dynamic_method,  dynamic_minClusterSize)
        if(length(arg_added))
        {
          tree$args[names(arg_added)] <- arg_added
        }
      }
    }
  } else if (inherits(tree, "hclust"))
  {
    cur.tree <- tree
  } else{
    stop(paste0("This function is designed to work either on outputs ",
                "from hclu_hierarclust() or hclust objects."), 
         call. = FALSE)
  }
  
  output_cut_height <- NULL
  output_n_clust <- NULL
  if(dynamic_tree_cut){
    clusters <- data.frame(site = tree$algorithm$final.tree$labels,
                           cluster = dynamicTreeCut::cutreeDynamic(
                             tree$algorithm$final.tree,
                             method = dynamic_method,
                             minClusterSize = dynamic_minClusterSize,
                             distM = dist_matrix,
                             ...))
    
    # Set NAs for unassigned sites
    if(any(clusters$cluster == 0)){
      message(paste0("Some sites were not assigned to any cluster. ",
                     "They will have a NA in the cluster data.frame."))
      clusters$cluster[which(clusters$cluster == 0)] <- NA
    }
    output_n_clust <- length(unique(stats::na.omit(clusters$cluster)))
  } else if(!is.null(n_clust))
  {
    n_clust <- n_clust[order(n_clust)]
    if(find_h){
      clusters <- data.frame(matrix(
        nrow = length(cur.tree$labels),
        ncol = length(n_clust) + 1,
        dimnames = list(cur.tree$labels,
                        c("name", paste0("k_", n_clust)))))
      clusters$name <- cur.tree$labels
      for(cur_n in n_clust){
        if(length(n_clust) < 10){
          message("Determining the cut height to reach ", cur_n, " groups...")
        }
        k <- 0
        h1 <- h_max
        h0 <- h_min
        h <- h0 + (h1 - h0) / 2
        # Algorithm to quickly find the height of cut corresponding to the
        # requested number of clusters
        iter <- 0
        max_iter <- 500
        while(k != cur_n &  h1 != h0 & iter < max_iter){
          h <- h0 + (h1 - h0) / 2
          cls <- stats::cutree(cur.tree, h = h)
          k <- max(cls)
          if(k < cur_n) {
            h1 <- h
          } else if (k > cur_n) {
            h0 <- h
          }
          iter <- iter + 1
        }
        if(length(n_clust) < 10){
          message(paste0("--> ", h))
        }
        if(k != cur_n) {
          warning(paste0("The requested number of cluster could not be found ",
                         "for k = ",
                         cur_n, 
                         ". Closest number found: ", 
                         k))
        }
        clusters[, paste0("k_", cur_n)] <- as.character(cls)
        output_cut_height <- c(output_cut_height, h)
        
        output_n_clust <- c(output_n_clust, k)
      }
      names(output_cut_height) <- paste0("k_", n_clust)
    } else {
      cls <- stats::cutree(cur.tree, k = n_clust)
      if(length(n_clust) == 1){
        clusters <- data.frame(names(cls),
                               cluster = cls)
      } else {
        clusters <- data.frame(rownames(cls),
                               cluster = cls)
      }

      names(clusters) <- c("name", paste0("k_", n_clust))
      output_cut_height <- "unknown"
      output_n_clust <- sapply(n_clust,
                               function(k, cl){
                                 length(unique(cl[, paste0("k_", k)]))
                               }, cl = clusters)
    }
    
  } else if(!is.null(cut_height)) {
    cut_height <- cut_height[order(cut_height, decreasing = TRUE)]
    cls <- stats::cutree(cur.tree,
                              h = cut_height)
    if(length(cut_height) == 1) {
      clusters <- data.frame(site = names(cls),
                             cluster = as.character(cls))
    } else {
      clusters <- data.frame(site = rownames(cls),
                             cluster = cls)
    }
    colnames(clusters) <- c("site", paste0("h_",  cut_height))
    output_n_clust <- sapply(cut_height,
                             function(h, cl){
                               length(unique(cl[, paste0("h_", h)]))
                             }, cl = clusters)
    names(output_n_clust) <- paste0("h_", cut_height)
    output_cut_height <- cut_height
  }
  
  # Note: show_hierarchy is stored in args but doesn't change behavior 
  # significantly for hierarchical clustering since cutree returns integers.
  # The hierarchy is already preserved through multiple partition columns.
  clusters <- knbclu(clusters, reorder = TRUE)
  
  if(inherits(tree, "bioregion.clusters")) {
    cur.tree$args$cut_height <- cut_height
    tree$clusters <- clusters
    tree$algorithm$output_n_clust <- output_n_clust
    tree$algorithm$output_cut_height <- output_cut_height
    
    tree$cluster_info <- data.frame(
      partition_name = names(tree$clusters)[2:length(tree$clusters),
                                            drop = FALSE],
      n_clust = output_n_clust)
    
    if(!is.null(n_clust)) {
      tree$cluster_info$requested_n_clust <- n_clust
      if(find_h) {
        tree$cluster_info$output_cut_height <- output_cut_height
      }
    } else if(!is.null(cut_height)) {
      tree$cluster_info$requested_cut_height <- cut_height
    }
    
    # Update hierarchical status based on number of partitions
    tree$inputs$hierarchical <- ifelse(ncol(tree$clusters) > 2,
                                       TRUE,
                                       FALSE)
    
    return(tree)
  } else if (inherits(tree, "hclust")){
    return(clusters)
  }
}
