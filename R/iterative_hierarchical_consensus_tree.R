iterative_consensus_tree <- function(
    dissim, 
    sites = unique(c(dissim[, 1],
                            dissim[, 2])), 
    index = colnames(dissim)[3],
    method = "average",
    depth = 1, 
    previous_height = Inf, 
    verbose = TRUE,
    n_runs = 100,
    max_remaining_size = length(sites),
    monotonicity_direction = c("top-down", "bottom-up")) {

  if(inherits(dissim, "bioregion.pairwise.metric") |
     inherits(dissim, "data.frame")) {
    dissim[, 3] <- dissim[, index]
    dissim <- stats::as.dist(
      net_to_mat(dissim[, 1:3], weight = TRUE, squared = TRUE, symmetrical = TRUE)
    )
  } 
  
  monotonicity_direction <- match.arg(monotonicity_direction)
  
  # Progress info because the algorithm is quite long to run
  if (verbose && interactive()) {
    if(max_remaining_size > 2) {
      cat(sprintf(
        "\rProcessing height: %.3f | Current branch size: %d | Max cluster size remaining: %d     ", 
        previous_height,
        length(sites),
        max(unlist(max_remaining_size))
      ))
    } else if(max_remaining_size == 2) {
      cat(sprintf(
        "\rProcessing height: %.3f | Current branch size: %d | All clusters successfully processed     ", 
        previous_height,
        length(sites)
      ))
    }
    
    utils::flush.console()  # Ensures the message is displayed in real-time
  }
  
  # Base case: If the subcluster is a single-site cluster, return it as a leaf node
  if (length(sites) == 1) {
    return(list(name = sites, height = 0, depth = depth, children = NULL))
  }
  
  # Step 1: subset dissimilarity matrix to current cluster only 
  dist_matrix <- subset_dist(dissim,
                             sites)
  
  # current_net <- net[net[, 1] %in% sites, ]
  # comat <- net_to_mat(current_net, weight = FALSE)
  # 
  # # Step 2: dissimilarity matrix
  # dissim <- dissimilarity(comat, metric = index)
  # dist_matrix <- stats::as.dist(
  #   net_to_mat(dissim, weight = TRUE, squared = TRUE, symmetrical = TRUE)
  # )
  
  # Step 3: find the majority vote for a binary split among many trees
  if(nrow(dist_matrix) > 2) { # Only if we have more than 2 sites
    subclusters <- stable_binary_split(dist_matrix, 
                                       method = method,
                                       n_runs = n_runs,
                                       binsplit = "tree")
  } else { # if 2 sites only we already have the subclusters
    subclusters <- list(attr(dist_matrix, "Labels")[1],
                        attr(dist_matrix, "Labels")[2])
  }
  
  # Ensure the smallest subcluster is the first (same rules as for hclust)
  if (length(subclusters[[1]]) > length(subclusters[[2]])) {
    subclusters <- list(subclusters[[2]], subclusters[[1]])
  }
  
  # Calculate height based on the specified method
  cluster1_sites <- subclusters[[1]]
  cluster2_sites <- subclusters[[2]]
  pairwise_distances <- as.matrix(dist_matrix)[cluster1_sites, cluster2_sites]
  
  # Calculate centroids for ward.D and ward.D2
  centroid1 <- colMeans(as.matrix(dist_matrix)[cluster1_sites, , drop = FALSE])
  centroid2 <- colMeans(as.matrix(dist_matrix)[cluster2_sites, , drop = FALSE])
  
  # Compute squared Euclidean distance between centroids
  centroid_distance <- sum((centroid1 - centroid2)^2)

  # Compute height based on the selected method
  calculated_height <- switch(
    method,
    "single" = min(pairwise_distances),
    "complete" = max(pairwise_distances),
    "average" = mean(pairwise_distances),
    "mcquitty" = mean(pairwise_distances),  # Note: This is acceptable for initial height
    "ward.D" = (length(cluster1_sites) * length(cluster2_sites)) / 
      (length(cluster1_sites) + length(cluster2_sites)) * sqrt(centroid_distance),
    "ward.D2" = (length(cluster1_sites) * length(cluster2_sites)) /  
      (length(cluster1_sites) + length(cluster2_sites)) * centroid_distance,
    "centroid" = centroid_distance, 
    "median" = 0.5 * centroid_distance,  
    stop("method argument is not valid")
  )
  
  # Step 4: height constraint (monotonic constraint of trees)
  if (monotonicity_direction == "top-down") {
    height <- min(calculated_height, previous_height)
  } else if (monotonicity_direction == "bottom-up") {
    # Height will be adjusted after processing subclusters
    height <- calculated_height
  }
  if (height < 0) height <- 0
  
  # If we are at the iteration where we work on the max cluster size,
  # then we update max_remaining_size based on the largest subcluster size
  # (always in position 2)
  if (length(sites) == max_remaining_size) {
    max_remaining_size <- length(subclusters[[2]])
  }
  
  # # Process each subcluster and collect their structures
  children <- list()
  for (i in seq_along(subclusters)) {
    subcluster <- subclusters[[i]]
    child_structure <- iterative_consensus_tree(
      dissim = dissim,
      sites = subcluster,
      # index = index,
      method = method,
      depth = depth + 1,
      previous_height = height,
      verbose = verbose,
      max_remaining_size = max_remaining_size,
      monotonicity_direction = monotonicity_direction)
    children[[i]] <- child_structure
  }
  
  # After processing subclusters, adjust height if monotonicity is bottom-up
  if (monotonicity_direction == "bottom-up") {
    max_child_height <- max(sapply(children, function(x) x$height))
    if (max_child_height > calculated_height) {
      height <- max_child_height
    } else {
      height <- calculated_height
    }
  }
  
  # Return the tree structure for the current split
  return(list(
    name = paste(sites, collapse = ","),
    height = height,
    depth = depth,
    children = children
  ))
}

subset_dist <- function(d, indices) {
  d_mat <- as.matrix(d)
  d_subset <- d_mat[indices, indices]
  as.dist(d_subset)
}


# Function to find the majority vote for a binary split among many trees
stable_binary_split <- function(dist.obj,
                                method = "average", 
                                n_runs = 100,
                                binsplit = "tree") {
  # Step 1: randomize distance matrices and generate trees
  randomtrees <- list()
  # no need to make more permutations than the max number of possible 
  # permutations so we set the max value to the lowest of n_runs or 
  # factorial(n_runs)
  for (run in 1:min(n_runs, factorial(length(dist.obj)))) {
    randomtrees[[run]] <- list()
    randomtrees[[run]]$dist.matrix <- .randomizeDistance(dist.obj)
    randomtrees[[run]]$hierartree <-  
      fastcluster::hclust(randomtrees[[run]]$dist.matrix, 
                          method = method)
  }
  
  # Step 2: cut trees at 2 clusters for each run
  trees <- lapply(randomtrees, function(trial) trial$hierartree)
  clusts <- lapply(trees, function(tree) {
    suppressMessages(cut_tree(
      stats::as.hclust(tree), n_clust = 2, find_h = FALSE))
  })
  
  # Step 3: ensure we have the same order of sites in each partition
  fixed_order <- sort(clusts[[1]]$ID) 
  df_list_ordered <- lapply(clusts, function(df) {
    df[match(fixed_order, df$ID), ]
  })
  
  # Step 4: make a data frame with all partitions to use the same functions as
  # in compare_partitions()
  partitions <- data.frame(ID = fixed_order)
  partitions <- cbind(partitions, lapply(df_list_ordered, `[[`, 2))
  colnames(partitions)[-1] <- paste0("run", seq_along(clusts))
  
  # Step 5: Calculate pairwise memberships for clustering consistency
  site_pw_memberships <- get_pairwise_membership(partitions[, -1])
  site_names <- do.call(rbind, strsplit(rownames(site_pw_memberships), "_"))
  
  # Step 6: Compute pairwise dissimilarity based on memberships
  # First as a network format as it is easier to handle
  pairwise_proportion <- data.frame(
    Site1 = fixed_order[as.numeric(site_names[, 1])],
    Site2 = fixed_order[as.numeric(site_names[, 2])],
    Weight = rowSums(!site_pw_memberships) / ncol(site_pw_memberships)
  ) # Note with use !site_pw_memberships, because we need dissimilarity here,
  # not similarity
  
  # Then we convert to dist object
  dist_pw_prop <- stats::as.dist(
    net_to_mat(pairwise_proportion, weight = TRUE, squared = TRUE, symmetrical = TRUE)
  )
  
  # Step 7: Split into two clusters form the dist matrix
  if (binsplit == "tree") {
    pw_tree <- fastcluster::hclust(dist_pw_prop, method = method)
    groups <- stats::cutree(pw_tree, k = 2)
  } else if (binsplit == "pam") {
    groups <- cluster::pam(dist_pw_prop, k = 2)$clustering
  }
  groups <- data.frame(Site = names(groups), cluster = groups)

  # tree_eval_cur <- tree_eval(pw_tree,
  #                            dist_pw_prop)$cophcor
  # cat ("Tree eval :", tree_eval_cur, "\n")
  # # if(tree_eval_cur < .5) {
  #   # cat ("Tree eval < .5 for tree ", paste0(sites))
  #   randomtrees2 <- list()
  #   for (run in 1:10) {
  #     randomtrees2[[run]] <- list()
  #     randomtrees2[[run]]$dist.matrix <- .randomizeDistance(dist_pw_prop)
  #     randomtrees2[[run]]$hierartree <-  
  #       fastcluster::hclust(randomtrees2[[run]]$dist.matrix, 
  #                           method = method)
  #     randomtrees2[[run]]$cophcor <- 
  #       tree_eval(randomtrees2[[run]]$hierartree,
  #                 randomtrees2[[run]]$dist.matrix)$cophcor
  #   }
  #   coph.coeffs <- sapply(1:10, function(x)
  #   {
  #     randomtrees2[[run]]$cophcor
  #   })
  #   
  #   message(paste0(" -- range of cophenetic correlation coefficients among
  #                    trials: ", round(min(coph.coeffs), 4),
  #                  " - ", round(max(coph.coeffs), 4)))
  #   
  #   # There might be multiple trees with the highest cophenetic correlation
  #   # coefficient, so we arbitrarily take the first one
  #   best.run <- which(coph.coeffs == max(coph.coeffs))[1]
  #   
  #   pw_tree <- randomtrees2[[best.run]]$hierartree
  # # }
  # It is not necessary to randomize the tree here 
  # because it does not impact the tree topology


  # png(paste0("temp/", as.numeric(Sys.time()), ".png"))
  # plot(pw_tree)
  # dev.off()
  
  return(split(groups$Site, groups$cluster))
}


reconstruct_hclust <- function(tree) {
  # Step 1: Assign negative IDs to leaf nodes based on alphabetical order
  get_leaves <- function(node) {
    if (is.null(node$children)) {
      return(node$name)
    } else {
      return(unlist(lapply(node$children, get_leaves)))
    }
  }
  
  leaf_labels <- sort(unique(get_leaves(tree)))
  n_leaves <- length(leaf_labels)
  
  # Map leaf labels to negative IDs based on alphabetical order
  leaf_ids <- stats::setNames(-seq_along(leaf_labels), leaf_labels)
  
  # Step 2: Traverse the tree and collect internal nodes
  internal_nodes <- list()
  node_counter <- 0
  
  collect_internal_nodes <- function(node) {
    if (is.null(node$children)) {
      # Leaf node
      node_id <- leaf_ids[node$name]
      return(node_id)
    } else {
      # Internal node
      left_id <- collect_internal_nodes(node$children[[1]])
      right_id <- collect_internal_nodes(node$children[[2]])
      
      node_counter <<- node_counter + 1
      temp_id <- paste0("n", node_counter)
      internal_nodes[[temp_id]] <<- list(
        left = left_id,
        right = right_id,
        height = node$height
      )
      return(temp_id)
    }
  }
  
  # Collect all internal nodes
  root_id <- collect_internal_nodes(tree)
  
  # Now, we have internal_nodes as a list with temporary IDs, and their left, right, height
  # We need to assign actual IDs to internal nodes, in order of increasing height
  
  # First, create a data frame of internal nodes
  internal_nodes_df <- data.frame(
    temp_id = names(internal_nodes),
    left = sapply(internal_nodes, function(x) x$left),
    right = sapply(internal_nodes, function(x) x$right),
    height = sapply(internal_nodes, function(x) x$height),
    stringsAsFactors = FALSE
  )
  
  # Sort internal nodes by height (increasing order)
  internal_nodes_df <- internal_nodes_df[order(internal_nodes_df$height, decreasing = FALSE), ]
  
  # Assign IDs to internal nodes starting from 1
  internal_node_ids <- seq_len(nrow(internal_nodes_df))
  names(internal_node_ids) <- internal_nodes_df$temp_id
  
  # Function to map IDs
  map_id <- function(x) {
    if (x %in% names(internal_node_ids)) {
      return(internal_node_ids[x])
    } else {
      return(as.integer(x))  # Should be negative IDs for leaves
    }
  }
  
  # Now, update left and right IDs in internal_nodes_df using the mapping
  internal_nodes_df$left <- sapply(internal_nodes_df$left, map_id)
  internal_nodes_df$right <- sapply(internal_nodes_df$right, map_id)
  
  # Build the merge matrix and height vector
  merge <- as.matrix(internal_nodes_df[, c("left", "right")])
  height <- internal_nodes_df$height
  
  # Compute order vector via depth-first traversal
  order <- integer(0)
  compute_order <- function(node) {
    if (is.null(node$children)) {
      leaf_index <- which(leaf_labels == node$name)
      order <<- c(order, leaf_index)
    } else {
      compute_order(node$children[[1]])
      compute_order(node$children[[2]])
    }
  }
  compute_order(tree)
  
  # Assemble the hclust object
  hclust_obj <- list(
    merge = merge,
    height = height,
    order = order,
    labels = leaf_labels,
    method = "Iterative Hierarchical Consensus Tree"
  )
  class(hclust_obj) <- "hclust"
  return(hclust_obj)
}

# Tests
# hc$height
# 
# b <- reconstruct_hclust(a)
# 
# plot(b)
# 
# cophenetic_debug(b)
# cutree(b, k=5)
# cutree(b, h=0.1)
# 
# 
# cophenetic_debug <- function(x) {
#   x <- as.hclust(x)
#   nobs <- length(x$order)
#   ilist <- vector("list", length = nobs)
#   out <- matrix(0, nrow = nobs, ncol = nobs)
#   for (i in 1:(nobs - 1)) {
#     inds <- x$merge[i, ]
#     ids1 <- if (inds[1L] < 0L) -inds[1L] else ilist[[inds[1L]]]
#     ids2 <- if (inds[2L] < 0L) -inds[2L] else ilist[[inds[2L]]]
#     if (is.null(ids1) || is.null(ids2)) {
#       print(paste("NULL encountered at row:", i))
#       print(inds)
#       stop("NULL values in cophenetic calculation")
#     }
#     ilist[[i]] <- c(ids1, ids2)
#     out[cbind(rep.int(ids1, rep.int(length(ids2), length(ids1))),
#               rep.int(ids2, length(ids1)))] <- x$height[i]
#   }
#   rownames(out) <- x$labels
#   as.dist(out + t(out))
# }
# 

# Previous working versions versions in case we need to go back to
# the code at some point:
# iterative_consensus_tree <- function(net, 
#                                      sites = unique(net[, 1]), 
#                                      index = NULL,
#                                      method = "average",
#                                      depth = 1, 
#                                      tree_structure = list(), 
#                                      previous_height = Inf, 
#                                      verbose = TRUE,
#                                      n_runs = 100,
#                                      max_remaining_size = length(sites)) {
# 
# 
#   # Progress info because the algorithm is quite long to run
#   if (verbose && interactive()) {
#     if(max_remaining_size > 2) {
#       cat(sprintf(
#         "\rProcessing height: %.3f | Current branch size: %d | Max cluster size remaining: %d     ", 
#         previous_height,
#         length(sites),
#         max(unlist(max_remaining_size))
#       ))
#     } else if(max_remaining_size == 2) {
#       cat(sprintf(
#         "\rProcessing height: %.3f | Current branch size: %d | All clusters successfully processed     ", 
#         previous_height,
#         length(sites)
#       ))
#     }
# 
#     utils::flush.console()  # Ensures the message is displayed in real-time
#   }
# 
#   
#   # Step 1: species composition matrix for the current cluster
#   current_net <- net[net[, 1] %in% sites, ]
#   comat <- net_to_mat(current_net, weight = FALSE)
#   
#   # Step 2: dissimilarity matrix
#   # TODO: implement the possibility of using other indices than the ones we have 
#   # in the package already, e.g. for phylobeta diversity
#   # TODO: implement correctly the formulas for abc and ABC
#   dissim <- dissimilarity(comat, metric = index)
#   dist_matrix <- stats::as.dist(
#     net_to_mat(dissim, weight = TRUE, squared = TRUE, symmetrical = TRUE)
#   )
#   
#   # Step 3: find the majority vote for a binary split among many trees
#   if(nrow(comat) > 2) { # Only if we have more than 2 sites
#     subclusters <- stable_binary_split(dist_matrix, 
#                                        method = method,
#                                        n_runs = n_runs,
#                                        binsplit = "tree")
#   } else { # if 2 sites only we already have the subclusters
#     subclusters <- list(rownames(comat)[1],
#                         rownames(comat)[2])
#   }
# 
#   
#   # Ensure the smallest subcluster is the first (same rules as for hclust)
#   if (length(subclusters[[1]]) > length(subclusters[[2]])) {
#     subclusters <- list(subclusters[[2]], subclusters[[1]])
#   }
#   
#   # Calculate height based on the specified method
#   cluster1_sites <- subclusters[[1]]
#   cluster2_sites <- subclusters[[2]]
#   pairwise_distances <- as.matrix(dist_matrix)[cluster1_sites, cluster2_sites]
#   
#   # Calculate centroids for ward.D and ward.D2
#   centroid1 <- colMeans(comat[cluster1_sites, , drop = FALSE])
#   centroid2 <- colMeans(comat[cluster2_sites, , drop = FALSE])
#   centroid_distance <- sum((centroid1 - centroid2)^2)
#   
#   # Compute height based on the selected method
#   # TODO: check if ward calculation align with calculations in hclust &
#   # fastclust
#   calculated_height <- switch(
#     method,
#     "single" = min(pairwise_distances),
#     "complete" = max(pairwise_distances),
#     "average" = mean(pairwise_distances),
#     "mcquitty" = mean(c(
#       mean(as.matrix(dist_matrix)[cluster1_sites, ]),
#       mean(as.matrix(dist_matrix)[cluster2_sites, ])
#     )),
#     "ward.D" = (length(cluster1_sites) * length(cluster2_sites)) / 
#       (length(cluster1_sites) + length(cluster2_sites)) * centroid_distance,
#     "ward.D2" = 0.5 * (length(cluster1_sites) * length(cluster2_sites)) / 
#       (length(cluster1_sites) + length(cluster2_sites)) * centroid_distance,
#     "centroid" = centroid_distance,
#     "median" = stats::median(pairwise_distances),
#     stop("method argument is not valid")
#   )
#   
#   # Step 4: height constraint (monotonic constraint of trees)
#   height <- min(calculated_height, previous_height)
#   if (height < 0) height <- 0
#   
#   # update tree_structure for the current split
#   tree_structure <- c(tree_structure, list(
#     list(
#       cluster = sites,
#       subclusters = subclusters,
#       height = height,
#       depth = depth
#     )
#   ))
#   
#   # If we are the iteration where we work on the max cluster size,
#   # then we update max_remaining_size based on the largest subcluster size
#   # (always in position 2)
#   if (length(sites) == max_remaining_size) {
#     max_remaining_size <- length(subclusters[[2]])
#   }
#   
#   # Next we process each subcluster
#   for (subcluster in subclusters) {
#     if (length(subcluster) == 1) {
#       # If the subcluster is a single-site cluster, add it as a leaf node
#       tree_structure <- c(tree_structure, list(list(leaf = subcluster)))
#     } else {
#       # Else if it is a normal cluster, recursively process it, passing the 
#       # current height as the previous height
#       tree_structure <- iterative_consensus_tree(
#         net = net, 
#         sites = subcluster, 
#         index = index,
#         method = method,  
#         depth = depth + 1, 
#         tree_structure = tree_structure, 
#         previous_height = height, 
#         verbose = verbose,
#         max_remaining_size = max_remaining_size)
#     }
#   }
# 
# 
# 
# 
#   # It works!! :)
#   return(tree_structure)
# }


# Function to make hclust object based on the output list we prepared above
# Much harder to make than it looks...
# reconstruct_hclust <- function(tree_structure) {
#   # Extract all labels and sort them alphabetically
#   labels <- sort(tree_structure[[1]]$cluster)
#   n <- length(labels)
#   
#   # We need to make sure that the labels do not contain our concatenation
#   # character or it may mess up the results
#   if (any(grepl("|", labels, fixed = TRUE))) {
#     stop(paste0("Error: Labels contain the separator character '", "|",
#                 "'. Please remove this character from labels"))
#   }
#   
#   
#   # Assign negative indices to leaf nodes based on alphabetical order
#   label_to_index <- stats::setNames(-seq_along(labels), labels)
#   
#   # Initialize variables
#   cluster_to_id <- list()  # Maps cluster keys to node IDs
#   merge_list <- list() # Merge object of hclust: matrix indicating how to 
#   # assemble nodes in the tree
#   height_list <- numeric() # List of heights for the hclust object, must be
#   # increasing
#   merge_index <- 1  # Index of the merge operation, corresponds to internal 
#   # node ID
#   
#   
#   # Create a mapping from cluster keys to nodes
#   cluster_to_node <- list()
#   for (node in tree_structure) {
#     if (!is.null(node$cluster)) {
#       key <- paste(sort(node$cluster), collapse = "|")
#       cluster_to_node[[key]] <- node
#     } else if (!is.null(node$leaf)) {
#       key <- node$leaf
#       cluster_to_node[[key]] <- node
#     }
#   }
#   
# 
#   # Assign IDs to leaf nodes
#   for (label in labels) {
#     cluster_key <- label
#     cluster_to_id[[cluster_key]] <- label_to_index[[label]]  # Negative ID
#   }
#   
#   # Collect all unique heights and sort them
#   node_heights <- sapply(tree_structure, function(node) {
#     if (!is.null(node$height)) {
#       return(node$height)
#     } else {
#       return(Inf)   # Assign infinite height to leaf nodes just for the next 
#       # loop
#     }
#   })
#   unique_heights <- sort(unique(node_heights))
#   
#   # Process nodes in order of increasing height
#   for (height in unique_heights) {
#     # Skip processing for leaf nodes (height == Inf)
#     if (height == Inf) {
#       next
#     }
#     
#     # Get all nodes at this height
#     nodes_at_height <- lapply(tree_structure, function(node) {
#       node_height <- ifelse(is.null(node$height), Inf, node$height)
#       if (node_height == height) {
#         return(node)
#       } else {
#         return(NULL)
#       }
#     })
#     nodes_at_height <- nodes_at_height[!sapply(nodes_at_height, is.null)]
#     
#     # Process nodes
#     pending_nodes <- nodes_at_height
#     while (length(pending_nodes) > 0) {
#       nodes_to_retry <- list()
#       for (node in pending_nodes) {
#         subclusters <- node$subclusters
#         left_subcluster <- subclusters[[1]]
#         right_subcluster <- subclusters[[2]]
#         
#         # Function to get ID of a subcluster
#         get_subcluster_id <- function(subcluster) {
#           if (length(subcluster) == 1 && subcluster %in% labels) {
#             # Leaf node
#             return(label_to_index[[subcluster]])  # Negative ID
#           } else {
#             key <- paste(sort(subcluster), collapse = "|")
#             id <- cluster_to_id[[key]]
#             return(id)
#           }
#         }
#         
#         left_id <- get_subcluster_id(left_subcluster)
#         right_id <- get_subcluster_id(right_subcluster)
#         
#         if (!is.null(left_id) && !is.null(right_id)) {
#           # Both subclusters have IDs, process the node
#           # Record the merge
#           merge_list[[merge_index]] <- c(left_id, right_id)
#           height_list[merge_index] <- height
#           
#           # Assign an ID to the current node (positive integer)
#           current_id <- merge_index
#           cluster_key <- paste(sort(node$cluster), collapse = "|")
#           cluster_to_id[[cluster_key]] <- current_id
#           
#           # Diagnostics to solve the issues we had
#           # cat("Assigned ID", current_id, "to cluster", cluster_key, "\n")
#           # cat("Merge", merge_index, ":", "Merging IDs", left_id, "and",
#           # right_id, "at height", height, "\n")
#           
#           merge_index <- merge_index + 1
#         } else {
#           # Subcluster IDs not yet assigned, defer processing
#           nodes_to_retry[[length(nodes_to_retry) + 1]] <- node
#         }
#       }
#       
#       if (length(nodes_to_retry) == length(pending_nodes)) {
#         stop("Cannot resolve dependencies among internal nodes at height ",
#              height,
#              "\nPlease contact us with your objects so that we can resolve",
#              " the bug.")
#       }
#       
#       pending_nodes <- nodes_to_retry
#     }
#   }
#   
#   # Convert the merge list to a matrix
#   merge_matrix <- do.call(rbind, merge_list)
#   
#   # Ensure that the height vector is in increasing order
#   if (any(diff(height_list) < 0)) {
#     stop("Heights must be in increasing order.")
#   }
#   
#   # Function to perform a depth-first traversal to get the order
#   order <- numeric()
#   traverse_order <- function(node_id) {
#     if (node_id < 0) {
#       # Leaf node
#       idx <- which(label_to_index == node_id)
#       order[length(order) + 1] <<- idx
#     } else {
#       # Internal node
#       merge_idx <- node_id
#       left_id <- merge_matrix[merge_idx, 1]
#       right_id <- merge_matrix[merge_idx, 2]
#       traverse_order(left_id)
#       traverse_order(right_id)
#     }
#   }
#   
#   # Start traversal from the last internal node
#   root_node_id <- merge_index - 1
#   traverse_order(root_node_id)
#   
#   # Construct the hclust object
#   hclust_obj <- list(
#     merge = merge_matrix,
#     height = height_list,
#     order = order,
#     labels = labels,
#     method = "Iterative Hierarchical Consensus Tree"
#   )
#   class(hclust_obj) <- "hclust"
#   
#   # FINALLY!
#   return(hclust_obj)
# }