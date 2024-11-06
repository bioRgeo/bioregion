iterative_consensus_tree <- function(net, 
                                     sites = unique(net[, 1]), 
                                     index = NULL,
                                     method = "average",
                                     depth = 1, 
                                     tree_structure = list(), 
                                     previous_height = Inf, 
                                     verbose = TRUE,
                                     n_runs = 100,
                                     max_remaining_size = length(sites)) {


  # Progress info because the algorithm is quite long to run
  if (verbose && interactive()) {
    cat(sprintf(
      "\rProcessing height: %.3f | Current branch size: %d | Max cluster size remaining: %d     ", 
      previous_height,
      length(sites),
      max(unlist(max_remaining_size))
    ))
    utils::flush.console()  # Ensures the message is displayed in real-time
  }

  
  # Step 1: species composition matrix for the current cluster
  current_net <- net[net[, 1] %in% sites, ]
  comat <- net_to_mat(current_net, weight = FALSE)
  
  # Step 2: dissimilarity matrix
  # TODO: implement the possibility of using other indices than the ones we have 
  # in the package already, e.g. for phylobeta diversity
  # TODO: implement correctly the formulas for abc and ABC
  dissim <- dissimilarity(comat, metric = index)
  dist_matrix <- stats::as.dist(
    net_to_mat(dissim, weight = TRUE, squared = TRUE, symmetrical = TRUE)
  )
  
  # Step 3: find the majority vote for a binary split among many trees
  subclusters <- stable_binary_split(dist_matrix, 
                                     method = method,
                                     n_runs = n_runs)
  
  # Ensure the smallest subcluster is the first (same rules as for hclust)
  if (length(subclusters[[1]]) > length(subclusters[[2]])) {
    subclusters <- list(subclusters[[2]], subclusters[[1]])
  }
  
  # Calculate height based on the specified method
  cluster1_sites <- subclusters[[1]]
  cluster2_sites <- subclusters[[2]]
  pairwise_distances <- as.matrix(dist_matrix)[cluster1_sites, cluster2_sites]
  
  # Calculate centroids for ward.D and ward.D2
  centroid1 <- colMeans(comat[cluster1_sites, , drop = FALSE])
  centroid2 <- colMeans(comat[cluster2_sites, , drop = FALSE])
  centroid_distance <- sum((centroid1 - centroid2)^2)
  
  # Compute height based on the selected method
  # TODO: check if ward calculation align with calculations in hclust &
  # fastclust
  calculated_height <- switch(
    method,
    "single" = min(pairwise_distances),
    "complete" = max(pairwise_distances),
    "average" = mean(pairwise_distances),
    "mcquitty" = mean(c(
      mean(as.matrix(dist_matrix)[cluster1_sites, ]),
      mean(as.matrix(dist_matrix)[cluster2_sites, ])
    )),
    "ward.D" = (length(cluster1_sites) * length(cluster2_sites)) / 
      (length(cluster1_sites) + length(cluster2_sites)) * centroid_distance,
    "ward.D2" = 0.5 * (length(cluster1_sites) * length(cluster2_sites)) / 
      (length(cluster1_sites) + length(cluster2_sites)) * centroid_distance,
    "centroid" = centroid_distance,
    "median" = stats::median(pairwise_distances),
    stop("method argument is not valid")
  )
  
  # Step 4: height constraint (monotonic constraint of trees)
  height <- min(calculated_height, previous_height)
  if (height < 0) height <- 0
  
  # update tree_structure for the current split
  tree_structure <- c(tree_structure, list(
    list(
      cluster = sites,
      subclusters = subclusters,
      height = height,
      depth = depth
    )
  ))
  
  # If we are the iteration where we work on the max cluster size,
  # then we update max_remaining_size based on the largest subcluster size
  # (always in position 2)
  if (length(sites) == max_remaining_size) {
    max_remaining_size <- length(subclusters[[2]])
  }
  
  # Next we process each subcluster
  for (subcluster in subclusters) {
    if (length(subcluster) == 1) {
      # If the subcluster is a single-site cluster, add it as a leaf node
      tree_structure <- c(tree_structure, list(list(leaf = subcluster)))
    } else {
      # Else if it is a normal cluster, recursively process it, passing the 
      # current height as the previous height
      tree_structure <- iterative_consensus_tree(
        net = net, 
        sites = subcluster, 
        index = index,
        method = method,  
        depth = depth + 1, 
        tree_structure = tree_structure, 
        previous_height = height, 
        verbose = verbose,
        max_remaining_size = max_remaining_size)
    }
  }




  # It works!! :)
  return(tree_structure)
}


# Function to find the majority vote for a binary split among many trees
stable_binary_split <- function(dist.obj,
                                method = "average", 
                                n_runs = 100) {
  
  # Step 1: randomize distance matrices and generate trees
  randomtrees <- list()
  for (run in 1:n_runs) {
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
  
  # Step 7: Create hierarchical tree from the pairwise dissimilarity matrix
  pw_tree <- fastcluster::hclust(dist_pw_prop, method = method)
  # Note we do not randomize here, but the tree topology is very conspicuous 
  # and so we probably do not need to randomize it
  groups <- stats::cutree(pw_tree, k = 2)
  groups <- data.frame(Site = names(groups), cluster = groups)
  
  return(split(groups$Site, groups$cluster))
}


# Function to make hclust object based on the output list we prepared above
# Much harder to make than it looks...
reconstruct_hclust <- function(tree_structure) {
  # Extract all labels and sort them alphabetically
  labels <- sort(tree_structure[[1]]$cluster)
  n <- length(labels)
  
  # We need to make sure that the labels do not contain our concatenation
  # character or it may mess up the results
  if (any(grepl("|", labels, fixed = TRUE))) {
    stop(paste0("Error: Labels contain the separator character '", "|",
                "'. Please remove this character from labels"))
  }
  
  
  # Assign negative indices to leaf nodes based on alphabetical order
  label_to_index <- stats::setNames(-seq_along(labels), labels)
  
  # Initialize variables
  cluster_to_id <- list()  # Maps cluster keys to node IDs
  merge_list <- list() # Merge object of hclust: matrix indicating how to 
  # assemble nodes in the tree
  height_list <- numeric() # List of heights for the hclust object, must be
  # increasing
  merge_index <- 1  # Index of the merge operation, corresponds to internal 
  # node ID
  
  
  # Create a mapping from cluster keys to nodes
  cluster_to_node <- list()
  for (node in tree_structure) {
    if (!is.null(node$cluster)) {
      key <- paste(sort(node$cluster), collapse = "|")
      cluster_to_node[[key]] <- node
    } else if (!is.null(node$leaf)) {
      key <- node$leaf
      cluster_to_node[[key]] <- node
    }
  }
  

  # Assign IDs to leaf nodes
  for (label in labels) {
    cluster_key <- label
    cluster_to_id[[cluster_key]] <- label_to_index[[label]]  # Negative ID
  }
  
  # Collect all unique heights and sort them
  node_heights <- sapply(tree_structure, function(node) {
    if (!is.null(node$height)) {
      return(node$height)
    } else {
      return(Inf)   # Assign infinite height to leaf nodes just for the next 
      # loop
    }
  })
  unique_heights <- sort(unique(node_heights))
  
  # Process nodes in order of increasing height
  for (height in unique_heights) {
    # Skip processing for leaf nodes (height == Inf)
    if (height == Inf) {
      next
    }
    
    # Get all nodes at this height
    nodes_at_height <- lapply(tree_structure, function(node) {
      node_height <- ifelse(is.null(node$height), Inf, node$height)
      if (node_height == height) {
        return(node)
      } else {
        return(NULL)
      }
    })
    nodes_at_height <- nodes_at_height[!sapply(nodes_at_height, is.null)]
    
    # Process nodes
    pending_nodes <- nodes_at_height
    while (length(pending_nodes) > 0) {
      nodes_to_retry <- list()
      for (node in pending_nodes) {
        subclusters <- node$subclusters
        left_subcluster <- subclusters[[1]]
        right_subcluster <- subclusters[[2]]
        
        # Function to get ID of a subcluster
        get_subcluster_id <- function(subcluster) {
          if (length(subcluster) == 1 && subcluster %in% labels) {
            # Leaf node
            return(label_to_index[[subcluster]])  # Negative ID
          } else {
            key <- paste(sort(subcluster), collapse = "|")
            id <- cluster_to_id[[key]]
            return(id)
          }
        }
        
        left_id <- get_subcluster_id(left_subcluster)
        right_id <- get_subcluster_id(right_subcluster)
        
        if (!is.null(left_id) && !is.null(right_id)) {
          # Both subclusters have IDs, process the node
          # Record the merge
          merge_list[[merge_index]] <- c(left_id, right_id)
          height_list[merge_index] <- height
          
          # Assign an ID to the current node (positive integer)
          current_id <- merge_index
          cluster_key <- paste(sort(node$cluster), collapse = "|")
          cluster_to_id[[cluster_key]] <- current_id
          
          # Diagnostics to solve the issues we had
          # cat("Assigned ID", current_id, "to cluster", cluster_key, "\n")
          # cat("Merge", merge_index, ":", "Merging IDs", left_id, "and",
          # right_id, "at height", height, "\n")
          
          merge_index <- merge_index + 1
        } else {
          # Subcluster IDs not yet assigned, defer processing
          nodes_to_retry[[length(nodes_to_retry) + 1]] <- node
        }
      }
      
      if (length(nodes_to_retry) == length(pending_nodes)) {
        stop("Cannot resolve dependencies among internal nodes at height ",
             height,
             "\nPlease contact us with your objects so that we can resolve",
             " the bug.")
      }
      
      pending_nodes <- nodes_to_retry
    }
  }
  
  # Convert the merge list to a matrix
  merge_matrix <- do.call(rbind, merge_list)
  
  # Ensure that the height vector is in increasing order
  if (any(diff(height_list) < 0)) {
    stop("Heights must be in increasing order.")
  }
  
  # Function to perform a depth-first traversal to get the order
  order <- numeric()
  traverse_order <- function(node_id) {
    if (node_id < 0) {
      # Leaf node
      idx <- which(label_to_index == node_id)
      order[length(order) + 1] <<- idx
    } else {
      # Internal node
      merge_idx <- node_id
      left_id <- merge_matrix[merge_idx, 1]
      right_id <- merge_matrix[merge_idx, 2]
      traverse_order(left_id)
      traverse_order(right_id)
    }
  }
  
  # Start traversal from the last internal node
  root_node_id <- merge_index - 1
  traverse_order(root_node_id)
  
  # Construct the hclust object
  hclust_obj <- list(
    merge = merge_matrix,
    height = height_list,
    order = order,
    labels = labels,
    method = "Iterative Hierarchical Consensus Tree"
  )
  class(hclust_obj) <- "hclust"
  
  # FINALLY!
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

