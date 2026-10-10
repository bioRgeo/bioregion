# Pure-R references for the C++ helpers of R/ihct.R (src/ihct_trees.cpp).
#
# They are written the slow and obvious way -- the sites of every node are
# carried around as a set and the heights are read off the dissimilarity matrix
# -- so that the fast versions can be checked against something that shares no
# code with them.
#
# A tree is the merge matrix and heights of an hclust object plus `leaf_site`,
# the site of the dissimilarity matrix (its row number) that each leaf stands
# for; see the notes at the top of R/ihct.R.

# the sites under every node of a tree
ref_node_members <- function(merge, leaf_site) {
  members <- vector("list", nrow(merge))
  for (k in seq_len(nrow(merge))) {
    members[[k]] <- unlist(lapply(merge[k, ], function(j) {
      if (j < 0) leaf_site[-j] else members[[j]]
    }))
  }
  members
}

# the two groups of sites the root of a tree joins
ref_top_division <- function(merge, leaf_site) {
  members <- ref_node_members(merge, leaf_site)
  lapply(merge[nrow(merge), ], function(j) if (j < 0) leaf_site[-j] else members[[j]])
}

# A tree as one string per node: its sites, its height and the number of site
# pairs it joins. Sorted, so that two trees can be compared whatever the order
# of their rows and whatever their leaf numbering.
ref_tree_key <- function(merge, height, pairs, leaf_site) {
  members <- ref_node_members(merge, leaf_site)
  sort(sprintf("%s|%.12f|%g",
               vapply(members, function(x) paste(sort(x), collapse = ","), ""),
               height, pairs))
}

# The same key for the tree that pruning a tree to the sites `keep` is TRUE for
# must produce: a node survives when both of the groups it joins keep at least
# one site, and its height is then the mean dissimilarity between what is left
# of them.
ref_prune_tree <- function(merge, leaf_site, keep, dist_mat) {
  members <- ref_node_members(merge, leaf_site)
  out <- character(0)
  for (k in seq_len(nrow(merge))) {
    sides <- lapply(merge[k, ], function(j) if (j < 0) leaf_site[-j] else members[[j]])
    sides <- lapply(sides, function(sites) sites[keep[sites]])
    if (!length(sides[[1]]) || !length(sides[[2]])) next
    out <- c(out, sprintf("%s|%.12f|%g",
                          paste(sort(c(sides[[1]], sides[[2]])), collapse = ","),
                          mean(dist_mat[sides[[1]], sides[[2]]]),
                          length(sides[[1]]) * length(sides[[2]])))
  }
  sort(out)
}
