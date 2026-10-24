#' Plot a hierarchical clustering bioregionalization
#'
#' Plots the output of a clustering function. For hierarchical clustering
#' ([hclu_hierarclust()], [hclu_diana()] and [cut_tree()]), the tree is drawn
#' with coloured bars under it showing the cluster of each site. Plots are not
#' available for the other clustering algorithms. Arguments starting with
#' `tree_` apply to trees only.
#'
#' @param x A `bioregion.clusters` object.
#'
#' @param tree_branch_partition Trees only. The partition
#' whose clusters give their colours to the branches of the tree: either a
#' partition number (`1` for the first partition, as in
#' `map_bioregions()`) or the name of a partition (e.g., `"K_8"`). By
#' default (`NULL`), the finest
#' partition (i.e. with the most clusters) is used.
#'
#' @param tree_labels Trees only. A `boolean` indicating
#' whether the site names should be written next to the bars (default `TRUE`).
#'
#' @param tree_flip Trees only. A `boolean` indicating
#' whether the order of the sites should be reversed, i.e., the tree drawn as
#' in a mirror (default `FALSE`). The tree remains the same: only the two
#' branches of each node are swapped.
#'
#' @param tree_horizontal Trees only. A `boolean` indicating
#' whether the tree should be drawn horizontally, with its root on the left and
#' the sites on the right, which makes site names easier to read (default
#' `FALSE`).
#'
#' @param ... Trees only. Further arguments passed to
#' [graphics::text()] to customize the site names (e.g., `cex` to change their
#' size).
#'
#' @details
#' **Tree plots (currently, only `hclu_` functions)**
#'
#' The tree is drawn with one coloured bar for each partition (i.e., each cut
#' of the tree), from coarsest to finest cuts. Each bar is labelled with
#' the name of its partition. In each bar, the colour next to a site shows the
#' cluster of that site. Sites that were not assigned to any cluster (`NA`)
#' are left blank.
#'
#' The branches that contain sites of a single cluster of
#' `tree_branch_partition` are drawn in the colour of that cluster. The other
#' branches, which join several clusters, are black. When the tree was cut at
#' given heights (or heights found for given numbers of clusters), dashed lines
#' show these heights.
#'
#' With the dynamic tree cut (`cut_tree(dynamic_tree_cut = TRUE)`), the
#' clusters are not obtained by cutting the whole tree at a single height,
#' and a cluster is not always a whole branch of the tree: it can be made of
#' neighbouring sites hanging from different branches. With
#' `dynamic_method = "hybrid"`, a cluster can even appear in several places
#' of the tree. The bars show these clusters as they are.
#'
#' The colours are those stored in `x` by [bioregion_colors()], so that the
#' tree and the maps made with [map_bioregions()] share the same colours. If
#' `x` has no colours, or if its colours do not match its partitions (e.g.,
#' after cutting the tree again), colours are assigned with
#' [bioregion_colors()] and its default parameters.
#'
#' The site names are written in the bottom margin of the plot, or in its right
#' margin if `tree_horizontal = TRUE`. Long site names may not fit in it:
#' either reduce their size (e.g., `cex = 0.5`), enlarge the margin with
#' [graphics::par()] before plotting (e.g., `par(mar = c(8, 4, 1, 1))`, or
#' `par(mar = c(4, 1, 4, 8))` for a horizontal tree), or remove them with
#' `tree_labels = FALSE`.
#'
#'
#' @seealso
#' For more details illustrated with a practical example,
#' see the vignette:
#' \url{https://biorgeo.github.io/bioregion/articles/a4_1_hierarchical_clustering.html}.
#'
#' Associated functions:
#' [hclu_hierarclust] [cut_tree] [bioregion_colors] [map_bioregions]
#'
#' @author
#' Boris Leroy (\email{leroy.boris@@gmail.com}) \cr
#' Maxime Lenormand (\email{maxime.lenormand@@inrae.fr}) \cr
#' Pierre Denelle (\email{pierre.denelle@@gmail.com})
#'
#' @examples
#' data(fishmat)
#' dissim <- dissimilarity(fishmat, metric = "Simpson")
#'
#' # One partition
#' tree1 <- hclu_hierarclust(dissim, n_clust = 8,
#'                           optimal_tree_method = "best", verbose = FALSE)
#' plot(tree1, cex = 0.3)
#'
#' # Several partitions, branches coloured with the partition with the fewest
#' # clusters (column 2 of tree2$clusters)
#' tree2 <- cut_tree(tree1, n_clust = c(4, 8, 16), verbose = FALSE)
#' plot(tree2, tree_branch_partition = 2, tree_labels = FALSE)
#'
#' # Same tree, mirrored
#' plot(tree2, tree_flip = TRUE, tree_labels = FALSE)
#'
#' # Dynamic tree cut
#' tree3 <- cut_tree(tree1, dynamic_tree_cut = TRUE, verbose = FALSE)
#' plot(tree3, tree_labels = FALSE)
#'
#'
#' # Same colours as on maps
#' tree3 <- bioregion_colors(tree2,
#'                           palette = "Antique")
#' plot(tree3, tree_labels = FALSE)
#' \dontrun{
#' data(fishsf)
#' map_bioregions(tree3, fishsf)
#' }
#'
#' # Horizontal tree, with room for the site names in the right margin
#' small_mat <- fishmat[1:30, ]
#' small_mat <- small_mat[, colSums(small_mat) > 0]
#' tree4 <- hclu_hierarclust(dissimilarity(small_mat, metric = "Simpson"),
#'                           n_clust = c(3, 6), optimal_tree_method = "best",
#'                           verbose = FALSE)
#' old_par <- par(mar = c(4, 1, 4, 7))
#' plot(tree4, tree_horizontal = TRUE, cex = 0.7)
#' par(old_par)
#'
#' @export
#' @method plot bioregion.clusters
plot.bioregion.clusters <- function(x,
                                    tree_branch_partition = NULL,
                                    tree_labels = TRUE,
                                    tree_flip = FALSE,
                                    tree_horizontal = FALSE,
                                    ...) {

  # Each type of algorithm has its own plot, made by an internal function.
  # Plots for other algorithms are to be added here, for instance:
  # - hclu_optics(): reachability plot
  # - nhclu_*(): sites in the space of the (dis)similarities
  # - netclu_*(): network
  # Arguments are named with the prefix of the type of plot they apply to
  # (e.g., tree_ for trees, which could later also include the hierarchical
  # trees of netclu_infomap()).
  if (x$name %in% c("hclu_hierarclust", "hclu_diana")) {
    plot_tree(x,
              tree_branch_partition = tree_branch_partition,
              tree_labels = tree_labels,
              tree_flip = tree_flip,
              tree_horizontal = tree_horizontal,
              ...)
  } else {
    stop(paste0("plot() is not available yet for the outputs of ", x$name,
                "()."),
         call. = FALSE)
  }

  invisible(NULL)
}


# Tree of hierarclust() or diana(), with one coloured bar per partition
plot_tree <- function(x, tree_branch_partition, tree_labels, tree_flip,
                      tree_horizontal, ...) {

  controls(args = tree_labels, data = NULL, type = "boolean")
  controls(args = tree_flip, data = NULL, type = "boolean")
  controls(args = tree_horizontal, data = NULL, type = "boolean")

  tree <- x$algorithm$final.tree
  if (x$name == "hclu_diana") {
    tree <- stats::as.hclust(tree)
  }
  # Reversing the order of the sites swaps the two branches of every node:
  # the tree is the same, drawn as in a mirror
  if (tree_flip) {
    tree$order <- rev(tree$order)
  }
  n <- length(tree$order)
  max_h <- max(tree$height)

  # Partitions, from the fewest to the most clusters, with sites in the order
  # of tree$labels
  partitions <- character(0)
  if (is.data.frame(x$clusters)) {
    clusters <- x$clusters[match(tree$labels, x$clusters[, 1]), ,
                           drop = FALSE]
    n_clust <- vapply(clusters[, -1, drop = FALSE],
                      function(cl) length(unique(stats::na.omit(cl))),
                      numeric(1))
    partitions <- colnames(clusters)[-1][order(n_clust)]

    if (is.null(tree_branch_partition)) {
      tree_branch_partition <- partitions[length(partitions)]
    } else {
      if (length(tree_branch_partition) != 1) {
        stop("tree_branch_partition must be of length 1.", call. = FALSE)
      }
      controls(args = tree_branch_partition, data = x$clusters,
               type = "input_partition_index")
      if (is.numeric(tree_branch_partition)) {
        tree_branch_partition <- colnames(x$clusters)[-1][tree_branch_partition]
      }
    }

    # Colours stored by bioregion_colors(), or new ones if they are missing
    # or do not match the partitions (the same as in map_bioregions())
    colors <- cluster_colors(x)
    # Colour of each site (rows) in each partition (columns)
    site_colors <- vapply(partitions, function(p) {
      colors[[p]]$color[match(clusters[, p], colors[[p]]$cluster)]
    }, character(n))
  }
  n_part <- length(partitions)

  # Heights of the bars, below the sites (h < 0)
  bar_h <- 0.045 * max_h
  gap <- 0.015 * max_h
  bars_bottom <- -gap - n_part * bar_h

  # Everything is positioned with two coordinates: the position along the
  # sites ("pos", from 1 to n, in the order of tree$order) and the height
  # ("h", the merge heights of the tree, bars below 0).
  # Vertical tree: pos on the x-axis, h on the y-axis.
  # Horizontal tree: h on the x-axis, reversed so that the root is on the
  # left, and pos on the y-axis, reversed so that the first site is on top.
  draw_segments <- function(pos0, h0, pos1, h1, ...) {
    if (tree_horizontal) {
      graphics::segments(h0, pos0, h1, pos1, ...)
    } else {
      graphics::segments(pos0, h0, pos1, h1, ...)
    }
  }
  draw_rect <- function(pos0, h0, pos1, h1, ...) {
    if (tree_horizontal) {
      graphics::rect(h0, pos0, h1, pos1, ...)
    } else {
      graphics::rect(pos0, h0, pos1, h1, ...)
    }
  }

  graphics::plot.new()
  h_axis <- pretty(c(0, max_h))
  h_title <- paste0(x$args$index, " dissimilarity")
  cut_height <- x$algorithm$output_cut_height
  if (tree_horizontal) {
    graphics::plot.window(xlim = c(max_h * 1.02, bars_bottom),
                          ylim = c(n + 0.5, 0.5),
                          xaxs = "i", yaxs = "i")
    graphics::axis(1, at = h_axis)
    graphics::mtext(h_title, side = 1, line = 2.5)
    if (is.numeric(cut_height)) {
      graphics::abline(v = cut_height, lty = 2, col = "grey55")
    }
  } else {
    graphics::plot.window(xlim = c(0.5, n + 0.5),
                          ylim = c(bars_bottom, max_h * 1.02),
                          xaxs = "i", yaxs = "i")
    graphics::axis(2, at = h_axis, las = 1)
    graphics::mtext(h_title, side = 2, line = 3)
    if (is.numeric(cut_height)) {
      graphics::abline(h = cut_height, lty = 2, col = "grey55")
    }
  }

  # Principles for building the tree: the position of a site is its rank in
  # tree$order, and the position of a merge is the mean of the positions of
  # its two children, at the merge height.
  # The tree is composed of segments based on pos and h and is constructed
  # from bottom to top. This allows us to reuse pos to draw colour rectangles
  # for clusters at the bottom.

  # Default colour: we put black to not mess up with bioregion_colors()'s grey
  # nuances for smaller clusters.
  # Black here either means the site has no cluster, or it is an
  # insignificant cluster as per user choice in bioregion_colors(), or the
  # tree is not coloured at all and all branches are black.
  default_col <- "black"
  site_pos <- numeric(n)
  site_pos[tree$order] <- seq_len(n)
  site_cl <- rep(NA_character_, n)
  site_col <- rep(default_col, n)
  if (n_part > 0) {
    site_cl <- clusters[, tree_branch_partition]
    site_col[!is.na(site_cl)] <- site_colors[!is.na(site_cl),
                                             tree_branch_partition]
  }
  node_pos <- numeric(n - 1)
  node_cl <- rep(NA_character_, n - 1)
  node_col <- rep(default_col, n - 1)
  # Loop through all merges
  for (i in seq_len(n - 1)) {
    child_pos <- child_h <- numeric(2)
    child_cl <- child_col <- character(2)
    # Loop through the two branches of the merge
    for (j in 1:2) {
      k <- tree$merge[i, j]
      # Negative means it's a site
      if (k < 0) {
        child_pos[j] <- site_pos[-k]
        child_h[j] <- 0 # So height 0
        child_cl[j] <- site_cl[-k]
        child_col[j] <- site_col[-k]
      } else { # Positive means it's another merge
        child_pos[j] <- node_pos[k]
        child_h[j] <- tree$height[k]
        child_cl[j] <- node_cl[k]
        child_col[j] <- node_col[k]
      }
    }
    # Merge at the average between children positions
    node_pos[i] <- mean(child_pos)
    if (!anyNA(child_cl) && child_cl[1] == child_cl[2]) {
      node_cl[i] <- child_cl[1]
      node_col[i] <- child_col[1]
    }
    # Each child is linked to its node by a line up to the height of the
    # node, then a line across to the position of the node
    draw_segments(child_pos, child_h, child_pos, tree$height[i],
                  col = child_col)
    draw_segments(child_pos, tree$height[i], node_pos[i], tree$height[i],
                  col = child_col)
  }

  # Rectangles below the tree for cluster colours
  # Loop through all partitions
  for (p in seq_len(n_part)) {
    # Define rectangle heights first
    bar_top <- -gap - (p - 1) * bar_h
    bar_bottom <- bar_top - 0.85 * bar_h

    cl <- clusters[tree$order, partitions[p]]
    col <- site_colors[tree$order, partitions[p]]
    # One rectangle per run of neighbouring sites of the same cluster
    # This method allows to have clusters split throughout the tree
    # (which can happen with dynamic tree cuts)
    # Each unassigned site is left blank
    run_id <- ifelse(is.na(cl), paste0("NA_", seq_len(n)), cl)
    # The next line checks where new clusters start among sites
    # by checking the list without first site vs. list without last site
    # If different = means start of a new cluster
    run_start <- which(c(TRUE, run_id[-1] != run_id[-n]))
    # ofc, cluster ends before the next start
    run_end <- c(run_start[-1] - 1, n)
    assigned <- !is.na(cl[run_start])
    # Positions of rectangles are at the positions of sites +- 0.5
    draw_rect(run_start[assigned] - 0.5, bar_bottom,
              run_end[assigned] + 0.5, bar_top,
              col = col[run_start[assigned]], border = "white", lwd = 0.6)
    draw_rect(0.5, bar_bottom, n + 0.5, bar_top, border = "grey70",
              lwd = 0.6)
    # Partition names: left of the bars, or above them if horizontal
    if (tree_horizontal) {
      graphics::text((bar_top + bar_bottom) / 2, 0.5, partitions[p],
                     srt = 90, adj = c(-0.1, 0.5), xpd = NA, cex = 0.75)
    } else {
      graphics::text(0.5, (bar_top + bar_bottom) / 2, partitions[p],
                     adj = c(1.1, 0.5), xpd = NA, cex = 0.75)
    }
  }

  # Site names: under the bars, or on their right if horizontal
  if (tree_labels) {
    if (tree_horizontal) {
      graphics::text(bars_bottom - 0.01 * max_h, seq_len(n),
                     tree$labels[tree$order], adj = c(0, 0.5), xpd = NA,
                     ...)
    } else {
      graphics::text(seq_len(n), bars_bottom - 0.01 * max_h,
                     tree$labels[tree$order], srt = 90, adj = c(1, 0.5),
                     xpd = NA, ...)
    }
  }
}
