#' Plot a bioregionalization
#'
#' Plots the clusters of a bioregionalization, whatever the algorithm used to
#' obtain them. Three types of plots are available: bars showing the size of
#' each cluster (all algorithms), trees (hierarchical clustering and nested
#' partitions, such as the levels of [netclu_infomap()]), and trees of the
#' clusters, with the bar of each cluster at the tip of its branch.
#'
#' @param x A `bioregion.clusters` object.
#'
#' @param plot_type A `character` string indicating the type of plot:
#' `"bars"`, `"tree"`, `"cluster_tree"` or `"auto"` (default). With `"auto"`,
#' the plot is a tree for [hclu_hierarclust()] and [hclu_diana()], a tree of
#' the clusters for hierarchical results (e.g., [netclu_infomap()] with
#' several levels), and bars otherwise. See Details.
#'
#' @param partitions The partitions to plot, given by their number (`1` for
#' the first partition, i.e. the first column of clusters in `x$clusters`, `2`
#' for the second, etc.) or by their names (e.g., `c("K_3", "K_8")`). By
#' default (`NULL`), all partitions are plotted. For trees, the partition with
#' the most clusters is the one at the tips of the tree.
#'
#' @param heights Trees only. A `character` string indicating the heights of
#' the branches of trees: `"tree"` (the heights of the tree of
#' [hclu_hierarclust()] or [hclu_diana()]), `"rank"` (the merges of that tree
#' in the same order, but evenly spaced; trees of clusters only), `"levels"`
#' (one step per partition) or `"auto"` (default: `"tree"` when the tree is
#' available, `"levels"` otherwise). See Details.
#'
#' @param collapse A `boolean` indicating whether the clusters at the tips of
#' a tree (`plot_type = "tree"`) should be collapsed into triangles (default
#' `FALSE`).
#'
#' @param branch_partition Trees only. The partition whose clusters give their
#' colors to the branches, given by its number or name (see `partitions`). By
#' default (`NULL`), the partition at the tips of the tree.
#'
#' @param color_partition Bars and trees of clusters only. The partition whose
#' clusters give their colors to the squares, given by its number or name (see
#' `partitions`). By default (`NULL`), each partition is colored by its own
#' clusters.
#'
#' @param items_per_unit Bars and trees of clusters only. The number of sites
#' drawn by each square and of species drawn by each circle: `"auto"`
#' (default; one, or more for large datasets), one number for both, or two
#' numbers (sites, species).
#'
#' @param bipartite_layout Bars of bipartite clusterings only. A `character`
#' string indicating how sites and species are shown: `"mirror"` (default;
#' sites above the baseline, species below), `"stack"` (species above the
#' sites, in the same bar) or `"side"` (sites and species in two bars, side by
#' side).
#'
#' @param max_clusters Bars only. The maximum number of bars (default 25).
#' The smallest clusters beyond this number are grouped in a last bar,
#' `"other"`.
#'
#' @param scales Bars of several partitions only. A `character` string
#' indicating whether the axis of the numbers of sites is the same in all
#' panels (`"fixed"`, default), so that panels can be compared, or adapted to
#' each panel (`"free"`).
#'
#' @param horizontal A `boolean` indicating whether the plot should be drawn
#' horizontally (default `FALSE`): bars from left to right, trees with their
#' root on the left.
#'
#' @param flip A `boolean` indicating whether the order of the bars or of the
#' branches of the trees should be reversed (default `FALSE`). A tree remains
#' the same: it is drawn as in a mirror.
#'
#' @param labels A `boolean` indicating whether the names of the sites
#' (trees) or of the clusters (bars and trees of clusters) should be written
#' (default `TRUE`).
#'
#' @param axes A `boolean` indicating whether the axis of the numbers of sites
#' (bars) or of the heights (trees) should be drawn (default `TRUE`).
#'
#' @param main A `character` string giving a title. For bars of several
#' partitions, one title per partition can be given. By default (`NULL`), bars
#' are titled with the name of their partition and trees have no title.
#'
#' @param xlab,ylab `character` strings giving the titles of the horizontal
#' and vertical axes. By default (`NULL`), the axis of the numbers of sites or
#' of the heights has a title, and the other axis has none. Use `""` to remove
#' a title.
#'
#' @param ... Further arguments passed to [graphics::text()] to customize the
#' names of the sites (trees) or of the clusters (bars and trees of clusters),
#' e.g., `cex` to change their size.
#'
#' @details
#' **Bars (`plot_type = "bars"`)**
#'
#' Each cluster is a bar made of small squares, one per site (or per group of
#' sites, see `items_per_unit`), and its number of sites is written at the end
#' of the bar. For bipartite clusterings, species are drawn as circles.
#' Clusters are sorted by decreasing number of sites, as in
#' [bioregion_colors()], and clusters of the same size by name (`2.1`, `2.2`,
#' ..., `2.10`). The smallest clusters beyond `max_clusters` are grouped in a
#' wide bar, `"other"`, as high as the highest bar: they are not one cluster,
#' so their total size is not compared with the other bars. With several
#' partitions, each partition is plotted in its own panel, with the same axis
#' of the numbers of sites in all panels (`scales = "fixed"`) or with an axis
#' adapted to each panel (`scales = "free"`). The squares are always square,
#' so their size can differ between panels. The axis does not apply to the
#' `"other"` bar, whose squares are spread over a wider bar.
#'
#' With `color_partition`, the squares take the colors of another partition.
#' With a partition with more clusters, the bars show how each cluster is
#' split; with a partition with fewer clusters, bars of the same color show
#' which clusters are grouped together (they are placed side by side).
#'
#' **Trees (`plot_type = "tree"`)**
#'
#' The tree is drawn with one colored bar for each partition, from the fewest
#' to the most clusters. In each bar, the color next to a site shows the
#' cluster of that site. Sites that were not assigned to any cluster (`NA`)
#' are left blank. The branches that contain sites of a single cluster of
#' `branch_partition` are drawn in the color of that cluster, the other
#' branches are black. When the tree was cut at given heights (or heights
#' found for given numbers of clusters), dashed lines show these heights.
#'
#' With `collapse = TRUE`, each cluster of the partition at the tips is drawn
#' as a triangle, as wide as its number of sites and as high as the height at
#' which its sites are joined. Clusters must then be single branches of the
#' tree: this is not always the case with the dynamic tree cut
#' (`cut_tree(dynamic_tree_cut = TRUE)`), whose clusters cannot be collapsed.
#'
#' **Trees of clusters (`plot_type = "cluster_tree"`)**
#'
#' The tree is drawn down to the clusters of the partition at the tips only,
#' and the bar of each cluster is drawn at the tip of its branch, with its
#' name and number of sites (and species). The other partitions are shown by
#' dashed lines (`heights = "tree"` or `"rank"`) or by the levels of the tree
#' (`heights = "levels"`).
#'
#' **Heights of the trees**
#'
#' - `"tree"`: the heights of the tree of [hclu_hierarclust()] or
#'   [hclu_diana()], i.e. the dissimilarities at which sites and clusters are
#'   joined. For trees of clusters, the height axis starts just below the
#'   lowest cut, to spread the merges between clusters.
#' - `"rank"`: the same tree, with its merges in the same order, but evenly
#'   spaced. Useful when the merges between clusters are crowded at the top
#'   of the tree. Trees of clusters only.
#' - `"levels"`: a tree rebuilt from the partitions, with one step per
#'   partition. A cluster split into more than two clusters at the next level
#'   is drawn as a rake, i.e. with all its branches at the same height. This
#'   is the only option for algorithms without a tree, such as
#'   [netclu_infomap()]. The partitions must be nested: each cluster of a
#'   partition must be included in a single cluster of the partitions with
#'   fewer clusters. Sites without a cluster (`NA`) in any of the partitions
#'   are left out of the tree.
#'
#' **Arguments that do not apply to a plot**
#'
#' Some arguments only apply to some types of plots (e.g., `max_clusters` to
#' bars, `collapse` to trees). When such an argument is changed from its
#' default value but does not apply to the plot, it is ignored with a
#' warning.
#'
#' **Colors**
#'
#' The colors are those stored in `x` by [bioregion_colors()], so that plots
#' and maps made with [map_bioregions()] share the same colors. If `x` has no
#' colors, or if its colors do not match its partitions (e.g., after cutting
#' the tree again), colors are assigned with [bioregion_colors()] and its
#' default parameters.
#'
#' **Names of sites and clusters**
#'
#' Site names are written in the bottom margin of trees, or in their right
#' margin if `horizontal = TRUE`. Long names may not fit in it: either reduce
#' their size (e.g., `cex = 0.5`), enlarge the margin with [graphics::par()]
#' before plotting (e.g., `par(mar = c(8, 4, 1, 1))`), or remove them with
#' `labels = FALSE`.
#'
#' @return
#' Nothing; the plot is drawn on the current graphics device.
#'
#' @seealso
#' For more details illustrated with a practical example,
#' see the vignette:
#' \url{https://biorgeo.github.io/bioregion/articles/a5_1_visualization.html}.
#'
#' Associated functions:
#' [hclu_hierarclust] [cut_tree] [netclu_infomap] [bioregion_colors]
#' [map_bioregions]
#'
#' @author
#' Boris Leroy (\email{leroy.boris@@gmail.com}) \cr
#' Maxime Lenormand (\email{maxime.lenormand@@inrae.fr}) \cr
#' Pierre Denelle (\email{pierre.denelle@@gmail.com})
#'
#' @examples
#' data(fishmat)
#' data(fishdf)
#' dissim <- dissimilarity(fishmat, metric = "Simpson")
#'
#' # Bars, for any algorithm
#' clust <- nhclu_pam(dissim, n_clust = c(4, 8))
#' plot(clust)
#' plot(clust, scales = "free")
#' plot(clust, partitions = 2, horizontal = TRUE, main = "PAM, 8 clusters")
#'
#' # Bars of a bipartite clustering: sites (squares) and species (circles)
#' clust_bip <- netclu_greedy(fishdf, bipartite = TRUE)
#' plot(clust_bip)
#' plot(clust_bip, bipartite_layout = "stack")
#'
#' # Tree of a hierarchical clustering, cut into 3, 6 and 12 clusters
#' tree <- hclu_hierarclust(dissim, n_clust = c(3, 6, 12),
#'                          optimal_tree_method = "best", verbose = FALSE)
#' plot(tree, labels = FALSE)
#'
#' # Same tree, with the 12 clusters collapsed
#' plot(tree, collapse = TRUE, labels = FALSE)
#'
#' # Tree of the 6 clusters of the second partition, with their sizes
#' plot(tree, plot_type = "cluster_tree", partitions = 1:2,
#'      horizontal = TRUE)
#'
#' # Same, with the merges evenly spaced
#' plot(tree, plot_type = "cluster_tree", partitions = 1:2,
#'      heights = "rank", horizontal = TRUE)
#'
#' # Bars of the 12 clusters, colored by the 3 clusters they belong to
#' plot(tree, plot_type = "bars", partitions = 3, color_partition = 1)
#'
#' \dontrun{
#' # Hierarchy of the levels of Infomap
#' clust_infomap <- netclu_infomap(fishdf, bipartite_version = TRUE,
#'                                 show_hierarchy = TRUE)
#' plot(clust_infomap, horizontal = TRUE)
#' plot(clust_infomap, plot_type = "tree", collapse = TRUE)
#' }
#'
#' @export
#' @method plot bioregion.clusters
plot.bioregion.clusters <- function(x,
                                    plot_type = "auto",
                                    partitions = NULL,
                                    heights = "auto",
                                    collapse = FALSE,
                                    branch_partition = NULL,
                                    color_partition = NULL,
                                    items_per_unit = "auto",
                                    bipartite_layout = "mirror",
                                    max_clusters = 25,
                                    scales = "fixed",
                                    horizontal = FALSE,
                                    flip = FALSE,
                                    labels = TRUE,
                                    axes = TRUE,
                                    main = NULL,
                                    xlab = NULL,
                                    ylab = NULL,
                                    ...) {

  # Trees of hclu_hierarclust() and hclu_diana() can be plotted before they
  # are cut, i.e. without clusters
  has_tree <- inherits(x, "bioregion.clusters") &&
    isTRUE(x$name %in% c("hclu_hierarclust", "hclu_diana"))
  has_clusters <- is.data.frame(x$clusters)
  if (!has_tree || has_clusters) {
    controls(args = NULL, data = x, type = "input_bioregionalization")
  }

  controls(args = plot_type, data = NULL, type = "character")
  if (!(plot_type %in% c("auto", "bars", "tree", "cluster_tree"))) {
    stop(paste0("plot_type must be one of 'auto', 'bars', 'tree' or ",
                "'cluster_tree'."), call. = FALSE)
  }
  controls(args = heights, data = NULL, type = "character")
  if (!(heights %in% c("auto", "tree", "rank", "levels"))) {
    stop("heights must be one of 'auto', 'tree', 'rank' or 'levels'.",
         call. = FALSE)
  }
  controls(args = collapse, data = NULL, type = "boolean")
  controls(args = horizontal, data = NULL, type = "boolean")
  controls(args = flip, data = NULL, type = "boolean")
  controls(args = labels, data = NULL, type = "boolean")
  controls(args = axes, data = NULL, type = "boolean")
  controls(args = scales, data = NULL, type = "character")
  if (!(scales %in% c("fixed", "free"))) {
    stop("scales must be one of 'fixed' or 'free'.", call. = FALSE)
  }
  if (!is.null(main)) {
    controls(args = main, data = NULL, type = "character_vector")
  }
  if (!is.null(xlab)) controls(args = xlab, data = NULL, type = "character")
  if (!is.null(ylab)) controls(args = ylab, data = NULL, type = "character")
  controls(args = bipartite_layout, data = NULL, type = "character")
  if (!(bipartite_layout %in% c("mirror", "stack", "side"))) {
    stop("bipartite_layout must be one of 'mirror', 'stack' or 'side'.",
         call. = FALSE)
  }
  controls(args = max_clusters, data = NULL, type = "strict_positive_integer")
  if (max_clusters < 2) {
    stop("max_clusters must be higher than 1.", call. = FALSE)
  }
  if (!identical(items_per_unit, "auto")) {
    if (!is.numeric(items_per_unit)) {
      stop(paste0("items_per_unit must be 'auto' or one or two strictly ",
                  "positive integers."), call. = FALSE)
    }
    controls(args = items_per_unit, data = NULL,
             type = "strict_positive_integer_vector")
    if (length(items_per_unit) > 2) {
      stop("items_per_unit must be of length 1 or 2.", call. = FALSE)
    }
  }

  # Uncut trees: only the tree itself can be plotted, without partitions
  if (!has_clusters) {
    if (plot_type %in% c("bars", "cluster_tree")) {
      stop(paste0("No clusters have been generated for your hierarchical ",
                  "tree, please extract clusters from the tree.\n",
                  "See ?hclu_hierarclust, ?hclu_diana or ?cut_tree."),
           call. = FALSE)
    }
    need_cut <- c(partitions = !is.null(partitions),
                  branch_partition = !is.null(branch_partition),
                  color_partition = !is.null(color_partition),
                  "heights = 'levels'" = heights == "levels")
    if (any(need_cut)) {
      stop(paste0(paste(names(need_cut)[need_cut], collapse = ", "),
                  " cannot be used before the tree is cut. ",
                  "See ?cut_tree."), call. = FALSE)
    }
  }

  # Partitions as names
  all_partitions <- if (has_clusters) colnames(x$clusters)[-1] else
    character(0)
  if (is.null(partitions)) {
    partitions <- all_partitions
  } else {
    controls(args = partitions, data = x$clusters,
             type = "input_partition_index")
    if (is.numeric(partitions)) partitions <- all_partitions[partitions]
  }
  if (!is.null(branch_partition)) {
    if (length(branch_partition) != 1) {
      stop("branch_partition must be of length 1.", call. = FALSE)
    }
    controls(args = branch_partition, data = x$clusters,
             type = "input_partition_index")
    if (is.numeric(branch_partition)) {
      branch_partition <- all_partitions[branch_partition]
    }
  }
  if (!is.null(color_partition)) {
    if (length(color_partition) != 1) {
      stop("color_partition must be of length 1.", call. = FALSE)
    }
    controls(args = color_partition, data = x$clusters,
             type = "input_partition_index")
    if (is.numeric(color_partition)) {
      color_partition <- all_partitions[color_partition]
    }
  }

  # Trees go from the fewest to the most clusters
  if (has_clusters) {
    n_clust <- vapply(partitions, function(p) {
      length(unique(stats::na.omit(x$clusters[, p])))
    }, numeric(1))
    tree_partitions <- partitions[order(n_clust)]
  } else {
    tree_partitions <- character(0)
  }
  nested <- length(tree_partitions) > 0 &&
    partitions_nested(x$clusters, tree_partitions)

  if (plot_type == "auto") {
    plot_type <- if (has_tree) {
      "tree"
    } else if (length(tree_partitions) > 1 && nested &&
               isTRUE(x$inputs$hierarchical)) {
      "cluster_tree"
    } else {
      "bars"
    }
  }

  # Arguments changed from their default but that do not apply to this plot
  bipartite <- has_clusters &&
    length(unique(attr(x$clusters, "node_type"))) > 1
  warn_ignored(plot_type,
               changed = c(heights = heights != "auto",
                           collapse = collapse,
                           branch_partition = !is.null(branch_partition),
                           color_partition = !is.null(color_partition),
                           items_per_unit = !identical(items_per_unit, "auto"),
                           bipartite_layout = bipartite_layout != "mirror",
                           max_clusters = max_clusters != 25,
                           scales = scales != "fixed"),
               bipartite = bipartite,
               n_partitions = length(partitions))

  if (plot_type == "bars") {
    plot_bars(x,
              partitions = partitions,
              color_partition = color_partition,
              items_per_unit = items_per_unit,
              bipartite_layout = bipartite_layout,
              max_clusters = max_clusters,
              scales = scales,
              horizontal = horizontal,
              flip = flip,
              labels = labels,
              axes = axes,
              main = main,
              xlab = xlab,
              ylab = ylab,
              ...)
    return(invisible(NULL))
  }

  if (heights == "auto") heights <- if (has_tree) "tree" else "levels"
  if (heights %in% c("tree", "rank") && !has_tree) {
    stop(paste0("heights = '", heights, "' needs the tree of ",
                "hclu_hierarclust() or hclu_diana(). Use heights = 'levels' ",
                "for the outputs of ", x$name, "()."), call. = FALSE)
  }
  if (heights == "rank" && plot_type != "cluster_tree") {
    stop("heights = 'rank' is only available for plot_type = 'cluster_tree'.",
         call. = FALSE)
  }
  if (heights == "levels" && !nested) {
    stop(paste0("The partitions are not nested: some clusters are split ",
                "between several clusters of a partition with fewer ",
                "clusters. Trees with heights = 'levels' can only be drawn ",
                "for nested partitions, e.g., the levels of netclu_infomap() ",
                "or several cuts of the same tree."), call. = FALSE)
  }
  if (is.null(branch_partition) && length(tree_partitions) > 0) {
    branch_partition <- tree_partitions[length(tree_partitions)]
  }

  plot_tree(x,
            partitions = tree_partitions,
            heights = heights,
            cluster_tree = plot_type == "cluster_tree",
            collapse = collapse,
            branch_partition = branch_partition,
            color_partition = color_partition,
            items_per_unit = items_per_unit,
            horizontal = horizontal,
            flip = flip,
            labels = labels,
            axes = axes,
            main = main,
            xlab = xlab,
            ylab = ylab,
            ...)

  invisible(NULL)
}

# Warning for the arguments changed from their default value (changed = TRUE)
# that do not apply to this type of plot, and are therefore ignored
warn_ignored <- function(plot_type, changed, bipartite, n_partitions) {
  used_by <- list(
    bars = c("color_partition", "items_per_unit", "bipartite_layout",
             "max_clusters", "scales"),
    tree = c("heights", "collapse", "branch_partition"),
    cluster_tree = c("heights", "branch_partition", "color_partition",
                     "items_per_unit"))
  ignored <- names(changed)[changed &
                              !(names(changed) %in% used_by[[plot_type]])]
  if (length(ignored) > 0) {
    several <- length(ignored) > 1
    warning(paste0(paste(ignored, collapse = ", "),
                   if (several) " do" else " does",
                   " not apply to plot_type = '", plot_type, "' and ",
                   if (several) "are" else "is", " ignored."),
            call. = FALSE)
  }
  if (plot_type == "bars") {
    if (changed[["bipartite_layout"]] && !bipartite) {
      warning(paste0("bipartite_layout only applies to bipartite ",
                     "clusterings (sites and species) and is ignored."),
              call. = FALSE)
    }
    if (changed[["scales"]] && n_partitions == 1) {
      warning(paste0("scales only applies to bars of several partitions ",
                     "and is ignored."), call. = FALSE)
    }
  }
}


# Bars -------------------------------------------------------------------------

# One panel per partition. Bars are drawn in numbers of sites along the bars
# (so that the axis shows numbers of sites) and in squares across the bars;
# each panel chooses the number of squares per row of its bars so that the
# squares are square and the bars fill the panel.
plot_bars <- function(x, partitions, color_partition, items_per_unit,
                      bipartite_layout, max_clusters, scales, horizontal,
                      flip, labels, axes, main, xlab, ylab, ...) {

  type <- attr(x$clusters, "node_type")
  bipartite <- length(unique(type)) > 1
  if (!bipartite) bipartite_layout <- "stack"
  k_unit <- units_per_item(items_per_unit, nrow(x$clusters))
  colors <- cluster_colors(x)

  content <- lapply(partitions, function(p) {
    fill_p <- if (is.null(color_partition)) p else color_partition
    bars_content(x$clusters[, p], x$clusters[, fill_p], type,
                 colors[[fill_p]], k_unit, max_clusters, bipartite_layout,
                 flip)
  })

  # Vertical bars side by side, horizontal bars one above the other
  if (length(partitions) > 1) {
    panels <- grDevices::n2mfrow(length(partitions))
    if (!horizontal) panels <- rev(panels)
    old_par <- graphics::par(mfrow = panels)
    on.exit(graphics::par(old_par))
  }

  # Length of the longest bars, in numbers of sites (species are converted
  # in the same units, so that squares and circles have the same size)
  up <- vapply(content, function(cnt) cnt$max_up * k_unit[["site"]],
               numeric(1))
  down <- vapply(content, function(cnt) cnt$max_down * k_unit[["site"]],
                 numeric(1))
  if (scales == "fixed") {
    up[] <- max(up)
    down[] <- max(down)
  }

  count_title <- if (!bipartite) {
    "Number of sites"
  } else if (bipartite_layout == "mirror") {
    "Number of sites (above) and species (below)"
  } else {
    "Number of sites and species"
  }
  main <- rep_len(if (is.null(main)) partitions else main, length(partitions))
  count_lab <- if (horizontal) xlab else ylab
  cluster_lab <- if (horizontal) ylab else xlab
  if (is.null(count_lab)) count_lab <- if (axes) count_title else ""

  for (j in seq_along(partitions)) {
    draw_bars(content[[j]], up = up[j], down = down[j], k_unit = k_unit,
              horizontal = horizontal, bipartite_layout = bipartite_layout,
              labels = labels, ...)
    if (axes) {
      bars_axis(up[j], down[j], k_unit, horizontal,
                mirror = bipartite_layout == "mirror")
    }
    graphics::title(xlab = if (horizontal) count_lab else cluster_lab,
                    ylab = if (horizontal) cluster_lab else count_lab)
    unit_txt <- units_text(k_unit, bipartite)
    if (!is.null(color_partition) && color_partition != partitions[j]) {
      unit_txt <- paste0(unit_txt, "; colors of ", color_partition)
    }
    graphics::mtext(main[j], side = 3, line = 1.2, adj = 0, cex = 0.9)
    graphics::mtext(unit_txt, side = 3, line = 0.2, adj = 0, cex = 0.7,
                    col = "grey30")
    if (bipartite) units_legend()
  }
}

# Clusters of a partition with their squares (sites) and circles (species).
# Clusters are sorted by decreasing number of sites (as in bioregion_colors())
# then by name; with the colors of another partition, clusters sharing their
# main color are put side by side.
bars_content <- function(cl, fill, type, fill_colors, k_unit, max_clusters,
                         bipartite_layout, flip) {
  cl <- as.character(cl)
  fill <- as.character(fill)
  fill_colors <- stats::setNames(fill_colors$color, fill_colors$cluster)
  keep <- !is.na(cl)

  n_sites <- tapply(type[keep] == "site", cl[keep], sum)
  clusters <- names(n_sites)
  clusters <- clusters[order(-n_sites, match(seq_along(clusters),
                                             natural_order(clusters)))]
  by_size <- clusters

  if (!identical(cl, fill)) {
    fill_n_sites <- tapply(type == "site", fill, sum)
    fill_rank <- names(fill_n_sites)[
      order(-fill_n_sites, match(seq_along(fill_n_sites),
                                 natural_order(names(fill_n_sites))))]
    main_fill <- vapply(clusters, function(k) {
      f <- table(fill[cl == k & !is.na(cl)])
      if (length(f) == 0) NA_character_ else names(f)[which.max(f)]
    }, character(1))
    clusters <- clusters[order(match(main_fill, fill_rank),
                               match(clusters, by_size))]
  }

  # The smallest clusters beyond max_clusters are grouped into a last bar,
  # under a name that differs from the names of the clusters (a cluster can
  # be named "other")
  n_other <- 0
  if (length(clusters) > max_clusters) {
    other <- by_size[max_clusters:length(by_size)]
    n_other <- length(other)
    other_key <- utils::tail(make.unique(c(by_size, "other")), 1)
    cl[cl %in% other] <- other_key
    clusters <- clusters[!(clusters %in% other)]
  }
  if (flip) clusters <- rev(clusters)
  if (n_other > 0) clusters <- c(clusters, other_key)
  is_other <- seq_along(clusters) == length(clusters) & n_other > 0

  units <- lapply(clusters, function(k) {
    cluster_units(which(cl == k), type, fill, fill_colors, k_unit)
  })
  stacks <- if (bipartite_layout == "stack") {
    lapply(units, function(u) list(all = do.call(rbind, u)))
  } else {
    units
  }
  up_types <- if (bipartite_layout == "mirror") "site" else
    c("all", "site", "species")
  stack_size <- function(s, t) {
    sum(vapply(s[names(s) %in% t], nrow, numeric(1)))
  }
  max_stack <- function(t) {
    max(c(0, unlist(lapply(stacks[!is_other], function(s) {
      vapply(s[names(s) %in% t], nrow, numeric(1))
    }))))
  }

  list(clusters = clusters, is_other = is_other, n_other = n_other,
       stacks = stacks, up_types = up_types,
       n_items = lapply(clusters, function(k) {
         table(factor(type[cl == k & !is.na(cl)], c("site", "species")))
       }),
       max_up = max_stack(up_types),
       max_down = if (bipartite_layout == "mirror") max_stack("species") else 0,
       other_up = if (n_other > 0) {
         s <- stacks[[which(is_other)]]
         if (bipartite_layout == "side") max(vapply(s, nrow, numeric(1))) else
           stack_size(s, up_types)
       } else {
         0
       },
       other_down = if (n_other > 0 && bipartite_layout == "mirror") {
         stack_size(stacks[[which(is_other)]], "species")
       } else {
         0
       },
       layout = bipartite_layout)
}

# Position of the bars for bars w squares wide; "other" is wider, so that it
# is as high as the highest bar
bars_geometry <- function(content, w) {
  rows_up <- ceiling(content$max_up / w)
  rows_down <- ceiling(content$max_down / w)
  wj <- rep(w, length(content$clusters))
  if (content$n_other > 0) {
    wj[content$is_other] <- max(
      w, ceiling(content$other_up / max(1, rows_up)),
      if (content$other_down > 0) ceiling(content$other_down /
                                            max(1, rows_down)))
  }
  slot <- if (content$layout == "side") 2 * wj + 0.6 else wj
  # Half a gap before the first bar, to keep it away from the axis
  gap <- max(1, round(0.6 * w))
  list(w = w, wj = wj, slot = slot,
       offsets = gap / 2 + cumsum(c(0, utils::head(slot + gap, -1))),
       span = sum(slot) + (length(slot) - 0.5) * gap,
       rows_up = rows_up, rows_down = rows_down)
}

# Bars of one partition, up to up sites above the baseline and down below.
# Along the bars, a row of w squares is w * k_unit["site"] sites long; across
# the bars, a square is 1 long. w is the largest number of squares per row for
# which the bars fit in the panel with square squares.
draw_bars <- function(content, up, down, k_unit, horizontal,
                      bipartite_layout, labels, ...) {
  mirror <- bipartite_layout == "mirror"
  # Room for the numbers at the ends of the bars
  margin <- 0.07 * (up + down)
  hi <- up + margin
  lo <- if (down > 0) -(down + margin) else 0
  pin <- graphics::par("pin")
  ratio <- if (horizontal) pin[2] / pin[1] else pin[1] / pin[2]
  across_len <- function(w) ratio * (hi - lo) / (w * k_unit[["site"]])
  fit <- bars_geometry(content, 1)
  for (w in seq_len(max(1, content$max_up))[-1]) {
    g <- bars_geometry(content, w)
    q <- w * k_unit[["site"]]
    if (g$span > across_len(w)) break
    # The last, incomplete row must not reach the top of the panel
    if (g$rows_up * q > up + 0.5 * margin ||
        g$rows_down * q > down + 0.5 * margin) next
    fit <- g
  }
  q <- fit$w * k_unit[["site"]]
  span <- max(fit$span, across_len(fit$w))
  start <- if (mirror) 0.15 * q else 0

  graphics::plot.new()
  if (horizontal) {
    graphics::plot.window(xlim = c(lo, hi), ylim = c(-span, 0),
                          xaxs = "i", yaxs = "i")
  } else {
    graphics::plot.window(xlim = c(0, span), ylim = c(lo, hi),
                          xaxs = "i", yaxs = "i")
  }
  n_cl <- length(content$clusters)
  wj <- fit$wj
  pad <- 0.015 * (hi - lo)

  for (j in seq_len(n_cl)) {
    offset <- fit$offsets[j]
    s <- content$stacks[[j]]
    for (m in seq_along(s)) {
      u <- s[[m]]
      species <- names(s)[m] == "species"
      dir <- if (mirror && species) -1 else 1
      off <- if (bipartite_layout == "side" && species) {
        offset + wj[j] + 0.6
      } else {
        offset
      }
      xy <- stack_xy(nrow(u), wj[j], off, dir, start, horizontal, q = q)
      draw_units(xy$x, xy$y, u$col, u$shape,
                 size = if (horizontal) c(0.82 * q, 0.82) else
                   c(0.82, 0.82 * q))
      # Number of items at the end of the stack
      n_end <- dir * (start + ceiling(nrow(u) / wj[j]) * q + pad)
      n_lab <- if (bipartite_layout == "stack") sum(content$n_items[[j]]) else
        content$n_items[[j]][[names(s)[m]]]
      mid <- off + min(wj[j], nrow(u)) / 2
      if (horizontal) {
        graphics::text(n_end, -mid, n_lab, cex = 0.65, col = "grey30",
                       adj = c(if (dir > 0) 0 else 1, 0.5), xpd = NA)
      } else {
        graphics::text(mid, n_end, n_lab, cex = 0.65, col = "grey30",
                       adj = c(0.5, if (dir > 0) 0 else 1), xpd = NA)
      }
    }
    if (!labels) next
    # Cluster name, centred on the squares actually drawn
    name <- if (content$is_other[j]) {
      paste0("other (", content$n_other, " clusters)")
    } else {
      content$clusters[j]
    }
    used <- min(wj[j], max(vapply(s, nrow, numeric(1))))
    centre <- if (bipartite_layout == "side") {
      offset + fit$slot[j] / 2
    } else {
      offset + used / 2
    }
    if (horizontal) {
      plot_text(lo - pad, -centre, name, ...,
                defaults = list(adj = c(1, 0.5), xpd = NA, cex = 0.75))
    } else {
      rotate <- !content$is_other[j] &&
        (n_cl > 12 || max(nchar(content$clusters)) > 4)
      plot_text(centre, lo - pad, name, ...,
                defaults = list(srt = if (rotate) 90 else 0,
                                adj = if (rotate) c(1, 0.5) else c(0.5, 1),
                                xpd = NA, cex = 0.75))
    }
  }
  if (mirror) {
    if (horizontal) graphics::abline(v = 0, col = "grey60") else
      graphics::abline(h = 0, col = "grey60")
  }
}

# Axis of the numbers of sites, and of species below the baseline in the
# mirror layout
bars_axis <- function(up, down, k_unit, horizontal, mirror) {
  at <- pretty(c(0, up))
  at <- at[at <= up + 0.07 * (up + down)]
  labels <- at
  if (mirror && down > 0) {
    # Species are drawn in numbers of sites: k_unit["species"] species take
    # the length of k_unit["site"] sites
    to_sites <- k_unit[["site"]] / k_unit[["species"]]
    species_at <- pretty(c(0, down / to_sites), n = 3)
    species_at <- species_at[species_at > 0 &
                               species_at * to_sites <= down]
    at <- c(-species_at * to_sites, at)
    labels <- c(species_at, labels)
  }
  graphics::axis(if (horizontal) 1 else 2, at = at, labels = labels, las = 1)
}


# Trees ------------------------------------------------------------------------

# Every tree is an hclust object: the tree of hclu_hierarclust() or
# hclu_diana(), or a tree rebuilt from nested partitions (heights = "levels").
# The same loop over the merges draws the full tree, the tree with collapsed
# clusters, and the tree of the clusters (cluster_tree = TRUE) with their bars.
plot_tree <- function(x, partitions, heights, cluster_tree, collapse,
                      branch_partition, color_partition, items_per_unit,
                      horizontal, flip, labels, axes, main, xlab, ylab,
                      ...) {

  if (heights == "levels") {
    tree <- nested_tree(x$clusters, partitions)
  } else {
    tree <- x$algorithm$final.tree
    if (x$name == "hclu_diana") tree <- stats::as.hclust(tree)
  }
  if (flip) tree$order <- rev(tree$order)
  n <- length(tree$labels)
  n_part <- length(partitions)
  type <- rep("site", n)
  if (n_part > 0) {
    rows <- match(tree$labels, x$clusters[, 1])
    clusters <- x$clusters[rows, , drop = FALSE]
    type <- attr(x$clusters, "node_type")[rows]
    colors <- cluster_colors(x)
    tips <- as.character(clusters[, partitions[n_part]])
    branch_cl <- as.character(clusters[, branch_partition])
  } else {
    tips <- branch_cl <- rep(NA_character_, n)
  }
  bipartite <- length(unique(type)) > 1
  collapse <- collapse && !cluster_tree && n_part > 0

  # Cluster of each merge: the cluster of all the sites under it, or NA
  tips_merge <- merge_clusters(tree, tips)
  if ((collapse || (cluster_tree && heights != "levels")) &&
      !clusters_are_branches(tree, tips, tips_merge)) {
    stop(paste0("The clusters of ", partitions[n_part], " span distinct ",
                "branches across the tree (as can happen with the dynamic ",
                "tree cut), so they cannot be ",
                if (collapse) "collapsed." else
                  "drawn as the tips of a tree of clusters."),
         call. = FALSE)
  }
  if (cluster_tree && anyNA(tips)) {
    stop(paste0("Some sites have no cluster in ", partitions[n_part],
                ", so the tree of its clusters cannot be drawn."),
         call. = FALSE)
  }

  # heights = "rank": the merges between clusters keep their order but are
  # evenly spaced; the merges inside clusters are put at 0
  between <- is.na(tips_merge)
  cut_h <- rep(NA_real_, n_part)
  # Heights of cut, in the order of the columns of x$clusters
  all_cut_h <- x$algorithm$output_cut_height
  if (heights != "levels" && is.numeric(all_cut_h) &&
      length(all_cut_h) == ncol(x$clusters) - 1) {
    cut_h <- unname(all_cut_h[match(partitions, colnames(x$clusters)[-1])])
  }
  if (heights == "rank") {
    steps <- sort(unique(tree$height[between]))
    cut_h <- vapply(cut_h, function(h) sum(steps < h) + 0.5, numeric(1))
    tree$height[between] <- match(tree$height[between], steps)
    tree$height[!between] <- 0
  }
  max_h <- max(tree$height)
  rank <- numeric(n)
  rank[tree$order] <- seq_len(n)

  # Coordinates: pos along the tips, h along the heights. Trees of clusters
  # are drawn in the units of the squares of their bars.
  if (cluster_tree) {
    terminals <- unique(tips[tree$order])
    fill <- if (is.null(color_partition)) tips else
      as.character(clusters[, color_partition])
    fill_colors <- colors[[if (is.null(color_partition)) partitions[n_part]
                           else color_partition]]
    fill_colors <- stats::setNames(fill_colors$color, fill_colors$cluster)
    k_unit <- units_per_item(items_per_unit, n)
    units <- lapply(terminals, function(k) {
      do.call(rbind, cluster_units(which(tips == k), type, fill, fill_colors,
                                   k_unit))
    })
    n_units <- vapply(units, nrow, numeric(1))
    n_t <- length(terminals)
    pin <- graphics::par("pin")
    aspect <- if (horizontal) pin[2] / pin[1] else pin[1] / pin[2]
    w <- if (horizontal) sqrt(1.9 * max(n_units) * aspect / (1.6 * n_t)) else
      sqrt(1.9 * aspect * max(n_units) / (1.6 * n_t))
    w <- max(1, min(round(w), max(n_units)))
    gap <- max(1, round(0.6 * w))
    offsets <- (seq_len(n_t) - 1) * (w + gap)
    centre <- offsets + pmin(w, n_units) / 2
    len <- ceiling(max(n_units) / w)
    span <- max(offsets) + w
    # Real heights start just below the lowest cut, to spread the merges
    # between clusters
    lower <- 0
    if (heights == "tree" && any(between)) {
      lower <- if (all(is.na(cut_h))) min(tree$height[between]) else
        min(cut_h, na.rm = TRUE)
      lower <- max(0, lower - 0.05 * (max_h - lower))
    }
    h_range <- if (max_h > lower) max_h - lower else 1
    base <- 0.3
    tree_len <- max(8, 0.7 * len)
    h_map <- function(h) base + (pmax(h, lower) - lower) / h_range * tree_len
    leaf_pos <- centre[match(tips, terminals)]
    to_xy <- if (horizontal) {
      function(p, h) list(x = -h_map(h), y = -p)
    } else {
      function(p, h) list(x = p, y = h_map(h))
    }
  } else {
    lower <- 0
    leaf_pos <- rank
    to_xy <- if (horizontal) {
      function(p, h) list(x = h, y = p)
    } else {
      function(p, h) list(x = p, y = h)
    }
  }
  draw_segments <- function(p0, h0, p1, h1, ...) {
    a <- to_xy(p0, h0)
    b <- to_xy(p1, h1)
    graphics::segments(a$x, a$y, b$x, b$y, ...)
  }

  # Window and height axis
  bar_h <- 0.045 * max_h
  gap_b <- 0.015 * max_h
  n_bars <- n_part + bipartite
  bars_bottom <- -gap_b - n_bars * bar_h
  graphics::plot.new()
  if (cluster_tree) {
    top <- h_map(max_h) * 1.02
    if (horizontal) {
      graphics::plot.window(xlim = c(-top, len * 1.15), ylim = c(-span, 0),
                            asp = 1)
    } else {
      graphics::plot.window(xlim = c(0, span), ylim = c(-len * 1.1, top),
                            asp = 1)
    }
  } else if (horizontal) {
    graphics::plot.window(xlim = c(max_h * 1.02, bars_bottom),
                          ylim = c(n + 0.5, 0.5), xaxs = "i", yaxs = "i")
  } else {
    graphics::plot.window(xlim = c(0.5, n + 0.5),
                          ylim = c(bars_bottom, max_h * 1.02),
                          xaxs = "i", yaxs = "i")
  }
  axis_side <- if (horizontal) 1 else 2
  axis_at <- function(h) {
    xy <- to_xy(rep(0, length(h)), h)
    if (horizontal) xy$x else xy$y
  }
  if (axes && heights == "tree") {
    at <- pretty(c(lower, max_h))
    at <- at[at >= lower & at <= max_h]
    graphics::axis(axis_side, at = axis_at(at), labels = at, las = 1)
  } else if (axes && heights == "levels") {
    lev_h <- n_part - seq_len(n_part) + 1
    graphics::axis(axis_side, at = axis_at(lev_h), labels = partitions,
                   las = 1, cex.axis = 0.8)
  }
  height_lab <- if (horizontal) xlab else ylab
  if (is.null(height_lab)) {
    height_lab <- if (axes && heights == "tree") {
      paste0(x$args$index, " dissimilarity")
    } else {
      ""
    }
  }
  other_lab <- if (horizontal) ylab else xlab
  graphics::title(xlab = if (horizontal) height_lab else other_lab,
                  ylab = if (horizontal) other_lab else height_lab)
  if (!is.null(main)) {
    graphics::mtext(main, side = 3, line = 1.2, adj = 0, cex = 0.9)
  }
  # Cuts of the partitions
  for (k in which(!is.na(cut_h))) {
    at <- axis_at(cut_h[k])
    usr <- graphics::par("usr")
    if (horizontal) {
      graphics::abline(v = at, lty = 2, col = "grey55")
      graphics::text(at, usr[4], partitions[k], srt = 90, adj = c(-0.1, 0.5),
                     cex = 0.65, col = "grey40", xpd = NA)
    } else {
      graphics::abline(h = at, lty = 2, col = "grey55")
      graphics::text(usr[2], at, partitions[k], adj = c(1, -0.3), cex = 0.65,
                     col = "grey40", xpd = NA)
    }
  }

  # Merges. The position of a site is its rank in tree$order (or the centre of
  # the bar of its cluster), and the position of a merge is the mean of the
  # positions of its two children. Each child is linked to its merge by a line
  # up to the height of the merge, then a line across to its position.
  # Branches holding the sites of a single cluster of branch_partition take
  # its color; black otherwise.
  site_col <- rep("black", n)
  if (n_part > 0) {
    branch_colors <- stats::setNames(colors[[branch_partition]]$color,
                                     colors[[branch_partition]]$cluster)
    site_col[!is.na(branch_cl)] <- branch_colors[branch_cl[!is.na(branch_cl)]]
  }
  branch_merge <- merge_clusters(tree, branch_cl)
  m <- nrow(tree$merge)
  node_pos <- node_min <- node_max <- numeric(m)
  node_col <- rep("black", m)
  if (n_part > 0) {
    node_col[!is.na(branch_merge)] <- branch_colors[
      branch_merge[!is.na(branch_merge)]]
  }
  hide_inside <- collapse || cluster_tree
  for (i in seq_len(m)) {
    child_pos <- child_h <- child_min <- child_max <- numeric(2)
    child_col <- character(2)
    for (j in 1:2) {
      k <- tree$merge[i, j]
      if (k < 0) {
        child_pos[j] <- leaf_pos[-k]
        child_h[j] <- 0
        child_col[j] <- site_col[-k]
        child_min[j] <- child_max[j] <- rank[-k]
      } else {
        child_pos[j] <- node_pos[k]
        child_h[j] <- tree$height[k]
        child_col[j] <- node_col[k]
        child_min[j] <- node_min[k]
        child_max[j] <- node_max[k]
      }
    }
    node_pos[i] <- mean(child_pos)
    node_min[i] <- min(child_min)
    node_max[i] <- max(child_max)
    if (hide_inside && !is.na(tips_merge[i])) next
    h <- rep(tree$height[i], 2)
    draw_segments(child_pos, child_h, child_pos, h, col = child_col)
    draw_segments(child_pos, h, rep(node_pos[i], 2), h, col = child_col)
  }

  # Highest merge of each collapsed cluster: a triangle down to its sites, or
  # a line down to its bar
  if (hide_inside) {
    parent <- integer(m)
    for (i in seq_len(m)) parent[tree$merge[i, tree$merge[i, ] > 0]] <- i
    tops <- which(!is.na(tips_merge) &
                    (parent == 0 | is.na(tips_merge[pmax(parent, 1)])))
    tips_colors <- stats::setNames(colors[[partitions[n_part]]]$color,
                                   colors[[partitions[n_part]]]$cluster)
    for (i in tops) {
      if (cluster_tree) {
        draw_segments(node_pos[i], 0, node_pos[i], tree$height[i],
                      col = node_col[i])
      } else {
        xy <- to_xy(c(node_pos[i], node_min[i] - 0.45, node_max[i] + 0.45),
                    c(tree$height[i], 0, 0))
        col <- tips_colors[[tips_merge[i]]]
        graphics::polygon(xy$x, xy$y, col = col, border = col)
      }
    }
  }

  if (cluster_tree) {
    for (j in seq_along(terminals)) {
      u <- units[[j]]
      xy <- stack_xy(nrow(u), w, offsets[j], 1, 0, horizontal)
      if (!horizontal) xy$y <- -xy$y
      draw_units(xy$x, xy$y, u$col, u$shape)
      n_j <- table(factor(type[tips == terminals[j]], c("site", "species")))
      lab <- if (bipartite) {
        sprintf("%s  (%d | %d)", terminals[j], n_j[[1]], n_j[[2]])
      } else {
        sprintf("%s  (%d)", terminals[j], n_j[[1]])
      }
      end <- ceiling(nrow(u) / w) + 0.4
      if (!labels) next
      if (horizontal) {
        plot_text(end, -centre[j], lab, ...,
                  defaults = list(adj = c(0, 0.5), cex = 0.7, xpd = NA))
      } else {
        plot_text(centre[j], -end, lab, ...,
                  defaults = list(adj = c(1, 0.5), srt = 90, cex = 0.7,
                                  xpd = NA))
      }
    }
    unit_txt <- units_text(k_unit, bipartite)
    if (bipartite) unit_txt <- paste0(unit_txt, "; (sites | species)")
    if (heights == "rank") unit_txt <- paste0(unit_txt,
                                              "; merges evenly spaced")
    graphics::mtext(unit_txt, side = 3, line = 0.2, adj = 0, cex = 0.7,
                    col = "grey30")
    if (bipartite) units_legend()
    return(invisible(NULL))
  }

  # Bars under the tree: one per partition, then sites and species. One
  # rectangle per run of neighbouring sites of the same cluster (a cluster can
  # be split across the tree with dynamic tree cuts); unassigned sites are
  # left blank.
  draw_rect <- function(p0, h0, p1, h1, ...) {
    a <- to_xy(p0, h0)
    b <- to_xy(p1, h1)
    graphics::rect(a$x, a$y, b$x, b$y, ...)
  }
  for (k in seq_len(n_bars)) {
    bar_top <- -gap_b - (k - 1) * bar_h
    bar_bottom <- bar_top - 0.85 * bar_h
    if (k <= n_part) {
      v <- as.character(clusters[tree$order, partitions[k]])
      col <- stats::setNames(colors[[partitions[k]]]$color,
                             colors[[partitions[k]]]$cluster)[v]
      name <- partitions[k]
    } else {
      v <- type[tree$order]
      col <- ifelse(v == "site", "grey25", "grey80")
      name <- "Type"
    }
    run_id <- ifelse(is.na(v), paste0("NA_", seq_len(n)), v)
    run_start <- which(c(TRUE, run_id[-1] != run_id[-n]))
    run_end <- c(run_start[-1] - 1, n)
    assigned <- !is.na(v[run_start])
    draw_rect(run_start[assigned] - 0.5, bar_bottom, run_end[assigned] + 0.5,
              bar_top, col = col[run_start[assigned]], border = "white",
              lwd = 0.6)
    draw_rect(0.5, bar_bottom, n + 0.5, bar_top, border = "grey70", lwd = 0.6)
    if (horizontal) {
      graphics::text((bar_top + bar_bottom) / 2, 0.5, name, srt = 90,
                     adj = c(-0.1, 0.5), xpd = NA, cex = 0.75)
    } else {
      graphics::text(0.5, (bar_top + bar_bottom) / 2, name, adj = c(1.1, 0.5),
                     xpd = NA, cex = 0.75)
    }
  }

  # Horizontal trees: the names of the partitions are above the bars, on the
  # right, so the legend goes on the left, next to the root
  if (bipartite) {
    position <- if (!horizontal) "topright" else if (flip) "bottomleft" else
      "topleft"
    graphics::legend(position, legend = c("sites", "species"),
                     fill = c("grey25", "grey80"), border = NA, bty = "n",
                     cex = 0.75, horiz = TRUE,
                     inset = if (horizontal) 0.01 else c(0, -0.08), xpd = NA)
  }
  if (labels) {
    if (horizontal) {
      plot_text(bars_bottom - 0.01 * max_h, seq_len(n),
                tree$labels[tree$order], ...,
                defaults = list(adj = c(0, 0.5), xpd = NA))
    } else {
      plot_text(seq_len(n), bars_bottom - 0.01 * max_h,
                tree$labels[tree$order], ...,
                defaults = list(srt = 90, adj = c(1, 0.5), xpd = NA))
    }
  }
}

# Cluster of each merge of a tree: the cluster of all the sites under the
# merge, or NA when they belong to several clusters (or have none)
merge_clusters <- function(tree, cl) {
  out <- rep(NA_character_, nrow(tree$merge))
  for (i in seq_len(nrow(tree$merge))) {
    child_cl <- vapply(tree$merge[i, ], function(k) {
      if (k < 0) cl[-k] else out[k]
    }, character(1))
    if (!anyNA(child_cl) && child_cl[1] == child_cl[2]) out[i] <- child_cl[1]
  }
  out
}

# TRUE if each cluster is a single branch of the tree, i.e. if each cluster
# has a single highest merge (or is a single site)
clusters_are_branches <- function(tree, cl, merge_cl) {
  parent <- integer(nrow(tree$merge))
  site_parent <- integer(length(cl))
  for (i in seq_len(nrow(tree$merge))) {
    for (k in tree$merge[i, ]) {
      if (k > 0) parent[k] <- i else site_parent[-k] <- i
    }
  }
  parent_cl <- function(p) ifelse(p == 0, NA_character_, merge_cl[pmax(p, 1)])
  top_merges <- !is.na(merge_cl) &
    (is.na(parent_cl(parent)) | parent_cl(parent) != merge_cl)
  top_sites <- !is.na(cl) &
    (is.na(parent_cl(site_parent)) | parent_cl(site_parent) != cl)
  all(table(c(merge_cl[top_merges], cl[top_sites])) == 1)
}

# Function to test if each cluster of a partition is included in
# a single cluster of the
# previous partition (with fewer clusters)
partitions_nested <- function(clusters, partitions) {
  for (k in seq_along(partitions)[-1]) {
    ok <- stats::complete.cases(clusters[, partitions[c(k - 1, k)]])
    n_parents <- tapply(clusters[ok, partitions[k - 1]],
                        clusters[ok, partitions[k]],
                        function(v) length(unique(v)))
    if (any(n_parents > 1)) return(FALSE)
  }
  TRUE
}

# Tree rebuilt from nested partitions (e.g. map equation results)
# with one step of height per partition:
# root > clusters of the first partition > ... > sites. A cluster split into
# n clusters at the next level becomes n - 1 binary merges at the same height,
# merged two by two so that the merge is drawn near the middle of its
# branches. Clusters are ordered by decreasing size, then by name.
nested_tree <- function(clusters, partitions) {
  L <- length(partitions)
  clusters <- clusters[stats::complete.cases(clusters[, partitions,
                                                      drop = FALSE]), ,
                       drop = FALSE]
  labels <- as.character(clusters[, 1])
  cl <- lapply(partitions, function(p) as.character(clusters[, p]))
  keys <- list()
  for (k in seq_len(L)) {
    size <- stats::ave(seq_along(labels), cl[[k]], FUN = length)
    name_rank <- match(cl[[k]], unique(cl[[k]])[natural_order(unique(
      cl[[k]]))])
    keys <- c(keys, list(-size, name_rank))
  }
  site_order <- do.call(order, c(keys, list(labels)))

  merge <- matrix(0L, 0, 2)
  height <- numeric(0)
  join <- function(refs, h) {
    while (length(refs) > 1) {
      pairs <- split(refs, ceiling(seq_along(refs) / 2))
      refs <- vapply(pairs, function(pr) {
        if (length(pr) == 1) return(pr)
        merge <<- rbind(merge, pr)
        height <<- c(height, h)
        nrow(merge)
      }, numeric(1))
    }
    refs
  }
  # Sites are referred to as -i, merges by their row in merge
  children <- -site_order
  for (k in L:1) {
    groups <- cl[[k]][site_order]
    if (k == L) {
      members <- split(children, factor(groups, unique(groups)))
    } else {
      sub <- cl[[k + 1]][site_order]
      first <- !duplicated(sub)
      members <- split(refs_k[sub[first]], factor(groups[first],
                                                   unique(groups[first])))
    }
    refs_k <- vapply(members, join, numeric(1), h = L - k + 1)
  }
  join(refs_k, L + 1)
  storage.mode(merge) <- "integer"
  structure(list(merge = unname(merge), height = height, order = site_order,
                 labels = labels, method = "nested partitions"),
            class = "hclust")
}


# Squares and circles ----------------------------------------------------------

# Number of sites per square and of species per circle: one, or a round number
# when there would be more than 1000 squares and circles
units_per_item <- function(items_per_unit, n_items) {
  if (identical(items_per_unit, "auto")) {
    round_numbers <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000)
    items_per_unit <- round_numbers[
      which(n_items / round_numbers <= 1000)[1]]
  }
  stats::setNames(rep_len(items_per_unit, 2), c("site", "species"))
}

units_text <- function(k_unit, bipartite) {
  txt <- paste0("1 square = ", k_unit[["site"]],
                if (k_unit[["site"]] > 1) " sites" else " site")
  if (bipartite) {
    txt <- paste0(txt, ", 1 circle = ", k_unit[["species"]], " species")
  }
  txt
}

units_legend <- function() {
  graphics::legend("topright", legend = c("sites", "species"),
                   pch = c(15, 16), col = "grey40", bty = "n", cex = 0.75,
                   horiz = TRUE, inset = c(0, -0.08), xpd = NA)
}

# Squares and circles of one cluster: sites first, then species, each ordered
# by the cluster of fill they belong to (largest first); one per k_unit items
cluster_units <- function(idx, type, fill, fill_colors, k_unit) {
  fill_size <- table(fill[idx])
  fill_names <- names(fill_size)
  fill_rank <- match(fill[idx],
                     fill_names[order(-fill_size, match(seq_along(fill_names),
                                                        natural_order(
                                                          fill_names)))])
  idx <- idx[order(type[idx] != "site", fill_rank)]
  out <- lapply(c(site = "site", species = "species"), function(t) {
    i <- idx[type[idx] == t]
    i <- i[(seq_along(i) - 1) %% k_unit[[t]] == 0]
    col <- unname(fill_colors[fill[i]])
    col[is.na(col)] <- "grey85"
    data.frame(col = col, shape = rep(t, length(i)))
  })
  out[vapply(out, nrow, numeric(1)) > 0]
}

# Coordinates of the squares of a bar w squares wide, with rows q long. The
# bar grows upwards (or rightwards if horizontal; dir = -1 reverses it) from
# start, and offset is its position across.
stack_xy <- function(n, w, offset, dir, start, horizontal, q = 1) {
  i <- seq_len(n) - 1
  across <- offset + i %% w + 0.5
  along <- dir * (start + (i %/% w + 0.5) * q)
  if (horizontal) list(x = along, y = -across) else list(x = across, y = along)
}

# Squares (sites) and circles (species) centred on (x, y), size[1] wide and
# size[2] high
draw_units <- function(x, y, col, shape, size = 0.82) {
  size <- rep_len(size, 2)
  sq <- shape == "site"
  if (any(sq)) {
    graphics::rect(x[sq] - size[1] / 2, y[sq] - size[2] / 2,
                   x[sq] + size[1] / 2, y[sq] + size[2] / 2, col = col[sq],
                   border = NA)
  }
  if (any(!sq)) {
    a <- seq(0, 2 * pi, length.out = 25)
    px <- rbind(outer(cos(a) * size[1] / 2, x[!sq], "+"), NA)
    py <- rbind(outer(sin(a) * size[2] / 2, y[!sq], "+"), NA)
    graphics::polygon(as.vector(px), as.vector(py), col = col[!sq],
                      border = NA)
  }
}

# Text with default arguments that users can change with ...
plot_text <- function(x, y, labels, ..., defaults = list()) {
  do.call(graphics::text, c(list(x = x, y = y, labels = labels),
                            utils::modifyList(defaults, list(...))))
}
