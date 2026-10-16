# Workflow tests
#
# Unlike the other test files, this file does not mirror a source file: it
# runs the analyses the way users run them, from the site-species matrix to
# the final outputs, and checks every hand-off between two functions. Each step
# takes the real output of the previous function, never a hand-made object.
#
# The shared data (wf_comat, wf_dissim, wf_sim) and the checks
# (expect_valid_bioregion_clusters(), expect_output_contains(), ...) are in
# helper-workflow.R.

# Inputs -----------------------------------------------------------------------
wf_tree <- hclu_hierarclust(wf_dissim,
                            n_clust = 2:15,
                            n_runs = 10,
                            seed = 1,
                            verbose = FALSE)
wf_eval <- bioregionalization_metrics(wf_tree,
                                      eval_metrics = "all",
                                      dissimilarity = wf_dissim,
                                      comat = wf_comat)
wf_kmeans <- nhclu_kmeans(wf_dissim, n_clust = 2:8, seed = 1)
wf_pam <- nhclu_pam(wf_dissim, n_clust = 2:8, seed = 1)
wf_greedy <- netclu_greedy(wf_sim)
wf_louvain <- netclu_louvain(wf_sim, seed = 1)
wf_net <- mat_to_net(wf_comat, weight = TRUE)
wf_greedy_bip <- netclu_greedy(wf_net, bipartite = TRUE)
wf_louvain_bip <- netclu_louvain(wf_net, bipartite = TRUE, seed = 1)

all_metrics <- c("prop_between_dissim", "anosim", "mean_endemics",
                 "tot_endemics")

# Checks the output of bioregionalization_metrics() against the
# bioregionalization it evaluates
expect_valid_evaluation <- function(ev, bioregionalization) {
  expect_s3_class(ev, "bioregion.bioregionalization.metrics")
  expect_equal(names(ev), c("partition", "n_bioregions", all_metrics))
  expect_equal(ev$partition, bioregionalization$cluster_info$partition_name)
  expect_equal(ev$n_bioregions,
               unname(bioregionalization$cluster_info$n_clust))
}

# Checks the output of find_optimal_n() and its print()
expect_valid_optimal_n <- function(opt, ev, metrics, criterion) {
  expect_s3_class(opt, "bioregion.optimal.n")
  expect_equal(names(opt$optimal_n), metrics)
  # The optimal numbers of clusters are numbers of clusters that were
  # evaluated
  found <- unlist(opt$optimal_n)
  expect_true(all(found[!is.na(found)] %in% ev$n_bioregions))
  expect_output_contains(
    print(opt),
    c(paste("Evaluated metric(s): ", paste(metrics, collapse = " ")),
      paste("Criterion chosen to optimise the number of clusters: ",
            criterion)))
  expect_output_contains(print(opt), optimal_n_lines(opt), whole_line = TRUE)
}


# Workflow 1: hierarchical clustering -------------------------------------------
test_that("workflow 1: hierarchical clustering", {

  expect_valid_bioregion_clusters(wf_tree)
  expect_valid_evaluation(wf_eval, wf_tree)

  # find_optimal_n() with each criterion
  criteria <- c("elbow", "increasing_step", "decreasing_step", "cutoff",
                "breakpoints", "min", "max")
  for (criterion in criteria) {
    # cutoff only works with one metric
    metrics <- if (criterion == "cutoff") "prop_between_dissim" else all_metrics
    opt <- find_optimal_n(wf_eval,
                          metrics_to_use = metrics,
                          criterion = criterion,
                          plot = FALSE,
                          verbose = FALSE)
    expect_valid_optimal_n(opt, wf_eval, metrics, criterion)
  }

  # The output of find_optimal_n() as n_clust
  opt <- find_optimal_n(wf_eval, plot = FALSE, verbose = FALSE)
  n_opt <- sort(unique(unlist(opt$optimal_n)))

  cut_opt <- cut_tree(wf_tree, n_clust = opt, verbose = FALSE)
  cut_num <- cut_tree(wf_tree, n_clust = n_opt, verbose = FALSE)
  expect_equal(cut_opt$clusters, cut_num$clusters)
  expect_equal(cut_opt$cluster_info, cut_num$cluster_info)
  expect_equal(cut_opt$cluster_info$requested_n_clust, n_opt)
  expect_valid_bioregion_clusters(cut_opt)

  tree_opt <- hclu_hierarclust(wf_dissim,
                               n_clust = opt,
                               n_runs = 10,
                               seed = 1,
                               verbose = FALSE)
  expect_equal(tree_opt$clusters, cut_num$clusters)
  expect_equal(tree_opt$cluster_info, cut_num$cluster_info)
  expect_valid_bioregion_clusters(tree_opt)

  # Cuts at given heights
  cut_h <- cut_tree(wf_tree, cut_height = c(0.3, 0.6, 0.9), verbose = FALSE)
  expect_valid_bioregion_clusters(cut_h)
})


# Workflow 2: DIANA -------------------------------------------------------------
test_that("workflow 2: DIANA", {

  # DIANA trees of Simpson dissimilarities often have many branches merging at
  # height 1 (sites sharing no species), so several requested numbers of
  # clusters cannot be reached. Each of them must give a warning naming the
  # number of clusters actually returned.
  requested <- 2:15
  warns <- character(0)
  dia <- withCallingHandlers(
    hclu_diana(wf_dissim, n_clust = requested, verbose = FALSE),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_valid_bioregion_clusters(dia)
  # This dataset does have unreachable numbers of clusters with DIANA (if this
  # fails, the rest of the test no longer checks the warnings)
  expect_gt(length(warns), 0)

  # The partitions are in the order of the requested numbers of clusters:
  # each one has the requested number, or the number given in its warning
  pattern <- "^.*for k = ([0-9]+)\\. Closest number found: ([0-9]+)$"
  expect_true(all(grepl(pattern, warns)))
  warned_k <- as.numeric(sub(pattern, "\\1", warns))
  warned_found <- as.numeric(sub(pattern, "\\2", warns))
  expected_n <- requested
  expected_n[match(warned_k, requested)] <- warned_found
  expect_equal(dia$cluster_info$n_clust, expected_n)

  ev <- bioregionalization_metrics(dia,
                                   eval_metrics = "all",
                                   dissimilarity = wf_dissim,
                                   comat = wf_comat)
  expect_valid_evaluation(ev, dia)

  opt <- find_optimal_n(ev, plot = FALSE, verbose = FALSE)
  expect_valid_optimal_n(opt, ev, all_metrics, "elbow")

  # The output of find_optimal_n() as n_clust
  n_opt <- sort(unique(unlist(opt$optimal_n)))
  dia_opt <- hclu_diana(wf_dissim, n_clust = opt, verbose = FALSE)
  dia_num <- hclu_diana(wf_dissim, n_clust = n_opt, verbose = FALSE)
  expect_equal(dia_opt$clusters, dia_num$clusters)
  expect_equal(dia_opt$cluster_info, dia_num$cluster_info)
  expect_valid_bioregion_clusters(dia_opt)
})


# Workflow 3: non-hierarchical clustering ---------------------------------------
# The metrics and the search for an optimal number of clusters on partitions
# that are not nested
test_that("workflow 3: non-hierarchical clustering", {

  for (bioregionalization in list(wf_kmeans, wf_pam)) {
    expect_valid_bioregion_clusters(bioregionalization)
    expect_output_contains(print(bioregionalization),
                           "Partitions are not hierarchical")

    ev <- bioregionalization_metrics(bioregionalization,
                                     eval_metrics = "all",
                                     dissimilarity = wf_dissim,
                                     comat = wf_comat)
    expect_valid_evaluation(ev, bioregionalization)

    opt <- find_optimal_n(ev, plot = FALSE, verbose = FALSE)
    expect_valid_optimal_n(opt, ev, all_metrics, "elbow")
  }
})


# Workflow 4: network clustering ------------------------------------------------
test_that("workflow 4: network clustering", {

  # Unipartite networks of sites, and bipartite site-species networks
  for (bioregionalization in list(wf_greedy, wf_louvain,
                                  wf_greedy_bip, wf_louvain_bip)) {
    expect_valid_bioregion_clusters(bioregionalization)

    orderings <- "n_sites"
    if (bioregionalization$inputs$bipartite) {
      orderings <- c(orderings, "n_species", "n_both")
    }
    for (ordering in orderings) {
      colored <- bioregion_colors(bioregionalization,
                                  cluster_ordering = ordering)
      expect_valid_bioregion_clusters(colored)

      # Every cluster of every partition has a color
      partitions <- names(colored$clusters)[-1]
      expect_equal(names(colored$colors), partitions)
      for (p in partitions) {
        clusters <- unique(stats::na.omit(colored$clusters[[p]]))
        expect_true(all(clusters %in% colored$colors[[p]]$cluster))
      }
      expect_equal(colored$clusters_colors[[1]], colored$clusters[[1]])
      expect_output_contains(print(colored), "Color palette assigned")
    }
  }

  # Bipartite networks: both sites and species are clustered
  expect_setequal(attr(wf_greedy_bip$clusters, "node_type"),
                  c("site", "species"))
  expect_equal(sum(attr(wf_greedy_bip$clusters, "node_type") == "site"),
               nrow(wf_comat))
})


# Workflow 5: everything that uses a bioregionalization -------------------------
test_that("workflow 5: site_species_metrics()", {

  # One output with several partitions, one with a single partition and a site
  # assigned to no cluster (a site sharing no species with the others is
  # absent from the similarity network), and a bipartite one
  expect_true(anyNA(wf_greedy$clusters[[2]]))
  ssm_inputs <- list(
    tree = list(x = wf_tree, cluster_on = "site", similarity = wf_sim),
    greedy = list(x = wf_greedy, cluster_on = "site", similarity = wf_sim),
    bipartite = list(x = wf_greedy_bip, cluster_on = "both",
                     similarity = NULL))

  for (input in ssm_inputs) {
    ssm <- suppressWarnings(site_species_metrics(
      input$x,
      bioregion_metrics = "all",
      bioregionalization_metrics = "all",
      cluster_on = input$cluster_on,
      comat = wf_comat,
      similarity = input$similarity,
      verbose = FALSE))
    expect_s3_class(ssm, "bioregion.site.species.metrics")

    # print(), str() and summary() list every computed metric
    metrics <- ssm_metric_names(ssm)
    expect_output_contains(print(ssm), metrics)
    expect_output_contains(str(ssm), metrics)
    expect_output_contains(summary(ssm), ssm_metric_names(ssm, columns = TRUE))

    # summary() gives the number of bioregions (and chorotypes) of the first
    # partition, and the names of the top species
    cl <- input$x$clusters
    node_type <- attr(cl, "node_type")
    n_bioregions <- length(unique(stats::na.omit(cl[node_type == "site", 2])))
    clust_str <- paste0("(", n_bioregions, " bioregions")
    if (input$cluster_on == "both") {
      n_chorotypes <- length(unique(stats::na.omit(
        cl[node_type == "species", 2])))
      clust_str <- paste0(clust_str, ", ", n_chorotypes, " chorotypes)")
    } else {
      clust_str <- paste0(clust_str, ")")
    }
    if (is.data.frame(ssm[[1]])) {
      first <- ssm
      first_name <- "Single partition"
    } else {
      first <- ssm[[1]]
      first_name <- names(ssm)[1]
    }
    sb <- first$species_bioregions
    top_species <- sb$species[order(-sb$ind_val_occ)][1]
    expect_output_contains(
      summary(ssm),
      c(paste("Partition 1:", first_name, clust_str), top_species))
  }

  # A site assigned to no cluster gets no silhouette
  ssm <- site_species_metrics(wf_greedy,
                              bioregion_metrics = "all",
                              bioregionalization_metrics = "all",
                              comat = wf_comat,
                              similarity = wf_sim,
                              verbose = FALSE)
  unassigned <- wf_greedy$clusters[is.na(wf_greedy$clusters[[2]]), 1]
  sil <- ssm$site_bioregionalization
  expect_true(all(is.na(sil$silhouette[sil$site %in% unassigned])))
  expect_false(anyNA(sil$silhouette[!(sil$site %in% unassigned)]))
})

test_that("workflow 5: bioregion_metrics()", {

  for (bioregionalization in list(wf_tree, wf_greedy)) {
    bm <- bioregion_metrics(bioregionalization, comat = wf_comat)
    partitions <- names(bioregionalization$clusters)[-1]
    if (length(partitions) == 1) {
      bm <- list(bm)
      names(bm) <- partitions
    }
    expect_equal(names(bm), partitions)
    for (p in partitions) {
      clusters <- bioregionalization$clusters[[p]]
      expect_setequal(bm[[p]]$bioregion, unique(stats::na.omit(clusters)))
      expect_equal(sum(bm[[p]]$n_sites), sum(!is.na(clusters)))
    }
  }
})

test_that("workflow 5: compare_bioregionalizations()", {

  cmp <- compare_bioregionalizations(wf_kmeans,
                                     cor_frequency = TRUE,
                                     verbose = FALSE)
  partitions <- names(wf_kmeans$clusters)[-1]
  n_partitions <- length(partitions)
  expect_s3_class(cmp, "bioregion.bioregionalization.comparison")
  expect_equal(unname(cmp$inputs["number_items"]), nrow(wf_kmeans$clusters))
  expect_equal(unname(cmp$inputs["number_bioregionalizations"]), n_partitions)
  expect_equal(colnames(cmp$pairwise_membership), partitions)
  expect_equal(names(cmp$bioregionalization_freq_cor), partitions)
  expect_equal(nrow(cmp$bioregionalization_comparison),
               choose(n_partitions, 2))

  most_representative <- names(cmp$bioregionalization_freq_cor)[
    which.max(cmp$bioregionalization_freq_cor)]
  expect_output_contains(
    print(cmp),
    c(paste(n_partitions, "bioregionalizations compared"),
      paste(nrow(wf_kmeans$clusters), "items in the clustering"),
      "Requested indices: rand jaccard",
      most_representative))

  # Same result as with the data.frame of the partitions
  cmp_df <- compare_bioregionalizations(wf_kmeans$clusters[, -1],
                                        cor_frequency = TRUE,
                                        verbose = FALSE)
  expect_equal(cmp_df$bioregionalization_comparison,
               cmp$bioregionalization_comparison)

  # One index only
  cmp_rand <- compare_bioregionalizations(wf_kmeans,
                                          indices = "rand",
                                          verbose = FALSE)
  expect_output_contains(print(cmp_rand), "Requested indices: rand")
})

test_that("workflow 5: map_bioregions()", {

  map <- fishsf[fishsf$Site %in% rownames(wf_comat), ]
  out <- map_bioregions(wf_tree,
                        map = map,
                        plot = FALSE,
                        map_as_output = TRUE)
  partitions <- names(wf_tree$clusters)[-1]
  expect_s3_class(out, "sf")
  expect_equal(nrow(out), nrow(wf_comat))
  expect_true(all(partitions %in% names(out)))
  # Each site keeps its clusters
  idx <- match(out$Site, wf_tree$clusters[[1]])
  for (p in partitions) {
    expect_equal(as.character(out[[p]]),
                 as.character(wf_tree$clusters[[p]][idx]))
  }
})


# Workflow 6: old names still accepted ------------------------------------------
# NEWS states that the old metric names are still accepted as inputs
test_that("workflow 6: old metric names are still accepted", {

  ev_old <- bioregionalization_metrics(
    wf_tree,
    eval_metrics = c("pc_distance", "avg_endemism", "tot_endemism"),
    dissimilarity = wf_dissim,
    comat = wf_comat)
  ev_new <- bioregionalization_metrics(
    wf_tree,
    eval_metrics = c("prop_between_dissim", "mean_endemics", "tot_endemics"),
    dissimilarity = wf_dissim,
    comat = wf_comat)
  expect_equal(ev_old, ev_new)

  opt_old <- find_optimal_n(wf_eval,
                            metrics_to_use = "pc_distance",
                            plot = FALSE,
                            verbose = FALSE)
  opt_new <- find_optimal_n(wf_eval,
                            metrics_to_use = "prop_between_dissim",
                            plot = FALSE,
                            verbose = FALSE)
  expect_equal(opt_old, opt_new)

  # A user's own data.frame with a column literally named pc_distance keeps
  # this name: the conversion only applies to names absent from the columns
  user_df <- data.frame(partition = wf_eval$partition,
                        n = wf_eval$n_bioregions,
                        pc_distance = wf_eval$prop_between_dissim)
  opt_user <- find_optimal_n(user_df,
                             metrics_to_use = "pc_distance",
                             plot = FALSE,
                             verbose = FALSE)
  expect_equal(names(opt_user$optimal_n), "pc_distance")
  expect_equal(unname(opt_user$optimal_n), unname(opt_new$optimal_n))
  expect_output_contains(print(opt_user), optimal_n_lines(opt_user),
                         whole_line = TRUE)
})
