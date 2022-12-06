# Preamble code ----------------------------------------------------------------
dissim <- dissimilarity(vegemat, metric = "all")

# User-defined number of clusters
tree1 <- hclu_hierarclust(dissim, n_clust = 2:50, index = "Simpson")

a <- partition_metrics(tree1,
                       dissimilarity = dissim,
                       net = vegedf,
                       species_col = "Species",
                       site_col = "Site",
                       eval_metric = c("tot_endemism",
                                       "avg_endemism",
                                       "pc_distance",
                                       "anosim"))

# Tests for valid outputs -----------------------------------------------------
optim_a <- find_optimal_n(a, plot = FALSE)

test_that("number of columns in output", {
  
  expect_equal(optim_a$optimal_nb_clusters$tot_endemism, 8L)
  expect_equal(optim_a$optimal_nb_clusters$avg_endemism, 8L)
  expect_equal(optim_a$optimal_nb_clusters$pc_distance, 14L)
  expect_equal(optim_a$optimal_nb_clusters$anosim, 17L)
  
})

# Tests for invalid inputs ----------------------------------------------------
test_that("error messages with wrong inputs", {
  expect_error(
    find_optimal_n("zz"),
    "partitions should be the output object from partition_metrics()or a data.frame",
    fixed = TRUE)
  
  expect_error(
    find_optimal_n(a, metrics_to_use = "zz"),
    "metrics_to_use should exist in the evaluation table",
    fixed = TRUE)
  
  expect_error(
    find_optimal_n(a, criterion = "zz"),
    "criterion must be one of elbow, increasing_step, decreasing_step,
           min, max, cutoff or mars",
    fixed = TRUE)
  
  expect_error(
    find_optimal_n(a, criterion = "increasing_step", step_quantile = "zz"),
    "step_quantile must be a numeric in the ]0,1[ interval. See help of
           the function.",
    fixed = TRUE)
  
  expect_error(
    find_optimal_n(a, criterion = "increasing_step", step_levels = "zz"),
    "step_levels must be a positive integer. See help of the
             function.",
    fixed = TRUE)
  
  expect_error(
    find_optimal_n(a, criterion = "increasing_step", step_round_above = "zz"),
    "step_round_above must be a boolean. See help of the function.",
    fixed = TRUE)
  
  expect_error(
    find_optimal_n(a, criterion = "mars", mars_breakpoints = "zz"),
    "mars_breakpoints should be 'all' or 'increasing' or 'decreasing'",
    fixed = TRUE)
  
  expect_error(
    find_optimal_n(a, plot = "zz"),
    "plot should be a Boolean.",
    fixed = TRUE)
})

# metric_cutoffs = c(.5, .75, .9, .95, .99, .999) # not tested but error message
# exists
