# Inputs -----------------------------------------------------------------------
dissim <- dissimilarity(fishmat, metric = "all")

set.seed(1)
tree1 <- hclu_hierarclust(dissim, 
                          n_clust = 2:48, 
                          index = "Simpson",
                          optimal_tree_method = "best")

a <- bioregionalization_metrics(tree1,
                                dissimilarity = dissim,
                                net = fishdf,
                                species_col = "Species",
                                site_col = "Site",
                                eval_metric = c("tot_endemism",
                                                "avg_endemism",
                                                "pc_distance",
                                                "anosim"))

a1 <- bioregionalization_metrics(tree1,
                                dissimilarity = dissim,
                                net = fishdf,
                                species_col = "Species",
                                site_col = "Site",
                                eval_metric = c("tot_endemism",
                                                "avg_endemism",
                                                "pc_distance"))

partdf <- data.frame(K = c(3,5,6,2,3), 
                     n_clusters = c(3,5,6,23,6),
                     metric1 = c(NA,0.2,0.6,0.7,13),
                     metric2 = c("1",0.2,0.7,12,8))

partdf1 <- data.frame(K = c(3,5), 
                      n_cluster = c(3,5),
                      metric1 = c(NA,0.2),
                      metric2 = c(1,0.2))

partdf2 <- data.frame(K = c(3,5), 
                   n_clusters = c(3,5),
                   metric1 = c(NA,0.2),
                   metric2 = c(1,0.2))

partdf3 <- data.frame(K = c(3,5,6,2,3), 
                     n_clusters = c(3,5,6,23,6),
                     metric1 = c(NA,0.2,0.6,0.7,13),
                     metric2 = rep(1,5))

partdf4 <- data.frame(K = c(3,5,6,2,3), 
                      n_clusters = c(3,5,6,23,6),
                      metric1 = rep(1,5),
                      metric2 = rep(1,5))

a2 <- a
a2$evaluation_df$tot_endemism <- c(2, rep(1, dim(a2$evaluation_df)[1]-1))

a3 <- a
a3$evaluation_df <- a3$evaluation_df[1:8,]



# Tests for valid outputs ------------------------------------------------------
test_that("valid output", {
  
  optim_a <- find_optimal_n(a1, 
                            plot = FALSE)
  expect_equal(optim_a$optimal_nb_clusters$tot_endemism, 13)
  expect_equal(optim_a$optimal_nb_clusters$avg_endemism, 12)
  expect_equal(optim_a$optimal_nb_clusters$pc_distance, 13)
  
  #optim_a <- find_optimal_n(a1, 
  #                          plot = TRUE)
  #expect_equal(optim_a$optimal_nb_clusters$tot_endemism, 13)
  #expect_equal(optim_a$optimal_nb_clusters$avg_endemism, 12)
  #expect_equal(optim_a$optimal_nb_clusters$pc_distance, 13)
  
  #optim_a <- find_optimal_n(a1,
  #                          criterion = "breakpoints",
  #                          plot = TRUE)
  #expect_equal(optim_a$optimal_nb_clusters$tot_endemism, 11)
  #expect_equal(optim_a$optimal_nb_clusters$avg_endemism, 10)
  #expect_equal(optim_a$optimal_nb_clusters$pc_distance, 13)
  
  #optim_a <- find_optimal_n(a1, 
  #                          criterion = "max",
  #                          plot = FALSE)
  
  #optim_a <- find_optimal_n(a1, 
  #                          criterion = "min",
  #                          plot = FALSE)

})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    find_optimal_n("zz"),
    "^bioregionalizations should be the output object from")
  
  expect_error(
    find_optimal_n(partdf),
    "^Your bioregionalization data.frame contains non numeric")
  
  expect_error(
    find_optimal_n(partdf1),
    "^bioregionalizations should be the output object from")
  
  expect_error(
    find_optimal_n(partdf2),
    "^The number of bioregionalizations is too low")
  
  expect_error(
    find_optimal_n(a,
                   criterion =  c("zz","zz")),
    "criterion must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    find_optimal_n(a,
                   criterion =  TRUE),
    "criterion must be a character.",
    fixed = TRUE)
  
  expect_error(
    find_optimal_n(a,
                   criterion =  "zzb"),
    "^Please choose criterion from the following:")
  
  expect_error(
    find_optimal_n(a,
                   metrics_to_use =  c(1,TRUE)),
    "metrics_to_use must be a character.",
    fixed = TRUE)
  
  expect_error(
    find_optimal_n(a, 
                   metrics_to_use = "zz"),
    "metrics_to_use should exist in the evaluation table.",
    fixed = TRUE)
  
  expect_warning(
    find_optimal_n(partdf3,
                   plot = FALSE),
    "^Metrics metric2")
  
  expect_warning(
    expect_error(
      find_optimal_n(partdf4,
                     plot = FALSE),
      "^The selected bioregionalization metrics"),
    "^Metrics metric1")
  
  expect_warning(
   find_optimal_n(a2,
                  criterion = "breakpoints",
                  plot = FALSE),
   "^Metrics tot_endemism")
  
  expect_message(
    expect_warning(
      find_optimal_n(a3,
                     criterion = "increasing_step",
                     plot = FALSE),
      "^Criterion 'increasing_step' cannot work properly with "),
   "^...Caveat: be cautious with the interpretation of")
  
  expect_warning(
    find_optimal_n(a3,
                   criterion = "decreasing_step",
                   plot = FALSE),
    "^Criterion 'decreasing_step' cannot work properly with")
  
  expect_error(
    find_optimal_n(a,
                   criterion = "increasing_step",
                   step_quantile = c("zz","zz")),
    "step_quantile must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    find_optimal_n(a,
                   criterion = "increasing_step",
                   step_quantile = "zz"),
    "step_quantile must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    find_optimal_n(a,
                   criterion = "increasing_step",
                   step_quantile = -1.1),
    "step_quantile must be strictly higher than 0.",
    fixed = TRUE)  
  
  expect_error(
    find_optimal_n(a,
                   criterion = "increasing_step",
                   step_quantile = 0),
    "step_quantile must be strictly higher than 0.",
    fixed = TRUE)  
  
  expect_error(
    find_optimal_n(a,
                   criterion = "increasing_step",
                   step_quantile = 1),
    "step_quantile must be in the ]0,1[ interval.",
    fixed = TRUE)  
  
  expect_error(
    find_optimal_n(a,
                   criterion = "increasing_step",
                   step_round_above = c(TRUE, TRUE)),
    "step_round_above must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    find_optimal_n(a,
                   criterion = "increasing_step",
                   step_round_above = 1),
    "step_round_above must be a boolean.",
     fixed = TRUE)  
  
  expect_error(
    find_optimal_n(a,
                   criterion = "increasing_step",
                   step_levels = c(1,1)),
    "step_levels must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    find_optimal_n(a,
                   criterion = "increasing_step",
                   step_levels = "zz"),
    "step_levels must be numeric.",
    fixed = TRUE)
    
  expect_error(
    find_optimal_n(a,
                   criterion = "increasing_step",
                   step_levels = 1.1),
    "step_levels must be an integer.",
    fixed = TRUE)
  
  expect_error(
    find_optimal_n(a,
                   criterion = "increasing_step",
                   step_levels = -1),
    "step_levels must be higher than 0.",
    fixed = TRUE)
  
  expect_error(
    find_optimal_n(a, 
                   plot = c(TRUE, TRUE)),
    "plot must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    find_optimal_n(a, 
                   plot = "zz"),
    "plot must be a boolean.",
    fixed = TRUE)
  
  expect_warning(
    find_optimal_n(a,
                   criterion = "elbow",
                   plot = FALSE),
    "^The elbow method")
  
  expect_error(
    find_optimal_n(a,
                   criterion = "cutoff",
                   plot = FALSE),
    "^Criterion 'cutoff' should probably")
  
})

# metric_cutoffs = c(.5, .75, .9, .95, .99, .999) # not tested but error message
# exists
