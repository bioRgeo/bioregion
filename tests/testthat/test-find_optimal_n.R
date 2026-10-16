# Inputs -----------------------------------------------------------------------
dissim <- dissimilarity(fishmat, metric = "all")

# 43 clusters cannot be reached on this tree (two merges at 1/7 up to a
# rounding error), hence the warning of cut_tree()
bioreg <- suppressWarnings(hclu_hierarclust(dissim, 
                                            n_clust = 2:48, 
                                            index = "Simpson",
                                            seed = 1,
                                            optimal_tree_method = "best",
                                            verbose = FALSE))

met <- bioregionalization_metrics(bioreg,
                                  eval_metrics = "all",
                                  dissimilarity = dissim,
                                  comat = fishmat)

met1 <- bioregionalization_metrics(bioreg,
                                   eval_metrics = c("mean_endemics",
                                                    "tot_endemics",
                                                    "prop_between_dissim"),
                                   dissimilarity = dissim,
                                   comat = fishmat)

partdf <- data.frame(K = c(3,5,6,2,3),
                     n_clusters = c(3,5,6,23,6),
                     metric1 = c(NA,0.2,0.6,0.7,13),
                     metric2 = c("1",0.2,0.7,12,8))

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

met2 <- met
met2$tot_endemics <- c(2, rep(1, dim(met2)[1]-1))

met3 <- met
met3 <- met3[1:8,]

# Tests for valid outputs ------------------------------------------------------
test_that("valid output", {
  
  optim_a <- find_optimal_n(met,
                            plot = FALSE,
                            verbose = FALSE)
  expect_equal(optim_a$optimal_n$tot_endemics, 13)
  expect_equal(optim_a$optimal_n$mean_endemics, 12)
  expect_equal(optim_a$optimal_n$prop_between_dissim, 13)

  # optim_a <- find_optimal_n(met1,
  #                           plot = TRUE,
  #                           verbose = FALSE)
  # expect_equal(optim_a$optimal_n$tot_endemics, 12)
  # expect_equal(optim_a$optimal_n$mean_endemics, 13)
  # expect_equal(optim_a$optimal_n$prop_between_dissim, 13)
  # 
  # optim_a <- find_optimal_n(met1,
  #                           criterion = "breakpoints",
  #                           plot = TRUE,
  #                           verbose = FALSE)
  # expect_equal(optim_a$optimal_n$tot_endemics, 10)
  # expect_equal(optim_a$optimal_n$mean_endemics, 9)
  # expect_equal(optim_a$optimal_n$prop_between_dissim, 11)

  # optim_a <- find_optimal_n(met1,
  #                           criterion = "max",
  #                           plot = FALSE,
  #                           verbose = FALSE)

  # optim_a <- find_optimal_n(met1,
  #                           criterion = "min",
  #                           plot = FALSE,
  #                           verbose = FALSE)

})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    find_optimal_n("zz"),
    "^evaluation_df should be the output object from")

  expect_error(
    find_optimal_n(partdf),
    "^evaluation_df contains non numeric")

  expect_error(
    find_optimal_n(partdf2),
    "^The number of partitions is too low")

  expect_error(
    find_optimal_n(met,
                   criterion =  c("zz","zz")),
    "criterion must be of length 1.",
    fixed = TRUE)

  expect_error(
    find_optimal_n(met,
                   criterion =  TRUE),
    "criterion must be a character.",
    fixed = TRUE)

  expect_error(
    find_optimal_n(met,
                   criterion =  "zzb"),
    "^Please choose criterion from the following:")

  expect_error(
    find_optimal_n(met,
                   metrics_to_use =  c(1,TRUE)),
    "metrics_to_use must be a character.",
    fixed = TRUE)

  expect_error(
    find_optimal_n(met,
                   metrics_to_use = "zz"),
    "metrics_to_use should exist in evaluation_df.",
    fixed = TRUE)

  expect_warning(
    find_optimal_n(partdf3,
                   plot = FALSE,
                   verbose = FALSE),
    "^Metrics metric2")

  expect_warning(
    expect_error(
      find_optimal_n(partdf4,
                     plot = FALSE),
      "^The selected evaluation metrics"),
    "^Metrics metric1")

  expect_warning(
   find_optimal_n(met2,
                  criterion = "breakpoints",
                  plot = FALSE,
                  verbose = FALSE),
   "^Metrics tot_endemics")

  quietly(
    expect_message(
      expect_warning(
        find_optimal_n(met3,
                       criterion = "increasing_step",
                       plot = FALSE),
        "^Criterion 'increasing_step' cannot work properly with "),
     "^...Caveat: be cautious with the interpretation of")
  )

  quietly(
    expect_warning(
      find_optimal_n(met3,
                     criterion = "decreasing_step",
                     plot = FALSE),
      "^Criterion 'decreasing_step' cannot work properly with")
  )

  expect_error(
    find_optimal_n(met,
                   criterion = "increasing_step",
                   step_quantile = c("zz","zz"),
                   verbose = FALSE),
    "step_quantile must be of length 1.",
    fixed = TRUE)

  expect_error(
    find_optimal_n(met,
                   criterion = "increasing_step",
                   step_quantile = "zz",
                   verbose = FALSE),
    "step_quantile must be numeric.",
    fixed = TRUE)

  expect_error(
    find_optimal_n(met,
                   criterion = "increasing_step",
                   step_quantile = -1.1,
                   verbose = FALSE),
    "step_quantile must be strictly higher than 0.",
    fixed = TRUE)

  expect_error(
    find_optimal_n(met,
                   criterion = "increasing_step",
                   step_quantile = 0,
                   verbose = FALSE),
    "step_quantile must be strictly higher than 0.",
    fixed = TRUE)

  expect_error(
    find_optimal_n(met,
                   criterion = "increasing_step",
                   step_quantile = 1,
                   verbose = FALSE),
    "step_quantile must be in the ]0,1[ interval.",
    fixed = TRUE)

  expect_error(
    find_optimal_n(met,
                   criterion = "increasing_step",
                   step_round_above = c(TRUE, TRUE),
                   verbose = FALSE),
    "step_round_above must be of length 1.",
    fixed = TRUE)

  expect_error(
    find_optimal_n(met,
                   criterion = "increasing_step",
                   step_round_above = 1,
                   verbose = FALSE),
    "step_round_above must be a boolean.",
     fixed = TRUE)

  expect_error(
    find_optimal_n(met,
                   criterion = "increasing_step",
                   step_levels = c(1,1),
                   verbose = FALSE),
    "step_levels must be of length 1.",
    fixed = TRUE)

  expect_error(
    find_optimal_n(met,
                   criterion = "increasing_step",
                   step_levels = "zz",
                   verbose = FALSE),
    "step_levels must be numeric.",
    fixed = TRUE)

  expect_error(
    find_optimal_n(met,
                   criterion = "increasing_step",
                   step_levels = 1.1,
                   verbose = FALSE),
    "step_levels must be an integer.",
    fixed = TRUE)

  expect_error(
    find_optimal_n(met,
                   criterion = "increasing_step",
                   step_levels = -1,
                   verbose = FALSE),
    "step_levels must be higher than 0.",
    fixed = TRUE)

  expect_error(
    find_optimal_n(met,
                   plot = c(TRUE, TRUE)),
    "plot must be of length 1.",
    fixed = TRUE)

  expect_error(
    find_optimal_n(met,
                   plot = "zz"),
    "plot must be a boolean.",
    fixed = TRUE)

  quietly(
    expect_warning(
      find_optimal_n(met,
                     criterion = "elbow",
                     plot = FALSE),
      "^The elbow method")
  )

  quietly(
    expect_error(
      find_optimal_n(met,
                     criterion = "cutoff",
                     plot = FALSE),
      "^Criterion 'cutoff' should probably")
  )

  expect_error(
    find_optimal_n(met,
                   plot = FALSE,
                   verbose = 1),
    "verbose must be a boolean.",
    fixed = TRUE)

  expect_error(
    find_optimal_n(met,
                   plot = FALSE,
                   verbose = c(TRUE, FALSE)),
    "verbose must be of length 1.",
    fixed = TRUE)
  
})

# metric_cutoffs = c(.5, .75, .9, .95, .99, .999) # not tested but error message
# exists

# Tests for the min, max and cutoff criteria -----------------------------------
test_that("criteria min and max", {

  for (criterion in c("min", "max")) {
    fun <- match.fun(criterion)
    opt <- find_optimal_n(met,
                          criterion = criterion,
                          plot = FALSE,
                          verbose = FALSE)
    metrics <- colnames(met)[-c(1, 2)]
    expect_equal(names(opt$optimal_n), metrics)
    for (m in metrics) {
      is_optimal <- met[[m]] == fun(met[[m]])
      expect_equal(opt$optimal_n[[m]], unique(met$n_bioregions[is_optimal]),
                   info = paste(criterion, m))
      # Only the optimal partitions of this metric are flagged
      expect_equal(opt$evaluation_df[[paste0("optimal_n_", m)]], is_optimal,
                   info = paste(criterion, m))
    }
  }

  # Nothing printed when verbose = FALSE
  expect_silent(find_optimal_n(met, criterion = "max",
                               plot = FALSE, verbose = FALSE))
})

test_that("criterion cutoff", {

  vals <- met$prop_between_dissim
  cutoffs <- c(0.5, 0.75, 0.9)
  opt <- find_optimal_n(met,
                        criterion = "cutoff",
                        metrics_to_use = "prop_between_dissim",
                        metric_cutoffs = cutoffs,
                        plot = FALSE,
                        verbose = FALSE)
  expected <- met$n_bioregions[sapply(cutoffs, function(x) which(vals >= x)[1])]
  expect_equal(opt$optimal_n$prop_between_dissim, expected)

  # A cutoff never reached by the metric gives no partition
  opt <- find_optimal_n(met,
                        criterion = "cutoff",
                        metrics_to_use = "prop_between_dissim",
                        metric_cutoffs = c(0.5, max(vals) + 1),
                        plot = FALSE,
                        verbose = FALSE)
  expect_equal(opt$optimal_n$prop_between_dissim,
               met$n_bioregions[which(vals >= 0.5)[1]])

  # No cutoff reached
  opt <- find_optimal_n(met,
                        criterion = "cutoff",
                        metrics_to_use = "prop_between_dissim",
                        metric_cutoffs = max(vals) + 1,
                        plot = FALSE,
                        verbose = FALSE)
  expect_length(opt$optimal_n$prop_between_dissim, 0)
  expect_false(any(opt$evaluation_df$optimal_n_prop_between_dissim))
  expect_output(print(opt), "prop_between_dissim - none found")
})
