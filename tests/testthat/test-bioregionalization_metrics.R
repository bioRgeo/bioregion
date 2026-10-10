# Inputs -----------------------------------------------------------------------
data("fishmat")
data("vegemat")
data("vegedf")

install_binaries(verbose = FALSE)

comatneg <- vegemat
comatneg[1,1] <- -1

comatwnames1 <- vegemat
rownames(comatwnames1) <- NULL
colnames(comatwnames1) <- NULL
comatwnames2 <- vegemat
rownames(comatwnames2)[1] <- "enistanutsi"
colnames(comatwnames2)[1] <- "enistanutsi"

vegemat_shuff <- vegemat[sample(dim(vegemat)[1],dim(vegemat)[1]),
                         sample(dim(vegemat)[2],dim(vegemat)[2])]

vegesim <- similarity(vegemat, metric = c("Jaccard", "Simpson", "Sorensen"))
vegedissim <- similarity_to_dissimilarity(vegesim)
fishsim <- similarity(fishmat, metric = c("Jaccard", "Bray"))
fishdissim <- similarity_to_dissimilarity(fishsim)

vegedissimwnames <- vegedissim
vegedissimwnames[vegedissimwnames == "35"] <- "einuastnie"
vegedissimwnames <- vegedissim
vegedissimwnames[vegedissimwnames == "35"] <- "einuastnie"

vegesim_shuff <- vegesim[sample(dim(vegesim)[1],dim(vegesim)[1]),]
vegedissim_shuff <- vegedissim[sample(dim(vegedissim)[1],dim(vegedissim)[1]),]

cluinfo <- netclu_infomap(vegedf, 
                          seed = 1, 
                          bipartite = TRUE)

cluinfospe <- site_species_subset(cluinfo, node_type = "species")

cluhier <- hclu_hierarclust(fishdissim,
                            index = "Jaccard",
                            method = "average",
                            randomize = FALSE,
                            optimal_tree_method = "best",
                            n_clust = c(1,2,3),
                            cut_height = NULL,
                            find_h = TRUE,
                            h_max = 1,
                            h_min = 0,
                            verbose = FALSE)

# Tests for valid outputs ------------------------------------------------------
test_that("valid output", {
  
  eval <-     bioregionalization_metrics(cluhier,
                                         eval_metrics = c("prop_between_dissim", 
                                                          "anosim",
                                                          "mean_endemics",
                                                          "tot_endemics"), 
                                         dissimilarity = fishdissim,
                                         comat = fishmat,
                                         anosim_permutations = 999)
  
  eval

})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_warning(
    bioregionalization_metrics(cluinfo,
                               eval_metrics = c("prop_between_dissim", "anosim"), 
                               dissimilarity = vegedissim,
                               comat = NULL,
                               eval_metric = "isetnusi",
                               net = NULL,
                               site_col = NULL, 
                               species_col = NULL),
    "eval_metric is deprecated.",
    fixed = TRUE)
  
  expect_warning(
    bioregionalization_metrics(cluinfo,
                               eval_metrics = c("prop_between_dissim", "anosim"), 
                               dissimilarity = vegedissim,
                               comat = NULL,
                               net = "insetnu",
                               site_col = NULL, 
                               species_col = NULL),
    "net is deprecated.",
    fixed = TRUE)
  
  expect_warning(
    bioregionalization_metrics(cluinfo,
                               eval_metrics = c("prop_between_dissim", "anosim"), 
                               dissimilarity = vegedissim,
                               comat = NULL,
                               net = NULL,
                               site_col = "ietis", 
                               species_col = NULL),
    "site_col is deprecated.",
    fixed = TRUE)
  
  expect_warning(
    bioregionalization_metrics(cluinfo,
                               eval_metrics = c("prop_between_dissim", "anosim"), 
                               dissimilarity = vegedissim,
                               comat = NULL,
                               net = NULL,
                               site_col = NULL, 
                               species_col = "ietis"),
    "species_col is deprecated.",
    fixed = TRUE)
  
  expect_error(
    bioregionalization_metrics(1),
    "bioregionalization must be a bioregion.clusters object.",
    fixed = TRUE)
  
  expect_error(
    bioregionalization_metrics(cluinfospe),
    "No bioregion are assigned to the site in bioregionalization.",
    fixed = TRUE)
  
  expect_error(
    bioregionalization_metrics(cluinfo,
                               eval_metrics = 1, 
                               dissimilarity,
                               dissimilarity_index = names(dissimilarity)[3], 
                               comat = NULL),
    "eval_metrics must be a character.",
    fixed = TRUE)
  
  expect_error(
    bioregionalization_metrics(cluinfo,
                               eval_metrics = c("prop_between_dissim", "instenusti"), 
                               dissimilarity,
                               dissimilarity_index = names(dissimilarity)[3], 
                               comat = NULL),
    "^One or several evaluation metrics chosen are not")
  
  expect_error(
    bioregionalization_metrics(cluinfo,
                               eval_metrics = "prop_between_dissim", 
                               dissimilarity = NULL,
                               comat = NULL),
    "At least dissimilarity or comat should be provided.",
    fixed = TRUE)
  
  expect_warning(
    bioregionalization_metrics(cluinfo,
                               eval_metrics = c("prop_between_dissim", "mean_endemics"),
                               dissimilarity = NULL,
                               comat = vegemat),
    "^Some metrics")
  
  expect_warning(
    bioregionalization_metrics(cluinfo,
                               eval_metrics = c("prop_between_dissim", "mean_endemics"), 
                               dissimilarity = vegedissim,
                               comat = NULL),
    "^Some metrics")
  
  expect_error(
    expect_warning(
      bioregionalization_metrics(cluinfo,
                                eval_metrics = "prop_between_dissim", 
                                dissimilarity = NULL,
                                comat = vegemat),
      "^Some metrics"),
    "At least one metric with the appropriate input should be specified.",
  fixed = TRUE)
  
  expect_error(
    bioregionalization_metrics(cluinfo,
                               eval_metrics = "prop_between_dissim", 
                               dissimilarity = 1,
                               comat = NULL),
    "^dissimilarity should be a bioregion.pairwise object created by")
  
  expect_error(
    bioregionalization_metrics(cluinfo,
                               eval_metrics = "prop_between_dissim", 
                               dissimilarity = vegedissim,
                               dissimilarity_index = c(1, 1),
                               comat = NULL),
    "dissimilarity_index must be of length 1.", 
    fixed = TRUE)
  
  expect_error(
    bioregionalization_metrics(cluinfo,
                               eval_metrics = "prop_between_dissim", 
                               dissimilarity = vegedissimwnames,
                               comat = NULL),
    "^Some sites are not found in dissimilarity:")
  
  expect_error(
    bioregionalization_metrics(cluinfo,
                               eval_metrics = "mean_endemics", 
                               dissimilarity = NULL,
                               comat = 1),
    "comat must be a matrix.", 
    fixed = TRUE)
  
  expect_error(
    bioregionalization_metrics(cluinfo,
                               eval_metrics = "mean_endemics", 
                               dissimilarity = NULL,
                               comat = comatneg),
    "Negative value(s) detected in comat!", 
    fixed = TRUE)
  
  expect_error(
    bioregionalization_metrics(cluinfo,
                               eval_metrics = "mean_endemics", 
                               dissimilarity = NULL,
                               comat = comatwnames1),
    "^Some sites are not found in comat:")
  
  expect_error(
    bioregionalization_metrics(cluinfo,
                               eval_metrics = "mean_endemics", 
                               dissimilarity = NULL,
                               comat = comatwnames2),
    "^Some sites are not found in comat:")
  
})
