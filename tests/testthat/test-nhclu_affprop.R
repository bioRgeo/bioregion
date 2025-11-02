# Inputs -----------------------------------------------------------------------
sim <- similarity(fishmat, 
                  metric = "all")
dissim <- similarity_to_dissimilarity(sim)
                  
df <- data.frame(ID1 = sim$Site1, ID2 = sim$Site2, W = sim$Simpson)

d <- dist(fishmat)
mat <- fishmat
rownames(mat) <- NULL
colnames(mat) <- NULL
d2 <- dist(mat)
d3 <- d
d3[1] <- "1"
d4 <- d
d4[1] <- NA

uni <- data.frame(
  Site1 = c("c", "b", "a"),
  Site2 = c("a", "c", "b"),
  Weight = c(10, 100, 1))

unina1 <- uni
unina1[1,1] <- NA

unina2 <- uni
unina2[1,3] <- NA

unichar <- uni
unichar$Weight <- unichar$Site1

uni2 <- data.frame(
  Site1 = c("a", "c", "a"),
  Site2 = c("a", "a", "b"),
  Weight = c(10, 100, 1))

uni3 <- data.frame(
  Site1 = c("c", "b", "a"),
  Site2 = c("a", "a", "b"),
  Weight = c(10, 100, 1))

uni4 <- data.frame(
  Site1 = c("c", "a", "a"),
  Site2 = c("a", "b", "b"),
  Weight = c(10, 100, 1))

# Tests for valid outputs ------------------------------------------------------
test_that("valid output", {
  
  clust <- nhclu_affprop(sim,
                         index = "Simpson",
                         seed = NULL,
                         p = NA,
                         q = NA,
                         maxits = 1000,
                         convits = 100,
                         lam = 0.9,
                         details = FALSE,
                         nonoise = FALSE,
                         K = NULL,
                         prc = NULL,
                         bimaxit = NULL,
                         exact = NULL,
                         algorithm_in_output = TRUE,
                         verbose = TRUE)
  expect_equal(inherits(clust, "bioregion.clusters"), TRUE)
  expect_equal(clust$name, "nhclu_affprop")
  expect_equal(clust$args$index, "Simpson")
  expect_equal(clust$args$seed, NULL)
  expect_equal(clust$args$p, NA)
  expect_equal(clust$args$q, NA)
  expect_equal(clust$args$maxits, 1000)
  expect_equal(clust$args$convits, 100)
  expect_equal(clust$args$lam, 0.9)
  expect_equal(clust$args$details, FALSE)
  expect_equal(clust$args$nonoise, FALSE)
  expect_equal(clust$args$algorithm_in_output, TRUE)
  expect_equal(clust$args$verbose, TRUE)
  expect_equal(clust$inputs$bipartite, FALSE)
  expect_equal(clust$inputs$weight, TRUE)
  expect_equal(clust$inputs$pairwise, TRUE)
  expect_equal(clust$inputs$pairwise_metric, "Simpson")
  expect_equal(clust$inputs$dissimilarity, FALSE)
  expect_equal(clust$inputs$nb_sites, 338)
  expect_equal(clust$inputs$hierarchical, FALSE)
  expect_equal(clust$inputs$data_type, "occurrence")
  expect_equal(clust$inputs$node_type, "site")
  expect_equal(sum(attr(clust$clusters, "node_type")=="site"), 
               dim(clust$clusters)[1])
  
  clust <- nhclu_affprop(sim,
                         index = 7,
                         seed = NULL,
                         p = NA,
                         q = NA,
                         maxits = 1000,
                         convits = 100,
                         lam = 0.9,
                         details = FALSE,
                         nonoise = FALSE,
                         K = NULL,
                         prc = NULL,
                         bimaxit = NULL,
                         exact = NULL,
                         algorithm_in_output = TRUE,
                         verbose = TRUE)
  expect_equal(clust$args$index, 7)
  expect_equal(clust$inputs$pairwise_metric, "Bray")
  
  clust <- nhclu_affprop(sim,
                         index = "Simpson",
                         seed = NULL,
                         p = NA,
                         q = NA,
                         maxits = 1000,
                         convits = 100,
                         lam = 0.9,
                         details = FALSE,
                         nonoise = TRUE,
                         K = NULL,
                         prc = NULL,
                         bimaxit = NULL,
                         exact = NULL,
                         algorithm_in_output = TRUE,
                         verbose = TRUE)
  expect_equal(clust$cluster_info[1,1], "K_16")
  expect_equal(clust$cluster_info[1,2], 16)
  
  clust <- nhclu_affprop(sim,
                         index = "Simpson",
                         seed = NULL,
                         p = NA,
                         q = NA,
                         maxits = 1000,
                         convits = 100,
                         lam = 0.9,
                         details = FALSE,
                         nonoise = TRUE,
                         K = NULL,
                         prc = NULL,
                         bimaxit = NULL,
                         exact = NULL,
                         algorithm_in_output = TRUE,
                         verbose = FALSE)
  expect_equal(clust$cluster_info[1,1], "K_16")
  expect_equal(clust$cluster_info[1,2], 16)
  
  clust <- nhclu_affprop(sim,
                         index = "Simpson",
                         seed = 2,
                         p = NA,
                         q = NA,
                         maxits = 1000,
                         convits = 100,
                         lam = 0.9,
                         details = FALSE,
                         nonoise = FALSE,
                         K = NULL,
                         prc = NULL,
                         bimaxit = NULL,
                         exact = NULL,
                         algorithm_in_output = TRUE,
                         verbose = TRUE)
  expect_equal(clust$cluster_info[1,1], "K_16")
  expect_equal(clust$cluster_info[1,2], 16)
  
  clust <- nhclu_affprop(sim,
                         index = "Simpson",
                         seed = 2,
                         p = NA,
                         q = NA,
                         maxits = 1000,
                         convits = 100,
                         lam = 0.9,
                         details = FALSE,
                         nonoise = FALSE,
                         K = NULL,
                         prc = NULL,
                         bimaxit = NULL,
                         exact = NULL,
                         algorithm_in_output = TRUE,
                         verbose = FALSE)
  expect_equal(clust$cluster_info[1,1], "K_16")
  expect_equal(clust$cluster_info[1,2], 16)
  
  clust <- nhclu_affprop(sim,
                         index = "Simpson",
                         seed = 2,
                         p = NA,
                         q = NA,
                         maxits = 1000,
                         convits = 100,
                         lam = 0.9,
                         details = FALSE,
                         nonoise = FALSE,
                         K = NULL,
                         prc = NULL,
                         bimaxit = NULL,
                         exact = NULL,
                         algorithm_in_output = TRUE)
  expect_equal(clust$cluster_info[1,1], "K_16")
  expect_equal(clust$cluster_info[1,2], 16)
  
  clust <- nhclu_affprop(d,
                         index = "Simpson",
                         seed = 2,
                         p = NA,
                         q = NA,
                         maxits = 1000,
                         convits = 100,
                         lam = 0.9,
                         details = FALSE,
                         nonoise = FALSE,
                         K = NULL,
                         prc = NULL,
                         bimaxit = NULL,
                         exact = NULL,
                         algorithm_in_output = TRUE,
                         verbose = FALSE)
  expect_equal(clust$cluster_info[1,1], "K_2")
  expect_equal(clust$cluster_info[1,2], 2)
  
  quietly(
    clust1 <- nhclu_affprop(sim,
                            index = "Simpson",
                            seed = 2,
                            p = NA,
                            q = NA,
                            maxits = 1000,
                            convits = 100,
                            lam = 0.9,
                            details = FALSE,
                            nonoise = FALSE,
                            K = 6,
                            prc = 1,
                            bimaxit = 5,
                            exact = TRUE,
                            algorithm_in_output = FALSE,
                            verbose = TRUE)
  )
  expect_equal(clust1$cluster_info[1,1], "K_6")
  expect_equal(clust1$cluster_info[1,2], 6)
  
  clust1 <- nhclu_affprop(sim,
                          index = "Simpson",
                          seed = 2,
                          p = NA,
                          q = NA,
                          maxits = 1000,
                          convits = 100,
                          lam = 0.9,
                          details = FALSE,
                          nonoise = FALSE,
                          K = 6,
                          prc = 1,
                          bimaxit = 5,
                          exact = TRUE,
                          algorithm_in_output = FALSE,
                          verbose = FALSE)
  expect_equal(clust1$cluster_info[1,1], "K_6")
  expect_equal(clust1$cluster_info[1,2], 6)
  
  clust2 <- nhclu_affprop(sim,
                          index = "Simpson",
                          seed = 2,
                          p = NA,
                          q = NA,
                          maxits = 1000,
                          convits = 100,
                          lam = 0.9,
                          details = FALSE,
                          nonoise = FALSE,
                          K = 6,
                          prc = 1,
                          bimaxit = 5,
                          exact = TRUE,
                          algorithm_in_output = FALSE,
                          verbose = FALSE)
  expect_equal(clust2$cluster_info[1,1], "K_6")
  expect_equal(clust2$cluster_info[1,2], 6)
  
  expect_equal(sum(clust1$clusters$K_6==clust2$clusters$K_6), 338)
  
  quietly(
    clust <- nhclu_affprop(sim,
                           index = "Simpson",
                           seed = NULL,
                           p = NA,
                           q = NA,
                           maxits = 1000,
                           convits = 100,
                           lam = 0.9,
                           details = FALSE,
                           nonoise = TRUE,
                           K = 6,
                           prc = 1,
                           bimaxit = 5,
                           exact = TRUE,
                           algorithm_in_output = FALSE,
                           verbose = TRUE)
  )
  expect_equal(clust$cluster_info[1,1], "K_6")
  expect_equal(clust$cluster_info[1,2], 6)
  
  clust <- nhclu_affprop(sim,
                         index = "Simpson",
                         seed = NULL,
                         p = NA,
                         q = NA,
                         maxits = 1000,
                         convits = 100,
                         lam = 0.9,
                         details = FALSE,
                         nonoise = TRUE,
                         K = 6,
                         prc = 1,
                         bimaxit = 5,
                         exact = TRUE,
                         algorithm_in_output = FALSE,
                         verbose = FALSE)
  expect_equal(clust$cluster_info[1,1], "K_6")
  expect_equal(clust$cluster_info[1,2], 6)
  
  # Test data_type with different similarity metrics
  # Create similarity objects with different metrics
  sim_occ <- similarity(fishmat, metric = c("Simpson", "Jaccard"))
  sim_abund <- similarity(fishmat, metric = "Bray")
  
  clust <- nhclu_affprop(sim_occ, 
                         index = "Simpson",
                         verbose = FALSE)
  expect_equal(clust$inputs$data_type, "occurrence")
  
  clust <- nhclu_affprop(sim_occ, 
                         index = "Jaccard",
                         verbose = FALSE)
  expect_equal(clust$inputs$data_type, "occurrence")
  
  clust <- nhclu_affprop(sim_abund, 
                         index = "Bray",
                         verbose = FALSE)
  expect_equal(clust$inputs$data_type, "abundance")
  
})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    nhclu_affprop("zz"),
    "^similarity is not a bioregion.pairwise object")
  
  expect_error(
    nhclu_affprop(dissim),
    "^similarity seems to be a dissimilarity object")
  
  expect_error(
    nhclu_affprop(sim[,1:2]),
    "similarity must be a data.frame with at least three columns.",
    fixed = TRUE)
  
  expect_error(
    nhclu_affprop(uni[,-3]),
    "similarity must be a data.frame with at least three columns.",
    fixed = TRUE)
  
  expect_error(
    nhclu_affprop(unina1),
    "NA(s) detected in similarity.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_affprop(unina2),
    "NA(s) detected in similarity.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_affprop(uni4),
    "The first two columns of similarity contain duplicated pairs of sites!",
    fixed = TRUE) 
  
  expect_error(
    nhclu_affprop(uni3),
    "^The first two columns of similarity contain") 
  
  expect_error(
    nhclu_affprop(uni2),
    "similarity contains rows with the same site on both columns!",
    fixed = TRUE) 
  
  expect_error(
    nhclu_affprop(unichar),
    "The weight column must be numeric.",
    fixed = TRUE)  
  
  expect_message(
    nhclu_affprop(d2),
    "No labels detected, they have been assigned automatically.",
    fixed = TRUE) 
  
  expect_error(
    nhclu_affprop(d3),
    "similarity must be numeric.",
    fixed = TRUE) 
  
  expect_error(
    nhclu_affprop(d4),
    "NA(s) detected in similarity.",
    fixed = TRUE) 
  
  expect_error(
    nhclu_affprop(sim, 
                  seed =  c("zz","zz")),
    "seed must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_affprop(sim, 
                  seed = "zz"),
    "seed must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_affprop(sim, 
                  seed = 1.1),
    "seed must be an integer.",
    fixed = TRUE) 
  
  expect_error(
    nhclu_affprop(sim, 
                  seed = -1),
    "seed must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    nhclu_affprop(sim, 
                  seed = 0),
    "seed must be strictly higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    nhclu_affprop(sim, 
                  p = c("1","2")),
    "p must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_affprop(sim, 
                  p = "1"),
    "p must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_affprop(sim, 
                  p = c(1,2)),
    "^If p is a vector")  
  
  expect_error(
    nhclu_affprop(sim, 
                   q =  c("zz","zz")),
    "q must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_affprop(sim, 
                   q = "zz"),
    "q must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_affprop(sim, 
                  q = -1),
    "q must be higher than 0.",
    fixed = TRUE)
  
  expect_error(
    nhclu_affprop(sim, 
                  q = 1.1),
    "q should be in the interval [0, 1].",
    fixed = TRUE)
  
  expect_error(
    nhclu_affprop(sim,  
                  maxits =  c("zz","zz")),
    "maxits must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_affprop(sim,  
                  maxits = "zz"),
    "maxits must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_affprop(sim,  
                  maxits = 1.1),
    "maxits must be an integer.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_affprop(sim,  
                  maxits = -1),
    "maxits must be higher than 0.",
    fixed = TRUE) 
  
  expect_error(
    nhclu_affprop(sim,  
                  convits =  c("zz","zz")),
    "convits must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_affprop(sim,  
                  convits = "zz"),
    "convits must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_affprop(sim,  
                  convits = 1.1),
    "convits must be an integer.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_affprop(sim,  
                  convits = -1),
    "convits must be higher than 0.",
    fixed = TRUE)
  
  expect_error(
    nhclu_affprop(sim, 
                  lam =  c("zz","zz")),
    "lam must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_affprop(sim, 
                  lam = "zz"),
    "lam must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_affprop(sim, 
                  lam = -1),
    "lam must be higher than 0.",
    fixed = TRUE)
  
  expect_error(
    nhclu_affprop(sim, 
                  lam = 0.4),
    "lam should be in the interval [0.5, 1).",
    fixed = TRUE)
  
  expect_error(
    nhclu_affprop(sim, 
                  lam = 1),
    "lam should be in the interval [0.5, 1).",
    fixed = TRUE)
  
  expect_error(
    nhclu_affprop(sim, 
                  details = 1),
    "details must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    nhclu_affprop(sim,
                  details = c("zz","zz")),
    "details must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    nhclu_affprop(sim, 
                  nonoise = 1),
    "nonoise must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    nhclu_affprop(sim,
                  nonoise = c("zz","zz")),
    "nonoise must be of length 1.",
    fixed = TRUE)
  
  expect_message(
    nhclu_affprop(sim,
                  seed = 1,
                  nonoise = TRUE),
    "^A random number generator")
  
  expect_error(
    nhclu_affprop(sim,  
                  K =  c("zz","zz")),
    "K must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_affprop(sim,  
                  K = "zz"),
    "K must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_affprop(sim,  
                  K = 1.1),
    "K must be an integer.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_affprop(sim,  
                  K = -1),
    "K must be strictly higher than 0.",
    fixed = TRUE)
  
  expect_error(
    nhclu_affprop(sim,  
                  K = 0),
    "K must be strictly higher than 0.",
    fixed = TRUE)
  
  expect_error(
    nhclu_affprop(sim,  
                  K = 1,
                  prc = NULL,
                  bimaxit = NULL,
                  exact = NULL),
    "^When K is not NULL, you need to define prc.")
  
  expect_error(
    nhclu_affprop(sim,  
                  K = 1,
                  prc = 1,
                  bimaxit = NULL,
                  exact = NULL),
    "When K is not NULL, you need to define bimaxit.",
    fixed = TRUE)
  
  expect_error(
    nhclu_affprop(sim,  
                  K = 1,
                  prc = 1,
                  bimaxit = 1,
                  exact = NULL),
    "When K is not NULL, you need to define exact.",
    fixed = TRUE)
  
  expect_message(
    nhclu_affprop(sim,  
                  K = NULL,
                  prc = 1),
    "^prc argument will be")
  
  expect_error(
    nhclu_affprop(sim, 
                  K = 1,
                  prc =  c("zz","zz"),
                  bimaxit = 1,
                  exact = TRUE),
    "prc must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_affprop(sim, 
                  K = 1,
                  prc = "zz",
                  bimaxit = 1,
                  exact = TRUE),
    "prc must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_affprop(sim, 
                  K = 1,
                  prc = -1,
                  bimaxit = 1,
                  exact = TRUE),
    "prc must be higher than 0.",
    fixed = TRUE)
  
  expect_error(
    nhclu_affprop(sim, 
                  K = 1,
                  prc = 100.1,
                  bimaxit = 1,
                  exact = TRUE),
    "prc should be in the interval [0, 100].",
    fixed = TRUE)
  

  quietly(
    expect_message(
      nhclu_affprop(sim,  
                    K = NULL,
                    prc = 1,
                    bimaxit = 1),
      "^bimaxit argument will be")
  )
  
  expect_error(
    nhclu_affprop(sim,  
                  K = 1,
                  prc = 1,
                  bimaxit =  c("zz","zz"),
                  exact = TRUE),
    "bimaxit must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_affprop(sim,  
                  K = 1,
                  prc = 1,
                  bimaxit = "zz",
                  exact = TRUE),
    "bimaxit must be numeric.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_affprop(sim,  
                  K = 1,
                  prc = 1,
                  bimaxit = 1.1,
                  exact = TRUE),
    "bimaxit must be an integer.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_affprop(sim,  
                  K = 1,
                  prc = 1,
                  bimaxit = -1,
                  exact = TRUE),
    "bimaxit must be strictly higher than 0.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_affprop(sim,  
                  K = 1,
                  prc = 1,
                  bimaxit = 0,
                  exact = TRUE),
    "bimaxit must be strictly higher than 0.",
    fixed = TRUE)  
  
  quietly(
    expect_message(
      nhclu_affprop(sim,  
                    K = NULL,
                    prc = 1,
                    bimaxit = 1,
                    exact = TRUE),
      "^exact argument will be")
  )
  
  expect_error(
    nhclu_affprop(sim,  
                  K = 1,
                  prc = 1,
                  bimaxit =  1,
                  exact = c("zz","zz")),
    "exact must be of length 1.",
    fixed = TRUE)  
  
  expect_error(
    nhclu_affprop(sim,  
                  K = 1,
                  prc = 1,
                  bimaxit = 1,
                  exact = "zz"),
    "exact must be a boolean.",
    fixed = TRUE)   

  expect_error(
    nhclu_affprop(sim, 
                  algorithm_in_output = 1),
    "algorithm_in_output must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    nhclu_affprop(sim,
                  algorithm_in_output = c("zz","zz")),
    "algorithm_in_output must be of length 1.",
    fixed = TRUE)
  
})
