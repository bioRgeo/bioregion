# Inputs -----------------------------------------------------------------------
data("fishmat")
data("vegemat")
data("vegedf")
data("vegesf")

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
fishsim <- similarity(fishmat, metric = c("Jaccard", "Bray"))

cluinfo <- netclu_infomap(vegedf, 
                          seed = 1, 
                          bipartite = TRUE)
cluinfospe <- site_species_subset(cluinfo, node_type = "species")

clulouv <- netclu_louvain(vegesim)

cluhier <- hclu_hierarclust(similarity_to_dissimilarity(fishsim),
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

clualt1 <- clulouv
clualt1$inputs <- clualt1$inputs[-8]
clualt2 <- clulouv
clualt2$inputs <- clualt2$inputs[-9]
clualt3 <- clulouv
clualt3$inputs$data_type <- "unsitneits"
clualt4 <- clulouv
clualt4$inputs$node_type <- "unsitneits"

# Tests for valid outputs ------------------------------------------------------
test_that("valid output", {
  
  # Check output 
  met <- bioregion_metrics(cluinfo, 
                           comat = vegemat,
                           map = vegesf)
  expect_equal(dim(met), c(8,7))
  
  # Check with shuffled comat matrix
  metshuff <- bioregion_metrics(cluinfo, 
                                comat = vegemat_shuff,
                                map = vegesf)
  expect_equal(identical(met, metshuff), TRUE)
  
  # Check without spatial coherence
  metnosf <- bioregion_metrics(cluinfo, 
                               comat = vegemat,
                               map = NULL)
  expect_equal(dim(metnosf), c(8,5))
  expect_equal(metnosf, met[, 1:5])
  
  # Check NbSites (with site_specis_metrics)
  met <- bioregion_metrics(cluhier, 
                           comat = fishmat,
                           map = fishsf)
  expect_equal(length(met), 3)
  
  metest <- site_species_metrics(cluhier,
                                 bioregion_metrics = "CoreTerms", 
                                 comat = fishmat,
                                 verbose = FALSE)
  expect_equal(sum(met$K_3$NbSites == metest$K_3$species_bioregions$n_b[1:3]), 
               3)

})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_warning(
    bioregion_metrics(cluinfo,
                      comat = vegemat,
                      map = NULL,
                      col_bioregion = "udinetaunitseu"),
    "col_bioregion is deprecated.",
    fixed = TRUE)
  
  expect_error(
    bioregion_metrics(1),
    "bioregionalization must be a bioregion.clusters object.",
    fixed = TRUE)
  
  expect_error(
    bioregion_metrics(clualt1),
    "^bioregionalization is a bioregion.cluster object but it has been altered")
  
  expect_error(
    bioregion_metrics(clualt2),
    "^bioregionalization is a bioregion.cluster object but it has been altered")
  
  expect_error(
    bioregion_metrics(clualt3),
    "^bioregionalization is a bioregion.cluster object but it has been altered")
  
  expect_error(
    bioregion_metrics(clualt4),
    "^bioregionalization is a bioregion.cluster object but it has been altered")
  
  expect_error(
    bioregion_metrics(cluinfospe),
    "^No bioregions are assigned to the sites in")
  
  expect_error(
    bioregion_metrics(cluinfo,
                      comat = 1),
    "comat must be a matrix.",
    fixed = TRUE)
  
  expect_error(
    bioregion_metrics(cluinfo,
                      comat = comatneg),
    "Negative value(s) detected in comat!",
    fixed = TRUE)
  
  expect_error(
    bioregion_metrics(cluinfo,
                      comat = comatwnames1),
    "Site ID in bioregionalization and comat do not match!",
    fixed = TRUE)

  expect_error(
    bioregion_metrics(cluinfo,
                      comat = comatwnames2),
    "Site ID in bioregionalization and comat do not match!",
    fixed = TRUE)

  expect_error(
    bioregion_metrics(cluinfo,
                      comat = vegemat,
                      map = "1"),
    "map must be a sf or terra object.",
    fixed = TRUE)

})
