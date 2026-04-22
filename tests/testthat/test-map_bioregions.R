# Inputs -----------------------------------------------------------------------
data("fishmat")
data("vegemat")
data("vegedf")

install_binaries(verbose = FALSE)

vegesim <- similarity(vegemat, metric = c("Simpson"))
fishsim <- similarity(fishmat, metric = c("Jaccard"))

cluinfo <- netclu_infomap(vegedf, 
                          seed = 1, 
                          bipartite = TRUE)
cluinfocol <- bioregion_colors(cluinfo)

cluinfospe <- site_species_subset(cluinfo, 
                                  node_type = "species")

cluhier <- hclu_hierarclust(similarity_to_dissimilarity(fishsim),
                            index = "Jaccard",
                            seed = 1,
                            method = "average",
                            randomize = FALSE,
                            optimal_tree_method = "best",
                            n_clust = c(1,2,3),
                            cut_height = NULL,
                            find_h = TRUE,
                            h_max = 1,
                            h_min = 0,
                            verbose = FALSE)
#cluhiercol <- bioregion_colors(cluhier,
#                               palette = "Vivid",
#                               cluster_ordering = "n_sites",
#                               cutoff_insignificant = NULL)

cluhier7 <- hclu_hierarclust(similarity_to_dissimilarity(fishsim),
                             index = "Jaccard",
                             seed = 1,
                             method = "average",
                             randomize = FALSE,
                             optimal_tree_method = "best",
                             n_clust = 3:7,
                             cut_height = NULL,
                             find_h = TRUE,
                             h_max = 1,
                             h_min = 0,
                             verbose = FALSE)
cluhiercol7 <- bioregion_colors(cluhier7,
                                palette = "Vivid",
                                cluster_ordering = "n_sites",
                                cutoff_insignificant = NULL)

# sf
data("fishsf")
data("vegesf")

vegenodf <- vegesf
class(vegenodf) <- "sf"
vegesf1 <- vegesf[,-1]
vegesfw <- vegesf[-(1:10),]
vegesfplus <- rbind(vegesf, vegesf)
vegesfplus$Site[729:dim(vegesfplus)[1]] <- 
  as.character(as.numeric(vegesfplus$Site[729:dim(vegesfplus)[1]]) + 10000)
vegesfplus <- vegesfplus[sample(dim(vegesfplus)[1]),]

# Tests for valid outputs ------------------------------------------------------
test_that("valid outputs", {
  
  # Check that returned map correspond to the bioregionalization
  map <- map_bioregions(cluinfo, 
                        map = vegesf, 
                        partition_index = NULL,
                        map_as_output = TRUE, 
                        plot = FALSE)
  expect_equal(sum(map$Site == cluinfo$clusters$ID[1:715]), 715)
  expect_equal(sum(map$K_8 == cluinfo$clusters$K_8[1:715]), 715)  
  
  map <- map_bioregions(cluinfo, 
                        map = vegesfplus, 
                        partition_index = NULL,
                        map_as_output = TRUE, 
                        plot = FALSE)
  expect_equal(sum(map$Site == cluinfo$clusters$ID[1:715]), 715)
  expect_equal(sum(map$K_8 == cluinfo$clusters$K_8[1:715]), 715)  
  
  # Check output with partition_index (character vs numeric)
  map1 <- map_bioregions(cluhier, 
                        map = fishsf, 
                        partition_index = c("K_2","K_3"),
                        map_as_output = TRUE, 
                        plot = FALSE)
  map2 <- map_bioregions(cluhier, 
                         map = fishsf, 
                         partition_index = c(3,4),
                         map_as_output = TRUE, 
                         plot = FALSE)
  #map3 <- map_bioregions(cluhiercol, 
  #                       map = fishsf, 
  #                       partition_index = c(3,4),
  #                       map_as_output = TRUE, 
  #                       plot = FALSE)
  expect_equal(identical(map1, map2), TRUE)
  #expect_equal(identical(map2, map3), TRUE)
  
  # Check plot (one cluster & no color)
  expect_no_error_plotless(
    map_bioregions(cluinfo, 
                   map = vegesf, 
                   partition_index = NULL,
                   map_as_output = FALSE, 
                   plot = TRUE)
  )
  
  # Check plot (several clusters & no color)
  expect_no_error_plotless(
    map_bioregions(cluinfocol, 
                   map = vegesf, 
                   partition_index = NULL,
                   map_as_output = FALSE, 
                   plot = TRUE)
  )
  
  # Check plot (one cluster & color)
  expect_no_error_plotless(
    map_bioregions(cluhier, 
                   map = fishsf, 
                   partition_index = NULL,
                   map_as_output = FALSE, 
                   plot = TRUE)
  )
  
  # Check plot (several clusters & color)
  #expect_no_error_plotless(
  #  map_bioregions(cluhiercol, 
  #                 map = fishsf, 
  #                 partition_index = NULL,
  #                 map_as_output = FALSE, 
  #                 plot = TRUE)
  #)
  
  # Check plot (7 clusters & no color)
  expect_no_error_plotless(
    map_bioregions(cluhier7,
                   map = fishsf,
                   partition_index = NULL,
                   map_as_output = FALSE,
                   plot = TRUE)
  )
  
  # Check plot (7 clusters & color)
  expect_no_error(
    map_bioregions(cluhiercol7,
                   map = fishsf,
                   partition_index = NULL,
                   map_as_output = FALSE,
                   plot = TRUE)
  )

})
 
# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    map_bioregions(cluinfo,
                   map_as_output = c(TRUE,1)),
    "map_as_output must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    map_bioregions(cluinfo,
                   map_as_output = 1),
    "map_as_output must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    map_bioregions(cluinfo,
                   plot = c(TRUE,1)),
    "plot must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    map_bioregions(cluinfo,
                   plot = 1),
    "plot must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    map_bioregions(cluinfo,
                   map_as_output = FALSE,
                   plot = FALSE),
    "^At least one argument among")
  
  expect_error(
    map_bioregions(1),
    "bioregionalization must be a bioregion.clusters object.",
    fixed = TRUE)
  
  expect_error(
    map_bioregions(cluinfospe),
    "^No bioregions are assigned to the sites in")
  
  expect_error(
    map_bioregions(cluinfo,
                   partition_index = c(1,1)),
    "^Duplicated values detected in")
  
  expect_error(
    map_bioregions(cluinfo,
                   partition_index = TRUE),
    "partition_index should be numeric or character.",
    fixed = TRUE)
  
  expect_error(
    map_bioregions(cluinfo,
                   partition_index = "ID"),
    "^If partition_index is a character")


  expect_error(
    map_bioregions(cluinfo,
                   partition_index = "K_2"),
    "^If partition_index is a character")    
  
  expect_error(
    map_bioregions(cluinfo,
                   partition_index = c(1.1,2)),
    "^If partition_index is numeric") 
  
  expect_error(
    map_bioregions(cluinfo,
                   partition_index = c(1,2)),
    "^partition_index should be strictly higher") 
  
  expect_error(
    map_bioregions(cluinfo,
                   partition_index = 3),
    "^partition_index should be lower or equal to") 
  
  expect_error(
    map_bioregions(cluinfo,
                   map = 1),
    "map must be a sf or terra object.",
    fixed = TRUE)
  
  expect_error(
    map_bioregions(cluinfo,
                   map = vegenodf),
    "map must be a sf data.frame.",
    fixed = TRUE)
  
  expect_error(
    map_bioregions(cluinfo,
                   map = vegesf1),
    "^map must be a sf data.frame with")
  
  expect_error(
    map_bioregions(cluinfo,
                   map = vegesfw),
    "^Some sites are not found in map")
  
})

# Tests for external packages --------------------------------------------------
test_that("terra", {
  
  skip_if_not_installed_quiet("terra")
  quietly(library(terra))
  
  # SpatVector
  vegesv <- terra::vect(vegesf)
  vegesv0 <- terra::vect(vegesf1)
  vegesvw <- terra::vect(vegesfw)
  
  expect_error(
    map_bioregions(cluinfo,
                   map = vegesv0),
    "^map must be a SpatVector with")
  
  expect_error(
    map_bioregions(cluinfo,
                   map = vegesvw),
    "^Some sites are not found in map")
  
  # Check that returned map correspond to the bioregionalization
  map <- map_bioregions(cluinfo, 
                        map = vegesv, 
                        partition_index = NULL,
                        map_as_output = TRUE, 
                        plot = FALSE)
  expect_equal(sum(map$Site == cluinfo$clusters$ID[1:715]), 715)
  expect_equal(sum(map$K_8 == cluinfo$clusters$K_8[1:715]), 715)  
  
  # SpatRaster
  r <- terra::rast(vegesv, resolution = 10000)
  vegesr <- terra::rasterize(vegesv, r, field = "Site")
  r <- terra::rast(vegesv0, resolution = 10000)
  vegesr0 <- terra::rasterize(vegesv0, r)
  r <- terra::rast(vegesvw, resolution = 10000)
  vegesrw <- terra::rasterize(vegesvw, r, field = "Site")
  
  expect_error(
    map_bioregions(cluinfo,
                   map = vegesr0),
    "^Some sites are not found in map")
  
  expect_error(
    map_bioregions(cluinfo,
                   map = vegesrw),
    "^Some sites are not found in map")
  
  # Check that returned map correspond to the bioregionalization
  map <- map_bioregions(cluinfo, 
                        map = vegesr, 
                        partition_index = NULL,
                        map_as_output = TRUE, 
                        plot = FALSE)
  expect_equal(sum(map$Site == cluinfo$clusters$ID[1:715]), 715)
  expect_equal(sum(map$K_8 == cluinfo$clusters$K_8[1:715]), 715) 
  
  
})
