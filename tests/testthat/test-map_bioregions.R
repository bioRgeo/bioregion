# Inputs -----------------------------------------------------------------------
data("fishmat")
data("fishsf")
net <- similarity(fishmat, metric = "Simpson")
clu <- netclu_greedy(net)

cludf <- clu$clusters
cludfNA <- cludf
cludfNA[1,1] <- NA

# Tests for valid outputs ------------------------------------------------------
test_that("valid outputs", {

  map <- map_bioregions(clusters = clu, 
                        geometry = fishsf, 
                        write_clusters = TRUE, 
                        plot = FALSE)
  expect_equal(nrow(map), 338L)
  expect_equal(ncol(map), 3L)
  expect_equal(class(map)[1], "sf")
  expect_equal(class(map)[2], "data.frame")
  
  map <- map_bioregions(clusters = clu, 
                        geometry = fishsf, 
                        write_clusters = TRUE, 
                        plot = TRUE)
  expect_equal(nrow(map), 338L)
  expect_equal(ncol(map), 3L)
  expect_equal(class(map)[1], "sf")
  expect_equal(class(map)[2], "data.frame")
   
})
 
# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
   
  expect_error(
    map_bioregions(clusters = "zz"),
   "If not a bioregion.clusters's object, clusters must be a data.frame.",
   fixed = TRUE)
  
  expect_error(
    map_bioregions(clusters = data.frame(ID = rep(1,5))),
    "clusters must be a data.frame with at least two columns.",
    fixed = TRUE)
  
  expect_error(
    expect_message(
      map_bioregions(clusters = data.frame(ID = rep(1,5), 
                                           ID2 = rep(1,5)),
                     geometry = fishsf),
      "Duplicated site ID detected!",
      fixed = TRUE),
    "The site of clusters should be included in the sites of geometry.",
    fixed = TRUE)
  
  expect_error(
    map_bioregions(clusters = cludfNA, 
                   geometry = fishsf, 
                   write_clusters = TRUE, 
                   plot = FALSE),
    "^NA")
  
  expect_error(
    map_bioregions(clusters = clu, 
                   geometry = "zz"),
    "It seems that the geometry used is not an sf object.",
    fixed = TRUE)

  expect_error(
    map_bioregions(clusters = clu, 
                   geometry = fishsf, 
                   write_clusters = NA, 
                   plot = FALSE),
    "write_clusters must be a boolean.",
    fixed = TRUE)

  expect_error(
    map_bioregions(clusters = clu, 
                   geometry = fishsf, 
                   write_clusters = TRUE, 
                   plot = NA),
    "plot must be a boolean.",
    fixed = TRUE)
  
})
