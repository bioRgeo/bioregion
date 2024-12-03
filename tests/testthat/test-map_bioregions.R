# # Preamble code ----------------------------------------------------------------
# data("fishmat")
# data("fishsf")
# net <- similarity(fishmat, metric = "Simpson")
# clu <- netclu_greedy(net)
# 
# # Tests for valid outputs -----------------------------------------------------
# map <- map_bioregions(clu, fishsf, write_clusters = TRUE, plot = FALSE)
# 
# test_that("number of columns in output", {
#   
#   expect_equal(nrow(map), 338L)
#   expect_equal(ncol(map), 3L)
#   expect_equal(class(map)[1], "sf")
#   expect_equal(class(map)[2], "data.frame")
#   
# })
# 
# # Tests for invalid inputs ----------------------------------------------------
# test_that("error messages with wrong inputs", {
#   expect_error(
#     map_bioregions("zz"),
#     "If not a bioregion.clusters's object, clusters must be a data.frame.",
#     fixed = TRUE)
#   
#   expect_error(
#     map_bioregions(clu, "zz"),
#     "It seems that the geometry used is not an sf object.",
#     fixed = TRUE)
#   
#   expect_error(
#     map_bioregions(clu, fishsf, write_clusters = NA, plot = FALSE),
#     "write_clusters must be a boolean.",
#     fixed = TRUE)
#   
#   expect_error(
#     map_bioregions(clu, fishsf, write_clusters = TRUE, plot = NA),
#     "plot must be a boolean.",
#     fixed = TRUE)
# })
