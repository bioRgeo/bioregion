# Tests bioregion.site.species.metrics
test_that("Test generic functions for bioregion.site.species.metrics", {
  
  # Create a site_species_metrics object for testing
  data("fishmat")
  data("vegemat")
  data("vegedf")
  
  vegesim <- similarity(vegemat, metric = c("Jaccard", "Simpson", "Sorensen"))
  
  cluinfo <- netclu_infomap(vegedf, 
                            seed = 1, 
                            bipartite = TRUE)
  
  ind <- site_species_metrics(bioregionalization = cluinfo,
                              bioregion_indices = "all",
                              bioregionalization_indices = "all",
                              data_type = "both",
                              node_type = "both",
                              comat = vegemat,
                              similarity = vegesim,
                              index = 3,
                              verbose = FALSE)
  
  # Test class
  expect_s3_class(ind, "bioregion.site.species.metrics")
  
  # Test print method
  expect_output(print(ind), "Site and species contribution metrics")
  expect_output(print(ind), "Input summary:")
  expect_output(print(ind), "Node type:")
  expect_output(print(ind), "Computed indices:")
  expect_output(print(ind), "Data preview:")
  expect_output(print(ind), "Access data with:")
  
  # Test print with n_preview = 0 (no data preview, only dimensions)
  expect_output(print(ind, n_preview = 0), "Available components:")
  
  # Test print with custom n_preview
  expect_output(print(ind, n_preview = 5), "Site and species contribution metrics")
  
  # Test str method
  expect_output(str(ind), "bioregion.site.species.metrics object")
  expect_output(str(ind), "Partitions:")
  expect_output(str(ind), "Node type:")
  
  # Test summary method
  expect_output(summary(ind), "Summary of site and species contribution metrics")
  expect_output(summary(ind), "Settings:")
  expect_output(summary(ind), "Number of partitions:")
  expect_output(summary(ind), "Node type:")
  
  # Test summary with show_top_contributors = FALSE
  expect_output(summary(ind, show_top_contributors = FALSE), 
                "Summary of site and species contribution metrics")
  
  # Test summary with custom n_top
  expect_output(summary(ind, n_top = 3), 
                "Summary of site and species contribution metrics")
  
  # Test attributes
  expect_true(!is.null(attr(ind, "n_partitions")))
  expect_true(!is.null(attr(ind, "node_type")))
  expect_equal(attr(ind, "node_type"), "both")
  
  # Test with hierarchical clustering (multiple partitions)
  fishsim <- similarity(fishmat, metric = c("Jaccard", "Bray"))
  
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
  
  ind_hier <- site_species_metrics(bioregionalization = cluhier,
                                   bioregion_indices = "all",
                                   bioregionalization_indices = "all",
                                   data_type = "auto",
                                   node_type = "site",
                                   comat = fishmat,
                                   similarity = fishsim,
                                   index = 3,
                                   verbose = FALSE)
  
  # Test class for hierarchical
  expect_s3_class(ind_hier, "bioregion.site.species.metrics")
  
  # Test print for multiple partitions
  expect_output(print(ind_hier), "Site and species contribution metrics")
  expect_output(print(ind_hier), "Partitions:")
  
  # Test summary for multiple partitions
  expect_output(summary(ind_hier), 
                "Summary of site and species contribution metrics")
  
  # Test summary with n_partitions limit
  expect_output(summary(ind_hier, n_partitions = 1), 
                "Summary of site and species contribution metrics")
  
  # Test str for multiple partitions
  expect_output(str(ind_hier), "bioregion.site.species.metrics object")
  
  # Test with single node type (site only)
  ind_site <- site_species_metrics(bioregionalization = cluinfo,
                                   bioregion_indices = "Rho",
                                   bioregionalization_indices = "P",
                                   data_type = "occurrence",
                                   node_type = "site",
                                   comat = vegemat,
                                   similarity = vegesim,
                                   index = 3,
                                   verbose = FALSE)
  
  expect_s3_class(ind_site, "bioregion.site.species.metrics")
  expect_output(print(ind_site), "Site and species contribution metrics")
  expect_output(str(ind_site), "bioregion.site.species.metrics object")
  expect_output(summary(ind_site), "Summary of site and species contribution metrics")
  
  # Test with single node type (species only)
  ind_species <- site_species_metrics(bioregionalization = cluinfo,
                                      bioregion_indices = "CoreTerms",
                                      bioregionalization_indices = NULL,
                                      data_type = "occurrence",
                                      node_type = "species",
                                      comat = vegemat,
                                      verbose = FALSE)
  
  expect_s3_class(ind_species, "bioregion.site.species.metrics")
  expect_output(print(ind_species), "Site and species contribution metrics")
  expect_output(str(ind_species), "bioregion.site.species.metrics object")
  expect_output(summary(ind_species), "Summary of site and species contribution metrics")
  
  # Test with abundance data type
  ind_abund <- site_species_metrics(bioregionalization = cluinfo,
                                    bioregion_indices = "Rho",
                                    bioregionalization_indices = "P",
                                    data_type = "abundance",
                                    node_type = "site",
                                    comat = vegemat,
                                    verbose = FALSE)
  
  expect_s3_class(ind_abund, "bioregion.site.species.metrics")
  expect_output(print(ind_abund), "Site and species contribution metrics")
  expect_output(print(ind_abund), "Bioregion indices \\(abundance\\):")
  
  # Test with similarity-based metrics only
  ind_sim <- site_species_metrics(bioregionalization = cluinfo,
                                  bioregion_indices = "MeanSim",
                                  bioregionalization_indices = "Silhouette",
                                  data_type = "occurrence",
                                  node_type = "site",
                                  comat = NULL,
                                  similarity = vegesim,
                                  index = 3,
                                  verbose = FALSE)
  
  expect_s3_class(ind_sim, "bioregion.site.species.metrics")
  expect_output(print(ind_sim), "Site and species contribution metrics")
  expect_output(print(ind_sim), "Similarity-based indices:")
  
})



