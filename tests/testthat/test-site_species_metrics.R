# Inputs -----------------------------------------------------------------------
data("fishmat")
data("vegemat")
data("vegedf")

quietly(install_binaries(verbose = FALSE))

comatneg <- vegemat
comatneg[1,1] <- -1

comatwnames1 <- vegemat
rownames(comatwnames1) <- NULL
colnames(comatwnames1) <- NULL
comatwnames2 <- vegemat
rownames(comatwnames2)[1] <- "enistanutsi"
colnames(comatwnames2)[1] <- "enistanutsi"

vegemat_bin <- vegemat
vegemat_bin[vegemat_bin > 0] = 1

vegemat_shuff <- vegemat[sample(dim(vegemat)[1],dim(vegemat)[1]),
                         sample(dim(vegemat)[2],dim(vegemat)[2])]

vegesim <- similarity(vegemat, metric = c("Jaccard", "Simpson", "Sorensen"))
fishsim <- similarity(fishmat, metric = c("Jaccard", "Bray"))

vegesimwnames <- vegesim
vegesimwnames[1,1] <- "einuastnie"

vegesim_shuff <- vegesim[sample(dim(vegesim)[1],dim(vegesim)[1]),]

cluinfo <- netclu_infomap(vegedf, 
                          seed = 1, 
                          bipartite = TRUE)

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

cluinfospe <- site_species_subset(cluinfo, node_type = "species")

cluinfona <- cluinfo
cluinfona$inputs$data_type <- NA
cluinfoocc <- cluinfo
cluinfoocc$inputs$data_type <- "occurrence"
cluinfoabund <- cluinfo
cluinfoabund$inputs$data_type <- "abundance"

clulouv <- netclu_louvain(vegesim)

clualt1 <- clulouv
clualt1$inputs <- clualt1$inputs[-8]
clualt2 <- clulouv
clualt2$inputs <- clualt2$inputs[-9]

# Tests for valid outputs ------------------------------------------------------
test_that("valid output", {
  
  ind <- site_species_metrics(bioregionalization = cluinfo,
                              bioregion_indices = "all",
                              bioregionalization_indices = "all",
                              data_type = "both",
                              node_type = "both",
                              comat = vegemat,
                              similarity = vegesim,
                              index = 3,
                              verbose = FALSE)
  expect_equal(dim(ind$species_bioregions), 
               c(dim(vegemat)[2]*cluinfo$cluster_info[1,2],20))
  expect_equal(colnames(ind$species_bioregions)[1:5], c("Species",
                                                        "Bioregion",
                                                        "n_sb","n_s","n_b"))
  expect_equal(colnames(ind$species_bioregions)[12:14], c("w_sb","w_s","w_b"))
  expect_equal(dim(ind$species_bioregionalization), c(dim(vegemat)[2],3))
  expect_equal(colnames(ind$species_bioregionalization), c("Species",
                                                           "P_occ", "P_abund"))
  expect_equal(dim(ind$site_clusters), 
               c(dim(vegemat)[1]*cluinfo$cluster_info[1,2],20))
  expect_equal(colnames(ind$site_clusters)[1:5], c("Site",
                                                   "Species_cluster",
                                                   "n_gc","n_g","n_c"))
  expect_equal(colnames(ind$site_clusters)[12:14], c("w_gc","w_g","w_c"))
  expect_equal(dim(ind$site_clustering), c(dim(vegemat)[1],3))
  expect_equal(colnames(ind$site_clustering), c("Site",
                                                "P_occ", "P_abund"))
  expect_equal(dim(ind$site_bioregions),
               c(dim(vegemat)[1]*cluinfo$cluster_info[1,2],4))
  expect_equal(colnames(ind$site_bioregions), c("Site",
                                                "Bioregion",
                                                "MeanSim","SdSim"))
  expect_equal(dim(ind$site_bioregionalization), c(dim(vegemat)[1],2))
  expect_equal(colnames(ind$site_bioregionalization), c("Site",
                                                        "Silhouette"))
  
  ind_abund <- site_species_metrics(bioregionalization = cluinfo,
                                    bioregion_indices = "all",
                                    bioregionalization_indices = "all",
                                    data_type = "auto",
                                    node_type = "both",
                                    comat = vegemat,
                                    similarity = vegesim,
                                    index = 3,
                                    verbose = FALSE)
  expect_equal(ind_abund$species_bioregions[,1:5], 
               ind$species_bioregions[,c(1,2,12:14)])
  
  
  ind_shuff <- site_species_metrics(bioregionalization = cluinfo,
                                    bioregion_indices = "all",
                                    bioregionalization_indices = "all",
                                    data_type = "both",
                                    node_type = "both",
                                    comat = vegemat_shuff,
                                    similarity = vegesim_shuff,
                                    index = 3,
                                    verbose = FALSE)
  expect_equal(ind,ind_shuff)
  
  ind <- site_species_metrics(bioregionalization = cluinfo,
                              bioregion_indices = "Rho",
                              bioregionalization_indices = NULL,
                              data_type = "occurrence",
                              node_type = "site",
                              comat = vegemat,
                              similarity = vegesim,
                              index = 3,
                              verbose = FALSE)
  expect_equal(dim(ind$species_bioregions)[2], 3)
  expect_equal(names(ind$species_bioregions)[3], "Rho_occ")
  
  ind <- site_species_metrics(bioregionalization = cluinfo,
                              bioregion_indices = "Rho",
                              bioregionalization_indices = NULL,
                              data_type = "abundance",
                              node_type = "site",
                              comat = vegemat,
                              similarity = vegesim,
                              index = 3,
                              verbose = FALSE)
  expect_equal(dim(ind$species_bioregions)[2], 3)
  expect_equal(names(ind$species_bioregions)[3], "Rho_abund")
  
  ind <- site_species_metrics(bioregionalization = cluinfo,
                              bioregion_indices = NULL,
                              bioregionalization_indices = "P",
                              data_type = "occurrence",
                              node_type = "site",
                              comat = vegemat,
                              similarity = vegesim,
                              index = 3,
                              verbose = FALSE)
  expect_equal(dim(ind$species_bioregionalization)[2], 2)
  expect_equal(names(ind$species_bioregionalization)[2], "P_occ")
  
  ind <- site_species_metrics(bioregionalization = cluinfo,
                              bioregion_indices = "Rho",
                              bioregionalization_indices = "P",
                              data_type = "abundance",
                              node_type = "site",
                              comat = vegemat,
                              similarity = vegesim,
                              index = 3,
                              verbose = FALSE)
  expect_equal(dim(ind$species_bioregions)[2], 3)
  expect_equal(names(ind$species_bioregions)[3], "Rho_abund")
  expect_equal(dim(ind$species_bioregionalization)[2], 2)
  expect_equal(names(ind$species_bioregionalization)[2], "P_abund")
  
  ind <- site_species_metrics(bioregionalization = cluinfo,
                              bioregion_indices = "MeanSim",
                              bioregionalization_indices = NULL,
                              data_type = "occurrence",
                              node_type = "site",
                              comat = vegemat,
                              similarity = vegesim,
                              index = 3,
                              verbose = FALSE)
  expect_equal(dim(ind$site_bioregions)[2], 3)
  expect_equal(names(ind$site_bioregions)[3], "MeanSim")

  ind <- site_species_metrics(bioregionalization = cluinfo,
                              bioregion_indices = "all",
                              bioregionalization_indices = "all",
                              data_type = "occurrence",
                              node_type = "site",
                              comat = vegemat,
                              similarity = vegesim,
                              index = 3,
                              verbose = FALSE)
  expect_equal(ind$site_clusters, NULL)
  expect_equal(ind$site_clustering, NULL)
  
  ind <- site_species_metrics(bioregionalization = cluinfo,
                              bioregion_indices = "MeanSim",
                              bioregionalization_indices = "Silhouette",
                              data_type = "occurrence",
                              node_type = "site",
                              comat = vegemat,
                              similarity = vegesim,
                              index = 3,
                              verbose = FALSE)
  expect_equal(dim(ind$site_bioregions)[2], 3)
  expect_equal(names(ind$site_bioregions)[3], "MeanSim")
  expect_equal(dim(ind$site_bioregionalization)[2], 2)
  expect_equal(names(ind$site_bioregionalization)[2], "Silhouette")
  
  ind <- site_species_metrics(bioregionalization = cluinfo,
                              bioregion_indices = "all",
                              bioregionalization_indices = NULL,
                              data_type = "auto",
                              node_type = "site",
                              comat = vegemat,
                              similarity = vegesim,
                              index = 3,
                              verbose = FALSE)
  expect_equal(ind$species_bioregionalization, NULL)
  expect_equal(ind$site_bioregionalization, NULL)
  expect_equal(ind$site_clusters, NULL)
  expect_equal(ind$site_clustering, NULL)
  
  ind <- site_species_metrics(bioregionalization = cluinfo,
                              bioregion_indices = "Fidelity",
                              bioregionalization_indices = NULL,
                              data_type = "auto",
                              node_type = "site",
                              comat = vegemat,
                              similarity = 1,
                              index = 3,
                              verbose = FALSE)
  expect_equal(ind$site_bioregions, NULL)
  expect_equal(ind$species_bioregionalization, NULL)
  expect_equal(ind$site_bioregionalization, NULL)
  expect_equal(ind$site_clusters, NULL)
  expect_equal(ind$site_clustering, NULL)
  
  ind <- site_species_metrics(bioregionalization = cluhier,
                              bioregion_indices = "all",
                              bioregionalization_indices = "all",
                              data_type = "auto",
                              node_type = "site",
                              comat = fishmat,
                              similarity = fishsim,
                              index = 3,
                              verbose = FALSE)
  expect_equal(length(ind), 3)
  
})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    site_species_metrics(cluinfo,
                         verbose = c(TRUE,1)),
    "verbose must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(cluinfo,
                         verbose = 1),
    "verbose must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(1),
    "bioregionalization must be a bioregion.clusters object.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(clualt1),
    "^bioregionalization is a bioregion.cluster object but it has been altered")
  
  expect_error(
    site_species_metrics(clualt2),
    "^bioregionalization is a bioregion.cluster object but it has been altered")

  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_indices = 1,
                         bioregionalization_indices = NULL),
    "bioregion_indices must be a character.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_indices = c("Rho","nistenu"),
                         bioregionalization_indices = NULL),
    "^One or several bioregion indices chosen are not")
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_indices = NULL,
                         bioregionalization_indices = 1),
    "bioregionalization_indices must be a character.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_indices = NULL,
                         bioregionalization_indices = c("P", "nistenu")),
    "^One or several bioregionalization indices chosen are not")

  quietly(  
  expect_warning(
    site_species_metrics(cluinfo,
                         bioregion_indices = c("Rho", "MeanSim"),
                         bioregionalization_indices = NULL,
                         comat = vegemat),
    "^Site-to-bioregion indices ")
  )

  expect_warning(
    site_species_metrics(cluinfo,
                         bioregion_indices = c("Rho", "MeanSim"),
                         bioregionalization_indices = NULL,
                         comat = NULL,
                         similarity = vegesim),
    "^Species-to-bioregion indices ") 
  
  quietly(
  expect_warning(
    site_species_metrics(cluinfo,
                         bioregion_indices = NULL,
                         bioregionalization_indices = "all",
                         comat = vegemat),
    "^Site-to-bioregionalization indices ")
  )
  
  expect_warning(
    site_species_metrics(cluinfo,
                         bioregion_indices = NULL,
                         bioregionalization_indices = "all",
                         comat = NULL,
                         similarity = vegesim),
    "^Species-to-bioregionalization indices ")
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_indices = NULL,
                         bioregionalization_indices = NULL,
                         comat = vegemat),
    "At least one type of indices with the appropriate inputs should be specified.",
    fixed = TRUE)

  expect_error(
    site_species_metrics(cluinfo,
                         node_type = c(TRUE,1),
                         comat = vegemat),
    "node_type must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(cluinfo,
                         node_type = 1,
                         comat = vegemat),
    "node_type must be a character.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(cluinfo,
                         node_type = "iueui",
                         comat = vegemat),
    "^Please choose node_type from the following")
  
  expect_error(
    site_species_metrics(cluinfospe,
                         bioregion_indices = "all",
                         bioregionalization_indices = NULL,
                         node_type = "site",
                         comat = vegemat,
                         similarity = vegesim),
    "^Site-to-bioregion metrics")
  
  expect_error(
    site_species_metrics(cluinfospe,
                         bioregion_indices = "all",
                         bioregionalization_indices = NULL,
                         node_type = "both",
                         comat = vegemat,
                         similarity = vegesim),
    "^Site-to-bioregion metrics")
  
  expect_error(
    site_species_metrics(cluinfospe,
                         bioregion_indices = "MeanSim",
                         bioregionalization_indices = NULL,
                         node_type = "site",                        
                         comat = vegemat,
                         similarity = vegesim),
    "^Site-to-bioregion metrics")
  
  expect_error(
    site_species_metrics(cluinfospe,
                         bioregion_indices = NULL,
                         bioregionalization_indices = "Silhouette",
                         node_type = "site",
                         comat = vegemat,
                         similarity = vegesim),
    "^Site-to-bioregion metrics")
  
  
  expect_error(
    site_species_metrics(clulouv,
                         node_type = "species",
                         comat = vegemat,
                         similarity = vegesim),
    "^Site-to-species_cluster metrics")
  
  expect_error(
    site_species_metrics(clulouv,
                         node_type = "both",
                         comat = vegemat,
                         similarity = vegesim),
    "^Site-to-species_cluster metrics")
  
  expect_error(
    site_species_metrics(cluinfospe,
                         bioregion_indices = "CoreTerms",
                         bioregionalization_indices = NULL,
                         node_type = "site",
                         comat = vegemat,
                         similarity = vegesim),
    "^Species-to-bioregion metrics")
  
  expect_error(
    site_species_metrics(cluinfospe,
                         bioregion_indices = NULL,
                         bioregionalization_indices = "P",
                         node_type = "site",
                         comat = vegemat,
                         similarity = vegesim),
    "^Species-to-bioregion metrics")
  
  expect_error(
    site_species_metrics(cluinfospe,
                         bioregion_indices = "CoreTerms",
                         bioregionalization_indices = NULL,
                         node_type = "both",
                         comat = vegemat,
                         similarity = vegesim),
    "^Species-to-bioregion metrics")
  
  expect_error(
    site_species_metrics(cluinfospe,
                         bioregion_indices = NULL,
                         bioregionalization_indices = "P",
                         node_type = "both",
                         comat = vegemat,
                         similarity = vegesim),
    "^Species-to-bioregion metrics")
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_indices = NULL,
                         bioregionalization_indices = "P",
                         comat=1),
    "comat must be a matrix.", 
    fixed = TRUE)

  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_indices = NULL,
                         bioregionalization_indices = "P",
                         comat = comatneg),
    "Negative value(s) detected in comat!", 
    fixed = TRUE)
  
  # expect_message(
  #   site_species_metrics(cluinfo,
  #                        comat = vegemat,
  #                        similarity = vegesim),
  #   "comat is based on abundance data.")
  # 
  # expect_message(
  #   site_species_metrics(cluinfo,
  #                        data_type = "occurrence",
  #                        comat = vegemat_bin,
  #                        similarity = vegesim),
  #   "comat is based on occurrence data.")
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_indices = NULL,
                         bioregionalization_indices = "P",
                         comat = vegemat,
                         data_type = c(TRUE,1)),
    "data_type must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_indices = NULL,
                         bioregionalization_indices = "P",
                         comat = vegemat,
                         data_type = 1),
    "data_type must be a character.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_indices = NULL,
                         bioregionalization_indices = "P",
                         comat = vegemat,
                         data_type = "iueui"),
    "^Please choose data_type from the following")
  
  expect_message(
    site_species_metrics(cluinfona,
                         comat = vegemat_bin,
                         bioregion_indices = "Rho",
                         bioregionalization_indices = NULL,
                         data_type = "auto"),
    "^No data type detected in bioregionalization and comat is based on occurence data")
  
  expect_message(
    site_species_metrics(cluinfona,
                         comat = vegemat,
                         bioregion_indices = "Rho",
                         bioregionalization_indices = NULL,
                         data_type = "auto"),
    "^No data type detected in bioregionalization and comat is based on abundance data")
  
  expect_message(
    site_species_metrics(cluinfoocc,
                         comat = vegemat_bin,
                         bioregion_indices = "Rho",
                         bioregionalization_indices = NULL,
                         data_type = "auto"),
    "^The bioregionalization is based on occurence data and comat is based on occurence data")
  
  expect_message(
    site_species_metrics(cluinfoocc,
                         comat = vegemat,
                         bioregion_indices = "Rho",
                         bioregionalization_indices = NULL,
                         data_type = "auto"),
    "^The bioregionalization is based on occurence data but note that even if comat is based on abundance")
  
  expect_message(
    site_species_metrics(cluinfoabund,
                         comat = vegemat_bin,
                         bioregion_indices = "Rho",
                         bioregionalization_indices = NULL,
                         data_type = "auto"),
    "^The bioregionalization is based on abundance data but comat is based on occurence data so occurrence")
  
  expect_message(
    site_species_metrics(cluinfoabund,
                         comat = vegemat,
                         bioregion_indices = "Rho",
                         bioregionalization_indices = NULL,
                         data_type = "auto"),
    "^The bioregionalization is based on abundance data and comat is based on abundance")
  
  #expect_warning(
  #  site_species_metrics(cluinfo,
  #                       comat = vegemat_bin,
  #                       bioregion_indices = "Rho",
  #                       bioregionalization_indices = NULL,
  #                       data_type = "abundance"),
  #  "^comat is based on occurence data so abundance-based indices won't be computed!")
  
  expect_error(
    site_species_metrics(cluinfo,
                         comat = comatwnames1,
                         verbose = FALSE),
    "Site ID in bioregionalization and comat do not match!", 
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(cluinfo,
                         comat = comatwnames2,
                         verbose = FALSE),
    "Site ID in bioregionalization and comat do not match!", 
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(cluinfo,
                         node_type = "species",
                         comat = comatwnames1,
                         verbose = FALSE),
    "Species ID in bioregionalization and comat do not match!", 
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(cluinfo,
                         node_type = "species",
                         comat = comatwnames2,
                         verbose = FALSE),
    "Species ID in bioregionalization and comat do not match!", 
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_indices = "all",
                         bioregionalization_indices = NULL,
                         comat = vegemat,
                         similarity = "1",
                         verbose = FALSE),
    "^similarity should be a bioregion.pairwise object created by")
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_indices = "MeanSim",
                         bioregionalization_indices = NULL,
                         comat = NULL,
                         similarity = "1"),
    "^similarity should be a bioregion.pairwise object created by")
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_indices = "MeanSim",
                         bioregionalization_indices = NULL,
                         comat = NULL,
                         similarity = vegesim,
                         index = c("1", 1)),
    "index must be of length 1.")
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_indices = "MeanSim",
                         bioregionalization_indices = NULL,
                         comat = NULL,
                         similarity = vegesim,
                         index = "zz"),
    "^If index is a character, it should be a column")
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_indices = "MeanSim",
                         bioregionalization_indices = NULL,
                         comat = NULL,
                         similarity = vegesim,
                         index = "Site1"),
    "^If index is a character, it should be a column")
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_indices = "MeanSim",
                         bioregionalization_indices = NULL,
                         comat = NULL,
                         similarity = vegesim,
                         index = 0.1),
    "If index is numeric, it should be an integer.")
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_indices = "MeanSim",
                         bioregionalization_indices = NULL,
                         comat = NULL,
                         similarity = vegesim,
                         index = 1),
    "index should be strictly higher than 2.")
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_indices = "MeanSim",
                         bioregionalization_indices = NULL,
                         comat = NULL,
                         similarity = vegesim,
                         index = 10),
    "index should be lower or equal to 5.")
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_indices = "MeanSim",
                         bioregionalization_indices = NULL,
                         comat = NULL,
                         similarity = vegesimwnames),
    "Site ID in bioregionalization and similarity do not match!")

})

# Tests for generic functions --------------------------------------------------
test_that("generic functions work correctly for bioregion.site.species.metrics", {
  
  # Create a site_species_metrics object for testing
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
