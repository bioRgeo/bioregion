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
                              bioregion_metrics = "all",
                              bioregionalization_metrics = "all",
                              data_type = "both",
                              cluster_on = "both",
                              comat = vegemat,
                              similarity = vegesim,
                              include_cluster = TRUE,
                              index = 3,
                              verbose = FALSE)
  expect_equal(dim(ind$species_bioregions), 
               c(dim(vegemat)[2]*cluinfo$cluster_info[1,2], 20))
  expect_equal(colnames(ind$species_bioregions)[1:5], c("Species",
                                                        "Bioregion",
                                                        "n_sb","n_s","n_b"))
  expect_equal(colnames(ind$species_bioregions)[12:14], c("w_sb","w_s","w_b"))
  expect_equal(dim(ind$species_bioregionalization), c(dim(vegemat)[2],3))
  expect_equal(colnames(ind$species_bioregionalization), c("Species",
                                                           "P_occ", "P_abund"))
  expect_equal(dim(ind$site_chorotypes), 
               c(dim(vegemat)[1]*cluinfo$cluster_info[1,2],20))
  expect_equal(colnames(ind$site_chorotypes)[1:5], c("Site",
                                                     "Chorotypes",
                                                     "n_gc","n_g","n_c"))
  expect_equal(colnames(ind$site_chorotypes)[12:14], c("w_gc","w_g","w_c"))
  expect_equal(dim(ind$site_chorological), c(dim(vegemat)[1],3))
  expect_equal(colnames(ind$site_chorological), c("Site",
                                                "P_occ", "P_abund"))
  expect_equal(dim(ind$site_bioregions),
               c(dim(vegemat)[1]*cluinfo$cluster_info[1,2],8))
  expect_equal(colnames(ind$site_bioregions), c("Site",
                                                "Bioregion",
                                                "Assigned",
                                                "Richness",
                                                "Rich_Endemics",
                                                "Prop_Endemics",
                                                "MeanSim","SdSim"))
  expect_equal(dim(ind$site_bioregionalization), c(dim(vegemat)[1],2))
  expect_equal(colnames(ind$site_bioregionalization), c("Site",
                                                        "Silhouette"))
  
  ind_abund <- site_species_metrics(bioregionalization = cluinfo,
                                    bioregion_metrics = "all",
                                    bioregionalization_metrics = "all",
                                    data_type = "auto",
                                    cluster_on = "both",
                                    comat = vegemat,
                                    similarity = vegesim,
                                    index = 3,
                                    verbose = FALSE)
  expect_equal(ind_abund$species_bioregions[,1:5], 
               ind$species_bioregions[,c(1,2,12:14)])
  
  
  ind_shuff <- site_species_metrics(bioregionalization = cluinfo,
                                    bioregion_metrics = "all",
                                    bioregionalization_metrics = "all",
                                    data_type = "both",
                                    cluster_on = "both",
                                    comat = vegemat_shuff,
                                    similarity = vegesim_shuff,
                                    include_cluster = TRUE,
                                    index = 3,
                                    verbose = FALSE)
  expect_equal(ind,ind_shuff)
  
  ind <- site_species_metrics(bioregionalization = cluinfo,
                              bioregion_metrics = "Rho",
                              bioregionalization_metrics = NULL,
                              data_type = "occurrence",
                              cluster_on = "site",
                              comat = vegemat,
                              similarity = vegesim,
                              index = 3,
                              verbose = FALSE)
  expect_equal(dim(ind$species_bioregions)[2], 3)
  expect_equal(names(ind$species_bioregions)[3], "Rho_occ")
  
  ind <- site_species_metrics(bioregionalization = cluinfo,
                              bioregion_metrics = "Rho",
                              bioregionalization_metrics = NULL,
                              data_type = "abundance",
                              cluster_on = "site",
                              comat = vegemat,
                              similarity = vegesim,
                              index = 3,
                              verbose = FALSE)
  expect_equal(dim(ind$species_bioregions)[2], 3)
  expect_equal(names(ind$species_bioregions)[3], "Rho_abund")
  
  ind <- site_species_metrics(bioregionalization = cluinfo,
                              bioregion_metrics = NULL,
                              bioregionalization_metrics = "P",
                              data_type = "occurrence",
                              cluster_on = "site",
                              comat = vegemat,
                              similarity = vegesim,
                              index = 3,
                              verbose = FALSE)
  expect_equal(dim(ind$species_bioregionalization)[2], 2)
  expect_equal(names(ind$species_bioregionalization)[2], "P_occ")
  
  ind <- site_species_metrics(bioregionalization = cluinfo,
                              bioregion_metrics = "Rho",
                              bioregionalization_metrics = "P",
                              data_type = "abundance",
                              cluster_on = "site",
                              comat = vegemat,
                              similarity = vegesim,
                              index = 3,
                              verbose = FALSE)
  expect_equal(dim(ind$species_bioregions)[2], 3)
  expect_equal(names(ind$species_bioregions)[3], "Rho_abund")
  expect_equal(dim(ind$species_bioregionalization)[2], 2)
  expect_equal(names(ind$species_bioregionalization)[2], "P_abund")
  
  ind <- site_species_metrics(bioregionalization = cluinfo,
                              bioregion_metrics = "MeanSim",
                              bioregionalization_metrics = NULL,
                              data_type = "occurrence",
                              cluster_on = "site",
                              comat = vegemat,
                              similarity = vegesim,
                              index = 3,
                              verbose = FALSE)
  expect_equal(dim(ind$site_bioregions)[2], 3)
  expect_equal(names(ind$site_bioregions)[3], "MeanSim")

  ind <- site_species_metrics(bioregionalization = cluinfo,
                              bioregion_metrics = "all",
                              bioregionalization_metrics = "all",
                              data_type = "occurrence",
                              cluster_on = "site",
                              comat = vegemat,
                              similarity = vegesim,
                              index = 3,
                              verbose = FALSE)
  expect_equal(ind$site_chorotypes, NULL)
  expect_equal(ind$site_chorological, NULL)
  
  ind <- site_species_metrics(bioregionalization = cluinfo,
                              bioregion_metrics = "MeanSim",
                              bioregionalization_metrics = "Silhouette",
                              data_type = "occurrence",
                              cluster_on = "site",
                              comat = vegemat,
                              similarity = vegesim,
                              index = 3,
                              verbose = FALSE)
  expect_equal(dim(ind$site_bioregions)[2], 3)
  expect_equal(names(ind$site_bioregions)[3], "MeanSim")
  expect_equal(dim(ind$site_bioregionalization)[2], 2)
  expect_equal(names(ind$site_bioregionalization)[2], "Silhouette")
  
  ind <- site_species_metrics(bioregionalization = cluinfo,
                              bioregion_metrics = "all",
                              bioregionalization_metrics = NULL,
                              data_type = "auto",
                              cluster_on = "site",
                              comat = vegemat,
                              similarity = vegesim,
                              index = 3,
                              verbose = FALSE)
  expect_equal(ind$species_bioregionalization, NULL)
  expect_equal(ind$site_bioregionalization, NULL)
  expect_equal(ind$site_chorotypes, NULL)
  expect_equal(ind$site_chorological, NULL)
  
  ind <- site_species_metrics(bioregionalization = cluinfo,
                              bioregion_metrics = "Fidelity",
                              bioregionalization_metrics = NULL,
                              data_type = "auto",
                              cluster_on = "site",
                              comat = vegemat,
                              similarity = 1,
                              index = 3,
                              verbose = FALSE)
  expect_equal(ind$site_bioregions, NULL)
  expect_equal(ind$species_bioregionalization, NULL)
  expect_equal(ind$site_bioregionalization, NULL)
  expect_equal(ind$site_chorotypes, NULL)
  expect_equal(ind$site_chorological, NULL)
  
  ind <- site_species_metrics(bioregionalization = cluhier,
                              bioregion_metrics = "all",
                              bioregionalization_metrics = "all",
                              data_type = "auto",
                              cluster_on = "site",
                              comat = fishmat,
                              similarity = fishsim,
                              index = 3,
                              verbose = FALSE)
  expect_equal(length(ind), 3)
  
  ind <- site_species_metrics(bioregionalization = cluinfo,
                              bioregion_metrics = "MeanSim",
                              bioregionalization_metrics = NULL,
                              data_type = "both",
                              cluster_on = "both",
                              comat = vegemat,
                              similarity = vegesim,
                              include_cluster = TRUE,
                              index = 3,
                              verbose = FALSE)
  expect_equal(dim(ind$site_bioregions)[2], 4)
  
  ind <- site_species_metrics(bioregionalization = cluinfo,
                              bioregion_metrics = "MeanSim",
                              bioregionalization_metrics = NULL,
                              data_type = "both",
                              cluster_on = "both",
                              comat = vegemat,
                              similarity = vegesim,
                              include_cluster = FALSE,
                              index = 3,
                              verbose = FALSE)
  expect_equal(dim(ind$site_bioregions)[2], 3)
  
  ind_base <- site_species_metrics(bioregionalization = cluinfo,
                                   bioregion_metrics = c("Richness","MeanSim"),
                                   bioregionalization_metrics = NULL,
                                   data_type = "both",
                                   cluster_on = "both",
                                   comat = vegemat,
                                   similarity = vegesim,
                                   include_cluster = FALSE,
                                   index = 3,
                                   verbose = FALSE)
  expect_equal(dim(ind_base$site_bioregions)[2], 4)
  expect_equal(identical(ind$site_bioregions, 
                         ind_base$site_bioregions[,-3]), 
               TRUE)
  
  
  
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
    site_species_metrics(cluinfo,
                         include_cluster = c(TRUE,1)),
    "include_cluster must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(cluinfo,
                         include_cluster = 1),
    "include_cluster must be a boolean.",
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
                         bioregion_metrics = 1,
                         bioregionalization_metrics = NULL),
    "bioregion_metrics must be a character.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_metrics = c("Rho","nistenu"),
                         bioregionalization_metrics = NULL),
    "^One or several bioregion metrics chosen are not")
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_metrics = NULL,
                         bioregionalization_metrics = 1),
    "bioregionalization_metrics must be a character.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_metrics = NULL,
                         bioregionalization_metrics = c("P", "nistenu")),
    "^One or several bioregionalization metrics chosen are not")
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_metrics = NULL,
                         bioregionalization_metrics = NULL,
                         comat = vegemat),
    "At least one type of metrics with the appropriate inputs should be specified.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_metrics = "all",
                         bioregionalization_metrics = "all",
                         comat = NULL,
                         similarity = NULL),
    "At least comat or similarity should be provided.", 
    fixed = TRUE)
  
  expect_warning(
    site_species_metrics(cluinfo,
                         bioregion_metrics = c("Rho", "MeanSim"),
                         bioregionalization_metrics = NULL,
                         comat = NULL,
                         similarity = vegesim),
    "^Some metrics ") 
  
  expect_warning(
    site_species_metrics(cluinfo,
                         bioregion_metrics = NULL,
                         bioregionalization_metrics = "all",
                         comat = NULL,
                         similarity = vegesim),
    "^Some metrics ")

  expect_warning(
    site_species_metrics(cluinfo,
                         bioregion_metrics = c("Rho", "MeanSim"),
                         bioregionalization_metrics = NULL,
                         data_type = "abundance",
                         comat = vegemat,
                         similarity = NULL),
    "^Some metrics ")

  expect_warning(
    site_species_metrics(cluinfo,
                         bioregion_metrics = NULL,
                         bioregionalization_metrics = "all",
                         data_type = "abundance",
                         comat = vegemat,
                         similarity = NULL),
    "^Some metrics ")
  
  expect_error(
    expect_warning(
      site_species_metrics(cluinfo,
                           bioregion_metrics = "MeanSim",
                           bioregionalization_metrics = NULL,
                           comat = vegemat,
                           similarity = NULL),
      "^Some metrics "),
    "At least one type of metrics with the appropriate inputs should be specified.",
    fixed = TRUE)
  
  expect_error(
    expect_warning(
      site_species_metrics(cluinfo,
                           bioregion_metrics = NULL,
                           bioregionalization_metrics = "Silhouette",
                           comat = vegemat,
                           similarity = NULL),
      "^Some metrics "),
    "At least one type of metrics with the appropriate inputs should be specified.",
    fixed = TRUE)
  
  expect_error(
    expect_warning(
      site_species_metrics(cluinfo,
                           bioregion_metrics = "Rho",
                           bioregionalization_metrics = NULL,
                           comat = NULL,
                           similarity = vegesim),
      "^Some metrics "),
    "At least one type of metrics with the appropriate inputs should be specified.",
    fixed = TRUE)
  
  expect_error(
    expect_warning(
      site_species_metrics(cluinfo,
                           bioregion_metrics = "Richness",
                           bioregionalization_metrics = NULL,
                           comat = NULL,
                           similarity = vegesim),
      "^Some metrics "),
    "At least one type of metrics with the appropriate inputs should be specified.",
    fixed = TRUE)
  
  expect_error(
    expect_warning(
      site_species_metrics(cluinfo,
                           bioregion_metrics = NULL,
                           bioregionalization_metrics = "P",
                           comat = NULL,
                           similarity = vegesim),
      "^Some metrics "),
    "At least one type of metrics with the appropriate inputs should be specified.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(cluinfo,
                         cluster_on = c(TRUE,1),
                         comat = vegemat),
    "cluster_on must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(cluinfo,
                         cluster_on = 1,
                         comat = vegemat),
    "cluster_on must be a character.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(cluinfo,
                         cluster_on = "iueui",
                         comat = vegemat),
    "^Please choose cluster_on from the following")
  
  expect_error(
    site_species_metrics(cluinfospe,
                         bioregion_metrics = "CoreTerms",
                         bioregionalization_metrics = NULL,
                         cluster_on = "site",
                         comat = vegemat,
                         similarity = vegesim),
    "^Species-to-bioregions/bioregionalization metrics")
  
  expect_error(
    site_species_metrics(cluinfospe,
                         bioregion_metrics = NULL,
                         bioregionalization_metrics = "P",
                         cluster_on = "site",
                         comat = vegemat,
                         similarity = vegesim),
    "^Species-to-bioregions/bioregionalization metrics")
  
  expect_error(
    site_species_metrics(cluinfospe,
                         bioregion_metrics = "CoreTerms",
                         bioregionalization_metrics = NULL,
                         cluster_on = "both",
                         comat = vegemat,
                         similarity = vegesim),
    "^Species-to-bioregions/bioregionalization metrics")
  
  expect_error(
    site_species_metrics(cluinfospe,
                         bioregion_metrics = NULL,
                         bioregionalization_metrics = "P",
                         cluster_on = "both",
                         comat = vegemat,
                         similarity = vegesim),
    "^Species-to-bioregions/bioregionalization metrics")
  
  expect_error(
    site_species_metrics(cluinfospe,
                         bioregion_metrics = "MeanSim",
                         bioregionalization_metrics = NULL,
                         cluster_on = "site",
                         comat = vegemat,
                         similarity = vegesim),
    "^Site-to-bioregions/bioregionalization metrics")
  
  expect_error(
    site_species_metrics(cluinfospe,
                         bioregion_metrics = NULL,
                         bioregionalization_metrics = "Silhouette",
                         cluster_on = "site",
                         comat = vegemat,
                         similarity = vegesim),
    "^Site-to-bioregions/bioregionalization metrics")
  
  expect_error(
    site_species_metrics(cluinfospe,
                         bioregion_metrics = "MeanSim",
                         bioregionalization_metrics = NULL,
                         cluster_on = "both",                        
                         comat = vegemat,
                         similarity = vegesim),
    "^Site-to-bioregions/bioregionalization metrics")
  
  expect_error(
    site_species_metrics(cluinfospe,
                         bioregion_metrics = NULL,
                         bioregionalization_metrics = "Silhouette",
                         cluster_on = "both",
                         comat = vegemat,
                         similarity = vegesim),
    "^Site-to-bioregions/bioregionalization metrics")
  
  expect_error(
    site_species_metrics(cluinfospe,
                         bioregion_metrics = "all",
                         bioregionalization_metrics = NULL,
                         cluster_on = "site",
                         comat = vegemat,
                         similarity = vegesim),
    "^Species/Site-to-bioregions/bioregionalization metrics")
  
  expect_error(
    site_species_metrics(cluinfospe,
                         bioregion_metrics = NULL,
                         bioregionalization_metrics = "all",
                         cluster_on = "site",
                         comat = vegemat,
                         similarity = vegesim),
    "^Species/Site-to-bioregions/bioregionalization metrics")
  
  expect_error(
    site_species_metrics(cluinfospe,
                         bioregion_metrics = "all",
                         bioregionalization_metrics = NULL,
                         cluster_on = "both",                        
                         comat = vegemat,
                         similarity = vegesim),
    "^Species/Site-to-bioregions/bioregionalization metrics")
  
  expect_error(
    site_species_metrics(cluinfospe,
                         bioregion_metrics = NULL,
                         bioregionalization_metrics = "all",
                         cluster_on = "both",
                         comat = vegemat,
                         similarity = vegesim),
    "^Species/Site-to-bioregions/bioregionalization metrics")
  
  expect_error(
    site_species_metrics(clulouv,
                         bioregion_metrics = "CoreTerms",
                         bioregionalization_metrics = NULL,
                         cluster_on = "species",
                         comat = vegemat,
                         similarity = vegesim),
    "^Site-to-chorotypes/chorological metrics")
  
  expect_error(
    site_species_metrics(clulouv,
                         bioregion_metrics = NULL,
                         bioregionalization_metrics = "P",
                         cluster_on = "species",
                         comat = vegemat,
                         similarity = vegesim),
    "^Site-to-chorotypes/chorological metrics")
  
  expect_error(
    site_species_metrics(clulouv,
                         bioregion_metrics = "CoreTerms",
                         bioregionalization_metrics = NULL,
                         cluster_on = "both",
                         comat = vegemat,
                         similarity = vegesim),
    "^Site-to-chorotypes/chorological metrics")
  
  expect_error(
    site_species_metrics(clulouv,
                         bioregion_metrics = NULL,
                         bioregionalization_metrics = "P",
                         cluster_on = "both",
                         comat = vegemat,
                         similarity = vegesim),
    "^Site-to-chorotypes/chorological metrics")

  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_metrics = NULL,
                         bioregionalization_metrics = "P",
                         comat=1),
    "comat must be a matrix.", 
    fixed = TRUE)

  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_metrics = NULL,
                         bioregionalization_metrics = "P",
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
                         bioregion_metrics = NULL,
                         bioregionalization_metrics = "P",
                         comat = vegemat,
                         data_type = c(TRUE,1)),
    "data_type must be of length 1.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_metrics = NULL,
                         bioregionalization_metrics = "P",
                         comat = vegemat,
                         data_type = 1),
    "data_type must be a character.",
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_metrics = NULL,
                         bioregionalization_metrics = "P",
                         comat = vegemat,
                         data_type = "iueui"),
    "^Please choose data_type from the following")
  
  expect_message(
    site_species_metrics(cluinfona,
                         comat = vegemat_bin,
                         bioregion_metrics = "Rho",
                         bioregionalization_metrics = NULL,
                         data_type = "auto"),
    "^No data type detected in bioregionalization and comat is based on occurence data")
  
  expect_message(
    site_species_metrics(cluinfona,
                         comat = vegemat,
                         bioregion_metrics = "Rho",
                         bioregionalization_metrics = NULL,
                         data_type = "auto"),
    "^No data type detected in bioregionalization and comat is based on abundance data")
  
  expect_message(
    site_species_metrics(cluinfoocc,
                         comat = vegemat_bin,
                         bioregion_metrics = "Rho",
                         bioregionalization_metrics = NULL,
                         data_type = "auto"),
    "^The bioregionalization is based on occurence data and comat is based on occurence data")
  
  expect_message(
    site_species_metrics(cluinfoocc,
                         comat = vegemat,
                         bioregion_metrics = "Rho",
                         bioregionalization_metrics = NULL,
                         data_type = "auto"),
    "^The bioregionalization is based on occurence data but note that even if comat is based on abundance")
  
  expect_message(
    site_species_metrics(cluinfoabund,
                         comat = vegemat_bin,
                         bioregion_metrics = "Rho",
                         bioregionalization_metrics = NULL,
                         data_type = "auto"),
    "^The bioregionalization is based on abundance data but comat is based on occurence data so occurrence")
  
  expect_message(
    site_species_metrics(cluinfoabund,
                         comat = vegemat,
                         bioregion_metrics = "Rho",
                         bioregionalization_metrics = NULL,
                         data_type = "auto"),
    "^The bioregionalization is based on abundance data and comat is based on abundance")
  
  #expect_warning(
  #  site_species_metrics(cluinfo,
  #                       comat = vegemat_bin,
  #                       bioregion_metrics = "Rho",
  #                       bioregionalization_metrics = NULL,
  #                       data_type = "abundance"),
  #  "^comat is based on occurence data so abundance-based metrics won't be computed!")
  
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
                         cluster_on = "species",
                         comat = comatwnames1,
                         verbose = FALSE),
    "Species ID in bioregionalization and comat do not match!", 
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(cluinfo,
                         cluster_on = "species",
                         comat = comatwnames2,
                         verbose = FALSE),
    "Species ID in bioregionalization and comat do not match!", 
    fixed = TRUE)
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_metrics = "all",
                         bioregionalization_metrics = NULL,
                         comat = vegemat,
                         similarity = "1",
                         verbose = FALSE),
    "^similarity should be a bioregion.pairwise object created by")
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_metrics = "MeanSim",
                         bioregionalization_metrics = NULL,
                         comat = NULL,
                         similarity = "1"),
    "^similarity should be a bioregion.pairwise object created by")
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_metrics = "MeanSim",
                         bioregionalization_metrics = NULL,
                         comat = NULL,
                         similarity = vegesim,
                         index = c("1", 1)),
    "index must be of length 1.")
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_metrics = "MeanSim",
                         bioregionalization_metrics = NULL,
                         comat = NULL,
                         similarity = vegesim,
                         index = "zz"),
    "^If index is a character, it should be a column")
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_metrics = "MeanSim",
                         bioregionalization_metrics = NULL,
                         comat = NULL,
                         similarity = vegesim,
                         index = "Site1"),
    "^If index is a character, it should be a column")
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_metrics = "MeanSim",
                         bioregionalization_metrics = NULL,
                         comat = NULL,
                         similarity = vegesim,
                         index = 0.1),
    "If index is numeric, it should be an integer.")
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_metrics = "MeanSim",
                         bioregionalization_metrics = NULL,
                         comat = NULL,
                         similarity = vegesim,
                         index = 1),
    "index should be strictly higher than 2.")
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_metrics = "MeanSim",
                         bioregionalization_metrics = NULL,
                         comat = NULL,
                         similarity = vegesim,
                         index = 10),
    "index should be lower or equal to 5.")
  
  expect_error(
    site_species_metrics(cluinfo,
                         bioregion_metrics = "MeanSim",
                         bioregionalization_metrics = NULL,
                         comat = NULL,
                         similarity = vegesimwnames),
    "Site ID in bioregionalization and similarity do not match!")

})
