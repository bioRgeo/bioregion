#' Calculate contribution metrics of sites and species to each bioregion
#'
#' This function computes metrics to assess the contribution of a given
#' species or site to each bioregion.
#'
#' @param bioregionalization A `bioregion.clusters` object.
#'
#' @param bioregion_metrics A `character` vector or a single `character` string
#' specifying the metrics to compute. Several metrics belonging to different 
#' categories are available: species to bioregions (`"Specificity"`, 
#' `"NSpecificity"`, `"Fidelity"`, `"NFidelity"`, `"IndVal"`, `"NIndVal"`, 
#' `"Rho"` and `"CoreTerms"`), site to bioregions (`"Richness"`, 
#' `"Rich_Endemics"`, `"Prop_Endemics"`, `"MeanSim"` and `"SdSim"`), or site 
#' to chorotypes (`"Specificity"`, `"NSpecificity"`, `"Fidelity"`, `"NFidelity"`, 
#' `"IndVal"`, `"NIndVal"`, `"Rho"` and `"CoreTerms"`)
#' when a cluster has been assigned to the species in 
#' `bioregionalization` (see Details). If `"all"` is specified, all 
#' metrics will be calculated.
#' 
#' @param bioregionalization_metrics A `character` vector or a single `character` 
#' string specifying the bioregionalization metrics to compute. Some aggregated 
#' metrics (such as the participation coefficient `"P"` or Silhouette index 
#' `"Silhouette"`) can be derived from the `bioregion_metrics` (see Details). 
#' If `"all"` is specified , all bioregionalization metrics will be calculated.
#' 
#' @param data_type A `character` string specifying which type of data should 
#' be considered to compute the related metrics (`"Specificity"`, 
#' `"NSpecificity"``"Fidelity"`, `"NFidelity"`, `"IndVal"`, `"NIndVal"`, 
#' `"Rho"`, `"CoreTerms"`): occurrences or abundances. By default (`"auto"`), 
#' the type of data is inferred from `bioregionalization` and/or the provided 
#' co-occurrence matrix (argument `comat`). Other available options are 
#' `"occurence"`, `"abundance"` or `"both"` (see Details).
#' 
#' @param cluster_on A `character` string specifying whether the related
#'  metrics (`"Specificity"`, `"NSpecificity"`, `"Fidelity"`, `"NFidelity"`, 
#'  `"IndVal"`, `"NIndVal"`, `"Rho"` and `"CoreTerms"`) 
#'  should be computed as contributions from species to bioregions (clusters 
#'  based on `"site"` by default), from site to chorotypes (clusters based on 
#'  `"species"`) or for `"both"`. The latter type of contribution is only 
#'  available if a cluster has been assigned to the species in 
#'  `bioregionalization` (see Details).
#' 
#' @param comat A co-occurrence `matrix` with sites as rows and species as
#' columns.
#' 
#' @param similarity The output object from [similarity()] or
#' [dissimilarity_to_similarity()].
#' 
#' @param include_cluster A `boolean` indicating whether to include an
#' additional column `Assigned` in the site-to-bioregions output, indicating for
#' each site whether the considered bioregion is its bioregion (`FALSE` by 
#' default).
#' 
#' @param index The name or number of the column to use as similarity. 
#' By default, the third column name of `similarity` is used.
#' 
#' @param verbose A `boolean` indicating whether to 
#' display progress messages. Set to `FALSE` to suppress these messages.
#' 
#' @return 
#' A `list` containing between one and six elements (listed below), depending on 
#' the selected metrics (`bioregion_metrics` and/or 
#' `bioregionalization_metrics`) and the type of cluster 
#' (bioregion when `cluster_on = "site"` and/or 
#' chorotypes when `cluster_on = "species"`).
#' \itemize{
#'   \item{**species_bioregions**: A `data.frame` containing the 
#'   species-to-bioregions metric(s) based on `comat`.}
#'   \item{**species_bioregionalization**: A `data.frame` containing the  
#'   species-to-bioregionalization metric(s) based on `comat`.}
#'   \item{**site_bioregions**: A `data.frame` containing the 
#'   site-to-bioregions metric(s) based on `comat` or `similarity`.}
#'   \item{**site_bioregionalization**: A `data.frame` containing the 
#'   site-to-bioregionalization metric(s) based on `comat` or `similarity`.}
#'   \item{**site_chorotypes**: A `data.frame` containing the  
#'   site-to-chorotypes metric(s) based on `comat`.}
#'   \item{**site_chorological**: A `data.frame` containing the 
#'   site-to-chorological classification metric(s) based on `comat`.}
#' }
#'
#' Note that if `bioregionalization` contains more than one partition 
#' (i.e., if `dim(bioregionalization$clusters) > 2`), a list of lists will be 
#' returned, with one sublist per partition.
#' 
#' @details
#' The **first type** of metrics provided by this function is based on the 
#' contribution of each species to a given bioregion (`cluster_on = "site"` or 
#' `cluster_on = "both"`). This is calculated from the bioregion assigned to 
#' each site (as defined in `bioregionalization`) and 
#' an associated site–species co-occurrence matrix `comat` (a strict match 
#' between site IDs is required). 
#' 
#' Occurrence-based metrics (`data_type = "occurrence"` or `data_type = "both"`) 
#' are derived from three core terms:
#' 
#' - **n_sb**: Number of sites belonging to bioregion **b** in which species 
#'   **s** is present.  
#' - **n_s**: Total number of sites in which species **s** is present.  
#' - **n_b**: Total number of sites belonging to bioregion **b**.
#' 
#' These **species-to-bioregions** metrics include:  
#' 
#' - [**Specificity**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#specificity-occurrence), 
#'   as described in De Cáceres M & Legendre P (2009).  
#' - [**NSpecificity**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#nspecificity-occurrence) 
#'   is the normalized version of **Specificity** accounting for bioregion size,  
#'   as described in De Cáceres M & Legendre P (2009).  
#' - [**Fidelity**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#fidelity-occurrence)  
#'   as described in De Cáceres M & Legendre P (2009).  
#' - [**IndVal**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#indval-occurrence),  
#'   as described in De Cáceres M & Legendre P (2009).  
#' - [**NIndVal**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#nindval-occurrence)
#'   is the normalized version of **IndVal** accounting for bioregion size,   
#'   as described in De Cáceres M & Legendre P (2009).  
#' - [**Rho**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#rho-occurrence),  
#'   as described in Lenormand *et al.* (2019).
#' 
#' **Species-to-bioregionalization** metrics can also be computed, such as:  
#' 
#' - [**Participation**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#participation-occurrence), 
#'   as described in Denelle *et al.* (2020).  
#' 
#' Abundance-weighted versions of these metrics (`data_type = "abundance"` or 
#' `data_type = "both"`) can also be derived using the following analogous core 
#' terms:  
#' 
#' - **w_sb**: Sum of abundances of species **s** in sites of bioregion **b**.  
#' - **w_s**: Total abundance of species **s**.  
#' - **w_b**: Total abundance of all species present in sites of bioregion **b**.  
#' 
#' These abundance-weighted terms correspond directly to their occurrence-based 
#' counterparts and allow computing abundance versions of the same 
#' metrics. Detailed formulas and examples are provided in the package vignette.
#' 
#' The **second type** of metrics provided by this function is based on the 
#' contribution of each site to a given bioregion (`cluster_on = "site"` or 
#' `cluster_on = "both"`). This is calculated from the bioregion assigned to 
#' each site (as defined in `bioregionalization`), an associated site–species 
#' co-occurrence matrix `comat` (for `"Richness"`, `"Rich_Endemics"` and 
#' `"Prop_Endemics"`, a strict match between site IDs is required) 
#' and/or a site–site similarity metric `similarity` (for `"MeanSim"` and 
#' `"SdSim"`, a strict match between site IDs is required).
#' 
#' These **site-to-bioregions** metrics include:  
#' 
#' - [**Richness**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#richness): 
#'   The number of species present in a site. 
#' - [**Rich_Endemics**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#richendemics): 
#'   The number of species present in a site that are endemic to a bioregion. 
#' - [**Prop_Endemics**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#propendemics): 
#'   The proportion of species present in a site that are endemic to a bioregion. 
#' - [**MeanSim**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#meansim): 
#'   The mean similarity of each site to the sites of every bioregion.  
#' - [**SdSim**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#sdsim): 
#'   The corresponding standard deviation of similarity values.  
#' 
#' **Site-to-bioregionalization** metrics can also be computed, such as:  
#' 
#' - [**Silhouette**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#silhouette), 
#'   as described in Rousseeuw (1987) (similarity-based version).  
#' 
#' The **third type** of metrics provided by this function is based on the 
#' contribution of each site to a given *chorotypes* (`cluster_on = "species"` 
#' or `cluster = "both"`, only when a 
#' cluster has been assigned to species in `bioregionalization`). 
#' This is calculated from the chorotype assigned to each species (as defined in 
#' `bioregionalization`) and an associated site–species co-occurrence matrix 
#' (a strict match between species IDs is required).
#' 
#' In this case, occurrence-based metrics (`data_type = "occurrence"` or 
#' `data_type = "both"`) are derived from three core terms:
#' 
#' - **n_gc**: Number of species belonging to chorotype **c** that are present 
#'   in site **g**.  
#' - **n_g**: Total number of species present in site **g**.  
#' - **n_c**: Total number of species belonging to chorotype **c**.  
#' 
#' As for the first type, all metrics (including their abundance-weighted 
#' versions when `data_type = "abundance"` or `data_type = "both"`) can be 
#' derived from these **site-to-chorotypes** relationships. Further 
#' details, mathematical definitions, and examples are provided in the package 
#' vignette.
#'   
#' @note If `data_type = "auto"`, the choice between occurrence- or abundance-
#' based metrics will be determined automatically from the input data, and a
#' message will explain the choice made.
#'   
#' @references
#' De Cáceres M & Legendre P (2009) Associations between species and groups of 
#' sites: indices and statistical inference. \emph{Ecology} 90, 3566--3574.
#'
#' Denelle P, Violle C & Munoz F (2020) Generalist plants are more competitive 
#' and more functionally similar to each other than specialist plants: insights 
#' from network analyses. \emph{Journal of Biogeography} 47, 1922–-1933.
#' 
#' Lenormand M, Papuga G, Argagnon O, Soubeyrand M, Alleaume S & Luque S (2019)
#' Biogeographical network analysis of plant species distribution in the 
#' Mediterranean region. \emph{Ecology and Evolution} 9, 237--250.
#'
#' Rousseeuw PJ (1987) Silhouettes: A graphical aid to the interpretation and 
#' validation of cluster analysis. \emph{Journal of Computational and Applied 
#' Mathematics} 20, 53--65.
#'
#' @seealso 
#' For more details illustrated with a practical example, 
#' see the vignette: 
#' \url{https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html}.
#' 
#' Associated functions: 
#' [bioregion_metrics] [bioregionalization_metrics]
#'  
#' @author
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) \cr
#' Boris Leroy (\email{leroy.boris@gmail.com}) \cr
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) 
#' 
#' @examples
#' data(fishmat)
#' 
#' fishsim <- similarity(fishmat, metric = "Jaccard")
#' 
#' bioregionalization <- hclu_hierarclust(similarity_to_dissimilarity(fishsim),
#'                                        index = "Jaccard",
#'                                        method = "average",
#'                                        randomize = TRUE,
#'                                        optimal_tree_method = "best",
#'                                        n_clust = c(1,2,3),
#'                                        verbose = FALSE)
#'                                      
#' ind <- site_species_metrics(bioregionalization = bioregionalization,
#'                              bioregion_metrics = "all",
#'                              bioregionalization_metrics = "all",
#'                              data_type = "auto",
#'                              cluster_on = "site",
#'                              comat = fishmat,
#'                              similarity = fishsim,
#'                              include_cluster = TRUE,
#'                              index = 3,
#'                              verbose = TRUE)
#'                                      
#' @importFrom stats aggregate sd
#' 
#' @export
site_species_metrics <- function(bioregionalization,
                                 bioregion_metrics = c("Specificity", 
                                                       "NSpecificity", 
                                                       "Fidelity", 
                                                       "IndVal", 
                                                       "NIndVal", 
                                                       "Rho"),
                                 bioregionalization_metrics = "P",
                                 data_type = "auto",
                                 cluster_on = "site",
                                 comat,
                                 similarity = NULL,
                                 include_cluster = FALSE,
                                 index = names(similarity)[3],
                                 verbose = TRUE){
  
  # Control verbose & include_cluster
  controls(args = verbose, data = NULL, type = "boolean")
  controls(args = include_cluster, data = NULL, type = "boolean")
  
  # Control bioregionalization
  controls(args = NULL, 
           data = bioregionalization, 
           type ="input_bioregionalization")
  
  # Extract node_type
  b_node_type <- bioregionalization$inputs$node_type 
  if(b_node_type != "species"){
    b_site <- bioregionalization$clusters[attr(bioregionalization$clusters, 
                                               "node_type") == "site", 1]
  }
  if(b_node_type != "site"){
    b_species <- bioregionalization$clusters[attr(bioregionalization$clusters, 
                                                  "node_type") == "species", 1]
  }
  
  # Extract data_type
  b_data_type <- bioregionalization$inputs$data_type 
  
  # List of acronyms
  #
  # s = species
  # g = site
  # b = bioregion
  #
  # sb = species-to-bioregions/bioregionalization (sb) and/or 
  #      site-to-chorotypes/chorological (gc) 
  #      according to cluster_on
  # gb = site-to-bioregions/bioregionalization 
  #
  # ind = bioregion metrics set by the user
  # sb_ind = species-to-bioregions (sb) or site-to-chorotypes (gc) metrics
  # gb_ind = site-to-bioregions metrics
  #
  # agind = bioregionalization metrics set by the user
  # sb_agind = species-to-bioregionalization (sb) or 
  #            site-to-chorological classification (gc) metrics
  # gb_agind = site-to-bioregionalization metrics
  #
  # comat_metrics = comat-based metrics
  # comat_gb_metrics = comat-based gb metrics
  # similarity_metrics = similarity-based metrics

  sb_ind <- c("Specificity", "NSpecificity", 
              "Fidelity", 
              "IndVal", "NIndVal", 
              "Rho", 
              "CoreTerms")
  gb_ind <- c("Richness", "Rich_Endemics", "Prop_Endemics", 
              "MeanSim", "SdSim")
  sb_agind <- c("P")
  gb_agind <- c("Silhouette")
  comat_metrics <- c("Specificity", "NSpecificity", "Fidelity", "IndVal", 
                     "NIndVal", "Rho", "CoreTerms", "Richness", "Rich_Endemics", 
                     "Prop_Endemics", "P")
  comat_gb_metrics <- c("Richness", "Rich_Endemics", "Prop_Endemics")
  similarity_metrics <- c("MeanSim", "SdSim", "Silhouette")
  
  # Control bioregion_metrics and set ind
  ind <- NULL
  if(!is.null(bioregion_metrics)){
    controls(args = bioregion_metrics, data = NULL, type = "character_vector")
    if ("all" %in% bioregion_metrics) {
      ind <- c(sb_ind, gb_ind)
    }else{
      ind <- bioregion_metrics
    }
    if (length(intersect(c(sb_ind, gb_ind), ind)) !=
        length(ind)) {
      stop(paste0("One or several bioregion metrics chosen are not", 
                  " available.\n",
                  "Please choose from the following:\n",
                  "Specificity, NSpecificity, Fidelity, IndVal, NIndVal, ",
                  "Rho, CoreTerms, Richness, Rich_Endemics, Prop_Endemics ", 
                  "MeanSim and SdSim"),
           call. = FALSE)
    }
  }
  
  # Control bioregionalization_metrics and set agind
  agind <- NULL
  if(!is.null(bioregionalization_metrics)){
    controls(args = bioregionalization_metrics, 
             data = NULL, 
             type = "character_vector")
    if ("all" %in% bioregionalization_metrics) {
      agind <- c(sb_agind, gb_agind)
    }else{
      agind <- bioregionalization_metrics
    }
    if (length(intersect(c(sb_agind, gb_agind), agind)) !=
        length(agind)) {
      stop(paste0("One or several bioregionalization metrics chosen are not", 
                  " available.\n",
                  "Please choose from the following:\n",
                  "P and Silhouette"),
           call. = FALSE)
    }
  }
  
  # Check if comat and/or similarity are needed and for potential conflicts
  metrics_needed <- c(ind, agind)
  comat_needed <- FALSE
  comat_gb_needed <- FALSE
  if(length(intersect(comat_metrics, metrics_needed))>0){
    comat_needed <- TRUE
  }
  if(length(intersect(comat_gb_metrics, metrics_needed))>0){
    comat_gb_needed <- TRUE
  }
  similarity_needed <- FALSE
  if(length(intersect(similarity_metrics, metrics_needed))>0){
    similarity_needed <- TRUE
  }
  
  if(is.null(comat) & is.null(similarity)){
    stop(paste0("At least comat or similarity should be provided."),
         call. = FALSE)
  }
  
  if(is.null(comat) & comat_needed){
    warning(paste0("Some metrics (", 
                   paste(intersect(metrics_needed, comat_metrics), collapse = ", "),
                   ") will be skipped because no co-occurrence matrix is provided."),
            call. = FALSE)
    ind <- setdiff(ind, comat_metrics)
    if(length(ind)==0){
      ind <- NULL
    }
    agind <- setdiff(agind, comat_metrics)
    if(length(agind)==0){
      agind <- NULL
    }
    comat_needed <- FALSE
    comat_gb_needed <- FALSE
  }
  
  if(is.null(similarity) & similarity_needed){
    warning(paste0("Some metrics (", 
                   paste(intersect(metrics_needed, similarity_metrics), collapse = ", "),
                   ") will be skipped because no similarity is provided."),
            call. = FALSE)
    ind <- setdiff(ind, similarity_metrics)
    if(length(ind)==0){
      ind <- NULL
    }
    agind <- setdiff(agind, similarity_metrics)
    if(length(agind)==0){
      agind <- NULL
    }
    similarity_needed <- FALSE
  }
  
  # Stop if no metrics
  if(is.null(ind) & is.null(agind)){
    stop(paste0("At least one type of metrics with the appropriate inputs ", 
                "should be specified."),
         call. = FALSE)
  }

  # Check type of metrics needed (sb, gb or both)
  type <- "both"
  if((length(intersect(ind, gb_ind)) == 0) & 
     (length(intersect(agind, gb_agind)) == 0)){
    type <- "sb"
  }
  if((length(intersect(ind, sb_ind)) == 0) & 
     (length(intersect(agind, sb_agind)) == 0)){
     type <- "gb"
  }
  
  # Control cluster_on
  controls(args = cluster_on, data = NULL, type = "character")
  if (!(cluster_on %in% c("site", "species", "both"))) {
    stop(paste0("Please choose cluster_on from the following:\n",
                "site, species or both."), 
          call. = FALSE)
  }
  
  # Control conflicts between b_node_type and cluster_on
  if(cluster_on != "species" & b_node_type == "species"){ 
    if(type == "sb"){ # species-to-bioregions (sb)
      stop(paste0("Species-to-bioregions/bioregionalization metrics are not ",
                  "available when no bioregion are assigned to the site in ",
                  "bioregionalization."
      ), call. = FALSE)
    } 
    if(type == "gb"){ # site-to-bioregions (gb)
      stop(paste0("Site-to-bioregions/bioregionalization metrics are not ",
                  "available when no bioregion are assigned to the site in ",
                  "bioregionalization."
      ), call. = FALSE)
    }
    if(type == "both"){ # both (gb & sb)
      stop(paste0("Species/Site-to-bioregions/bioregionalization metrics are ",
                  "not available when no bioregion are assigned to the site in ",
                  "bioregionalization."
      ), call. = FALSE)
    }  
  }
  if(cluster_on != "site" & b_node_type == "site"){
    if(type == "sb" | type == "both"){ # site-to-chorotypes (gc) 
      stop(paste0("Site-to-chorotypes/chorological metrics are not available ",
                  "when no chorotype are assigned to the species 
                  in bioregionalization."), 
           call. = FALSE)
    } 
  }
  
  # Control comat if needed
  if(comat_needed){
    
    controls(args = NULL, data = comat, type = "input_matrix")
    minco <- min(comat)
    if (minco < 0) {
      stop("Negative value(s) detected in comat!", 
           call. = FALSE)
    }
    
    # Attempt to check if comat is based on occurrence or abundance data
    bin <- all(unique(as.vector(comat)) %in% c(0, 1))
    #if(verbose){
    #  if(bin){
    #    message("comat is based on occurrence data.")
    #  }else{
    #    message("comat is based on abundance data.")
    #  }
    #}
    
    comat_site <- rownames(comat)
    comat_species <- colnames(comat)
    
    if(cluster_on != "species"){
      # Check that comat_site are in bioregionalization
      if(length(intersect(b_site, comat_site)) == length(b_site) &
         length(b_site) == length(comat_site)){
        #print("match!")
      }else{
        stop("Site ID in bioregionalization and comat do not match!", 
             call. = FALSE)
      }
      if(sum(b_site == comat_site) != length(b_site)){
        comat <- comat[match(b_site, comat_site),]
        comat_site <- rownames(comat)
      }
    }
    if(type != "gb"){
      if(cluster_on != "site"){
        # Check that comat_species are in bioregionalization
        if(length(intersect(b_species, comat_species)) == length(b_species) &
           length(b_species) == length(comat_species)){
          #print("match!")
        }else{
          stop("Species ID in bioregionalization and comat do not match!", 
               call. = FALSE)
        }
        if(sum(b_species == comat_species) != length(b_species)){
          comat <- comat[,match(b_species, comat_species)]
          comat_species <- colnames(comat)
        }
      }
    }
  }
  
  # Control data_type
  if(comat_needed){
    controls(args = data_type, data = NULL, type = "character")
    if (!(data_type %in% c("auto", "occurrence", "abundance", "both"))) {
      stop(paste0("Please choose data_type from the following:\n",
                  "auto, occurrence, abundance or both."), 
           call. = FALSE)
    }
    
    # auto
    if(data_type == "auto"){
      if(is.na(b_data_type)){
        if(bin){
          data_type <- "occurrence"
          if(verbose){
            message(paste0("No data type detected in bioregionalization and ",
                           "comat is based on occurence data so occurrence-",
                           "based metrics will be computed."))
          }
        }else{
          data_type <- "abundance"
          if(verbose){
            message(paste0("No data type detected in bioregionalization and ",
                           "comat is based on abundance data so abundance-",
                           "based metrics will be computed."))
          }
        }
      }else{
        if(b_data_type == "occurrence"){
          data_type <- "occurrence"
          if(bin){
            if(verbose){
              message(paste0("The bioregionalization is based on occurence data ",
                             "and comat is based on occurence data so occurrence-",
                             "based metrics will be computed."))
            }
          }else{
            if(verbose){
              message(paste0("The bioregionalization is based on occurence data ",
                             "but note that even if comat is based on abundance ", 
                             "data, occurrence-based metrics will be computed. ",
                             "Change data_type to change this beavior."))
            }
          }
          
        }
        if(b_data_type == "abundance"){
          if(bin){
            data_type <- "occurrence"
            if(verbose){
              message(paste0("The bioregionalization is based on abundance data ",
                             "but comat is based on occurence data so ",
                             "occurrence-based metrics will be computed."))
            }
            
          }else{
            data_type <- "abundance"
            if(verbose){
              message(paste0("The bioregionalization is based on abundance data ",
                             "and comat is based on abundance data so ",
                             "abundance-based metrics will be computed."))
            }
          }
        }
        
      }
    }
    #if(data_type == "abundance"  | data_type == "both"){
    #  if(bin){
    #    warning(paste0("comat is based on occurence data so abundance-based ",
    #                   "metrics won't be computed!"))
    #    data_type = "occurrence"
    #  }
    #}
  }
  
  # Control similarity if needed
  if(similarity_needed){
    
    controls(args = NULL, 
             data = similarity, 
             type = "input_conversion_similarity")
    controls(args = index, 
             data = similarity, 
             type = "input_net_index")
    
    similarity <- similarity
    similarity[,3] <- similarity[,index]
    similarity <- similarity[,1:3]
    similarity <- net_to_mat(similarity, 
                             weight = TRUE, 
                             squared = TRUE,
                             symmetrical = TRUE)
    
    sim_site <- rownames(similarity)
    
    if(length(intersect(b_site, sim_site)) == length(b_site) &
       length(b_site) == length(sim_site)){
      #print("match!")
    }else{
      stop("Site ID in bioregionalization and similarity do not match!", 
           call. = FALSE)
    }
    
    if(sum(b_site == sim_site) != length(b_site)){
      similarity <- similarity[match(b_site, sim_site), 
                               match(b_site, sim_site)]
      sim_site <- rownames(similarity)
    }
  }  
  
  # Loop over partitions
  bioregionalization <- bioregionalization$clusters
  nb_partitions <- dim(bioregionalization)[2] - 1
  
  output <- list()
  length(output) <- nb_partitions
  names(output) <- colnames(bioregionalization)[-1]
  
  for(k in 1:nb_partitions){
    
    # sb or gc
    if(type != "gb"){
      # sb
      if(cluster_on != "species"){
        sb <- sbgc(clusters = bioregionalization[attr(bioregionalization, 
                                                      "node_type") == "site",
                                                 (k+1)], 
                   bioregion_metrics = ind,
                   bioregionalization_metrics = agind,
                   comat = comat,
                   type = "sb",
                   data = data_type)
        
        if(!is.null(sb$bioregion1)){
          output[[k]]$species_bioregions <- sb$bioregion1
        }
        if(!is.null(sb$bioregion2)){
          output[[k]]$species_bioregionalization <- sb$bioregion2
        }
      }
      # gc
      if(cluster_on != "site"){
        gc <- sbgc(clusters = bioregionalization[attr(bioregionalization, 
                                                      "node_type") == "species",
                                                 (k+1)], 
                   bioregion_metrics = ind,
                   bioregionalization_metrics = agind,
                   comat = comat,
                   type = "gc",
                   data = data_type)
        
        if(!is.null(gc$bioregion1)){
          output[[k]]$site_chorotypes <- gc$bioregion1
        }
        if(!is.null(gc$bioregion2)){
          output[[k]]$site_chorological <- gc$bioregion2
        }
      }
    }
    
    #gb
    if(type != "sb"){
      
      gb <- gb(clusters = bioregionalization[attr(bioregionalization, 
                                                  "node_type") == "site",
                                             (k+1)],
               bioregion_metrics = ind,
               bioregionalization_metrics = agind,
               comat = comat,
               similarity = similarity,
               #data = data_type,
               include_cluster = include_cluster)
      
      if(!is.null(gb$bioregion1)){
        output[[k]]$site_bioregions <- gb$bioregion1
      }
      if(!is.null(gb$bioregion2)){
        output[[k]]$site_bioregionalization <- gb$bioregion2
      }
    }
    
  }
  
  # Return output
  if(nb_partitions == 1){
    output <- output[[1]]
  }
  
  attr(output, "n_partitions") <- nb_partitions
  attr(output, "cluster_on") <- if(type == "gb") "site" else cluster_on
  attr(output, "clustering_data_type") <- b_data_type
  attr(output, "index_data_type") <- if(type == "gb") NA else data_type
  attr(output, "has_similarity") <- type != "sb"
  attr(output, "has_comat") <- type != "gb"
  
  if(type != "gb") {
    sb_computed <- intersect(ind, sb_ind)
    if(data_type == "occurrence") {
      attr(output, "bioregion_metrics_occ") <- sb_computed
      attr(output, "bioregion_metrics_abd") <- character(0)
      attr(output, "bioregionalization_metrics_occ") <- intersect(agind, sb_agind)
      attr(output, "bioregionalization_metrics_abd") <- character(0)
    } else if(data_type == "abundance") {
      attr(output, "bioregion_metrics_occ") <- character(0)
      attr(output, "bioregion_metrics_abd") <- sb_computed
      attr(output, "bioregionalization_metrics_occ") <- character(0)
      attr(output, "bioregionalization_metrics_abd") <- intersect(agind, sb_agind)
    } else { # both
      attr(output, "bioregion_metrics_occ") <- sb_computed
      attr(output, "bioregion_metrics_abd") <- sb_computed
      attr(output, "bioregionalization_metrics_occ") <- intersect(agind, sb_agind)
      attr(output, "bioregionalization_metrics_abd") <- intersect(agind, sb_agind)
    }
  } else {
    attr(output, "bioregion_metrics_occ") <- character(0)
    attr(output, "bioregion_metrics_abd") <- character(0)
    attr(output, "bioregionalization_metrics_occ") <- character(0)
    attr(output, "bioregionalization_metrics_abd") <- character(0)
  }
  attr(output, "similarity_metrics") <- if(type != "sb") {
    c(intersect(ind, gb_ind), intersect(agind, gb_agind))
  } else character(0)
  
  class(output) <- c("bioregion.site.species.metrics", class(output))
  return(output)
  
}

