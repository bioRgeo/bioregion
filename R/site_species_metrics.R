#' Calculate metrics for sites and species relative to bioregions and chorotypes
#'
#' This function computes metrics that quantify how species and sites relate 
#' to clusters (bioregions or chorotypes). Depending on the type of clustering, 
#' metrics can measure how species are distributed across bioregions (site 
#' clusters), how sites relate to chorotypes (species clusters), or both.
#'
#' @param bioregionalization A `bioregion.clusters` object.
#'
#' @param bioregion_metrics A `character` vector or a single `character` string
#' specifying the metrics to compute for each cluster. Available metrics depend
#' on the type of clustering (see arg `cluster_on`):
#' \itemize{
#'   \item{**When sites are clustered into bioregions** (default case): 
#'   species-level metrics include `"Specificity"`, `"NSpecificity"`, 
#'   `"Fidelity"`, `"IndVal"`, `"NIndVal"`, `"Rho"`, and `"CoreTerms"`. 
#'   Site-level metrics include `"Richness"`, `"Rich_Endemics"`, 
#'   `"Prop_Endemics"`, `"MeanSim"`, and `"SdSim"`.}
#'   \item{**When species are clustered into chorotypes** (e.g., bipartite 
#'   network clustering): site-level metrics include `"Specificity"`, 
#'   `"NSpecificity"`, `"Fidelity"`, `"IndVal"`, `"NIndVal"`, `"Rho"`, 
#'   and `"CoreTerms"`.}
#' }
#' Use `"all"` to compute all available metrics. See Details for metric 
#' descriptions.
#' 
#' @param bioregionalization_metrics A `character` vector or a single `character` 
#' string specifying summary metrics computed across all clusters. These 
#' metrics assess how an entity (species or site) is distributed across the 
#' entire bioregionalization, rather than relative to each individual cluster:
#' \itemize{
#'   \item{`"P"`: Participation coefficient measuring how evenly a species or 
#'   site is distributed across clusters (0 = restricted to one cluster, 
#'   1 = evenly spread).}
#'   \item{`"Silhouette"`: How well a site fits its assigned bioregion compared 
#'   to the nearest alternative bioregion (requires similarity data).}
#' }
#' Use `"all"` to compute all available metrics.
#' 
#' @param data_type A `character` string specifying whether metrics should be 
#' computed based on presence/absence (`"occurrence"`) or abundance values 
#' (`"abundance"`). This affects how Specificity, Fidelity, IndVal, Rho and
#' CoreTerms are calculated:
#' \itemize{
#'   \item{`"auto"` (default): Automatically detected from input data (`bioregionalization` and/or `comat`).}
#'   \item{`"occurrence"`: Metrics based on presence/absence only.}
#'   \item{`"abundance"`: Metrics weighted by abundance values.}
#'   \item{`"both"`: Compute both versions of the metrics.}
#' }
#' 
#' @param cluster_on A `character` string specifying what was clustered in the 
#' bioregionalization, which determines what types of metrics can be computed:
#' \itemize{
#'   \item{`"site"` (default): Sites were clustered into bioregions. Metrics 
#'   describe how each **species** is distributed across bioregions.}
#'   \item{`"species"`: Species were clustered into chorotypes. Metrics describe 
#'   how each **site** relates to chorotypes. Only available when species have 
#'   been assigned to clusters (e.g., bipartite network clustering).}
#'   \item{`"both"`: Compute metrics for both perspectives. Only available 
#'   when both sites and species have cluster assignments.}
#' }
#' 
#' @param comat A site-species `matrix` with sites as rows and species as
#' columns. Values can be occurrence (1/0) or abundance. Required for most
#' metrics.
#' 
#' @param similarity A site-by-site similarity object from [similarity()] or
#' [dissimilarity_to_similarity()]. Required only for similarity-based metrics 
#' (`"MeanSim"`, `"SdSim"`, `"Silhouette"`).
#' 
#' @param include_cluster A `boolean` indicating whether to add an `Assigned` 
#' column in the output, marking `TRUE` for rows where the site belongs to the 
#' bioregion being evaluated. Useful for quickly identifying a site's own 
#' bioregion. Default is `FALSE`.
#' 
#' @param index The name or number of the column to use as similarity. 
#' By default, the third column name of `similarity` is used.
#' 
#' @param verbose A `boolean` indicating whether to 
#' display progress messages. Set to `FALSE` to suppress these messages.
#' 
#' @return 
#' A `list` containing one or more `data.frame` elements, depending on the 
#' selected metrics and clustering type:
#' 
#' **When sites are clustered (`cluster_on = "site"`):**
#' \itemize{
#'   \item{**species_bioregions**: Metrics for each species x bioregion 
#'   combination (e.g., Specificity, IndVal). One row per species x bioregion 
#'   pair.}
#'   \item{**species_bioregionalization**: Summary metrics for each species 
#'   across all bioregions (e.g., Participation coefficient). One row per 
#'   species.}
#'   \item{**site_bioregions**: Metrics for each site x bioregion combination 
#'   (e.g., MeanSim, Richness). One row per site x bioregion pair.}
#'   \item{**site_bioregionalization**: Summary metrics for each site 
#'   (e.g., Silhouette). One row per site.}
#' }
#' 
#' **When species are clustered (`cluster_on = "species"`):**
#' \itemize{
#'   \item{**site_chorotypes**: Metrics for each site x chorotype combination 
#'   (e.g., Specificity, IndVal). One row per site x chorotype pair.}
#'   \item{**site_chorological**: Summary metrics for each site across all 
#'   chorotypes (e.g., Participation coefficient). One row per site.}
#' }
#'
#' Note that if `bioregionalization` contains multiple partitions 
#' (i.e., if `dim(bioregionalization$clusters) > 2`), a nested list will be 
#' returned, with one sublist per partition.
#' 
#' @details
#' This function computes metrics that characterize the relationship between 
#' species, sites, and clusters. The available metrics depend on whether you 
#' clustered sites (into bioregions) or species (into chorotypes).
#' 
#' 
#' ## --- 1. Understanding the two perspectives ---
#' 
#' - **Bioregions** are clusters of sites with similar species composition.
#' - **Chorotypes** are clusters of species with similar distributions.
#' 
#' 
#' In general, the package is designed to cluster sites into bioregions. 
#' However, it is possible to group species into clusters. 
#' We call these species clusters 'chorotypes',
#' following conceptual definitions in the biogeographical literature, to
#' avoid any confusion in the calculation of metrics. 
#' 
#' In some cases, such as bipartite network clustering, both species and sites
#' receive the same clusters. We maintain the name distinction in the 
#' calculation of metrics - but remember that in this case 
#' BIOREGION IDs = CHOROTYPE IDs.
#' The `cluster_on` argument determines 
#' which perspective to use.
#' 
#' 
#' ## --- 2. Metrics when sites are clustered (`cluster_on = "site"` or `cluster_on = "both"`) ---
#' 
#' **Species-per-bioregion metrics** quantify how each species is distributed 
#' across bioregions. 
#' 
#' These metrics are derived from three core terms ([see the online vignette for a visual
#' diagram](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#metric-components)): 
#' 
#' - **n_sb**: Number of sites in bioregion **b** where species 
#' **s** is present
#' - **n_s**: Total number of sites in which species **s** is present.
#' - **n_b**: Total number of sites in bioregion **b**.
#' 
#' Abundance version of these core terms can also be calculated when 
#' `data_type = "abundance"` (or `data_type = "auto"` and 
#' `bioregionalization was based on abundance`):
#' 
#' - **w_sb**: Sum of abundances of species **s** in sites of bioregion **b**. 
#' - **w_s**: Total abundance of species **s**.  
#' - **w_b**: Total abundance of all species present in sites of bioregion **b**.
#' 
#' The species-per-bioregion metrics are (click on metric names to access formulas):
#' - [**Specificity**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#specificity-occurrence): 
#'   Fraction of a species' occurrences found in a given bioregion
#'   (De Cáceres & Legendre 2009). A value of 
#'   1 means the species occurs only in that bioregion.
#' - [**NSpecificity**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#nspecificity-occurrence): 
#'   Normalized specificity that accounts for differences in bioregion size
#'   (De Cáceres & Legendre 2009).
#' - [**Fidelity**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#fidelity-occurrence): 
#'   Fraction of sites in a bioregion where the species occurs
#'   (De Cáceres & Legendre 2009). A value of 1 
#'   means the species is present in all sites of that bioregion.
#' - [**IndVal**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#indval-occurrence): 
#'   Indicator Value = Specificity × Fidelity
#'   (De Cáceres & Legendre 2009). High values identify species 
#'   that are both restricted to and frequent within a bioregion.
#' - [**NIndVal**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#nindval-occurrence): 
#'   Normalized IndVal accounting for bioregion size
#'   (De Cáceres & Legendre 2009).
#' - [**Rho**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#rho-occurrence): 
#'   Standardized contribution index comparing observed vs. expected 
#'   co-occurrence under random association (Lenormand 2019).
#' - **CoreTerms**: Raw counts (n, n_b, n_s, n_sb) for custom calculations.
#' 
#' These metrics can be found in the output slot `species_bioregions`.
#' 
#' **Site-per-bioregion metrics** characterize sites relative to bioregions:
#' 
#' - [**Richness**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#diversity-endemicity-site-metrics): 
#'   Number of species in the site.
#' - [**Rich_Endemics**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#diversity-endemicity-site-metrics): 
#'   Number of species in the site that are endemic to one bioregion.
#' - [**Prop_Endemics**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#diversity-endemicity-site-metrics): 
#'   Proportion of endemic species in the site.
#' - [**MeanSim**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#meansim): 
#'   Mean similarity of a site to all sites in each bioregion.
#' - [**SdSim**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#sdsim): 
#'   Standard deviation of similarity values.
#' 
#' These metrics can be found in the output slot `site_bioregions`.
#' 
#' **Summary metrics across the whole bioregionalization:**
#' 
#' These metrics summarize how an entity (species or site) is distributed 
#' across all clusters, rather than in relation to each individual cluster.
#' 
#' *Species-level summary metric:*
#' - [**P**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#p-occurrence-1) 
#'   (Participation): Evenness of species distribution across bioregions
#'   (Denelle et al. 2020). Found in output slot `species_bioregionalization`.
#'
#' *Site-level summary metric:*
#' - [**Silhouette**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#silhouette): 
#'   How well a site fits its assigned bioregion vs. the nearest alternative 
#'   (Rousseeuw 1987). Found in output slot `site_bioregionalization`.
#' 
#' 
#' ## --- 3. Metrics when species are clustered (`cluster_on = "species"` or `cluster_on = "both"`) ---
#' 
#' **Site-per-chorotype metrics** quantify how each site relates to species 
#' clusters (chorotypes).
#' 
#' The same metrics as above (Specificity, Fidelity, 
#' IndVal, etc.) can be computed, but their interpretation is inverted. These
#' metrics are based on the following core terms:
#' 
#' - **n_gc**: Number of species belonging to chorotype **c** that are present 
#'   in site **g**.  
#' - **n_g**: Total number of species present in site **g**.  
#' - **n_c**: Total number of species belonging to chorotype **c**.
#' 
#' Abundance version of these core terms can also be calculated when 
#' `data_type = "abundance"` (or `data_type = "auto"` and 
#' `bioregionalization was based on abundance`).
#' 
#' Their interpretation changes, for example:
#'   
#' - **Specificity**: Fraction of a site's species belonging to a chorotype.
#' - **Fidelity**: Fraction of a chorotype's species present in the site.
#' - **IndVal**: Indicator value for site-chorotype associations.
#' - **P**: Evenness of sites across chorotypes
#' 
#'   
#' @note If `data_type = "auto"`, the choice between occurrence- or abundance-
#' based metrics will be determined automatically from the input data, and a
#' message will explain the choice made.
#'  
#' Strict matching between entity IDs (site and species IDs) in
#' `bioregionalization` and in `comat` / `similarity` is required.  
#'   
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
  # sb = species-per-bioregion/bioregionalization (sb) and/or 
  #      site-per-chorotypes/chorological (gc) 
  #      according to cluster_on
  # gb = site-per-bioregion/bioregionalization 
  #
  # ind = bioregion metrics set by the user
  # sb_ind = species-per-bioregion (sb) or site-per-chorotypes (gc) metrics
  # gb_ind = site-per-bioregion metrics
  #
  # agind = bioregionalization metrics set by the user
  # sb_agind = species-in-bioregionalization (sb) or 
  #            site-per-chorological classification (gc) metrics
  # gb_agind = site-in-bioregionalization metrics
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
  
  # # Control bioregion_metrics and set ind
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
    if(type == "sb"){ # species-per-bioregion (sb)
      stop(paste0("Species-per-bioregion or species-in-bioregionalization metrics are not ",
                  "available when no bioregion are assigned to the site in ",
                  "bioregionalization."
      ), call. = FALSE)
    } 
    if(type == "gb"){ # site-per-bioregion (gb)
      stop(paste0("Site-per-bioregion or site-in-bioregionalization metrics are not ",
                  "available when no bioregion are assigned to the site in ",
                  "bioregionalization."
      ), call. = FALSE)
    }
    if(type == "both"){ # both (gb & sb)
      stop(paste0("Species/Site-per-bioregion or -in-bioregionalization metrics are ",
                  "not available when no bioregion are assigned to the site in ",
                  "bioregionalization."
      ), call. = FALSE)
    }  
  }
  if(cluster_on != "site" & b_node_type == "site"){
    if(type == "sb" | type == "both"){ # site-per-chorotypes (gc) 
      stop(paste0("Site-per-chorotypes or species-in-chorological metrics are not available ",
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
        
        if(!is.null(sb$bioregion)){
          output[[k]]$species_bioregions <- sb$bioregion
        }
        if(!is.null(sb$bioregionalization)){
          output[[k]]$species_bioregionalization <- sb$bioregionalization
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
        
        if(!is.null(gc$bioregion)){
          output[[k]]$site_chorotypes <- gc$bioregion
        }
        if(!is.null(gc$bioregionalization)){
          output[[k]]$site_chorological <- gc$bioregionalization
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
      
      if(!is.null(gb$bioregion)){
        output[[k]]$site_bioregions <- gb$bioregion
      }
      if(!is.null(gb$bioregionalization)){
        output[[k]]$site_bioregionalization <- gb$bioregionalization
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

