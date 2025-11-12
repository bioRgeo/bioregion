#' Calculate contribution metrics of sites and species to each bioregion
#'
#' This function computes metrics to assess the contribution of a given
#' species or site to each bioregion.
#'
#' @param bioregionalization A `bioregion.clusters` object.
#'
#' @param bioregion_indices A `character` vector or a single `character` string
#' specifying the indices to compute. Several indices belonging to different 
#' categories are available: species to bioregions (`"Specificity"`, 
#' `"NSpecificity"`, `"Fidelity"`, `"NFidelity"`, `"IndVal"`, `"NIndVal"`, 
#' `"Rho"` and `"CoreTerms"`), site to bioregions (`"MeanSim"` and 
#' `"SdSim"`), or site to species clusters when a cluster/bioregion has been 
#' assigned to the species  in `bioregionalization` (see Details). If `"all"` 
#' is specified (default), all indices will be calculated.
#' 
#' @param bioregionalization_indices A `character` vector or a single `character` 
#' string specifying the bioregionalization indices to compute. Some aggregated 
#' indices (such as the participation coefficient or Silhouette index) can be 
#' derived from the `bioregion_indices` (see Details). If `"all"` is specified 
#' (default), all bioregionalization indices will be calculated.
#' 
#' @param data_type A `character` string specifying which type of data should 
#' be considered to compute the related indices (`"A"`, `"B"`,
#' `"IndVal"` and `"Rho"`): occurrences or abundances. By default (`"auto"`),
#' the type of data is inferred from `bioregionalization` and/or the provided
#' co-occurrence matrix (argument `comat`). Other available options are
#' `"occurence"`, `"abundance"` or `"both"` (see Details).
#' 
#' @param node_type A `character` string specifying whether the related
#'  indices (`"Specificity"`, `"NSpecificity"`, `"Fidelity"`, `"NFidelity"`, 
#'  `"IndVal"`, `"NIndVal"`, `"Rho"` and `"CoreTerms"`) 
#'  should be computed as contributions from species to bioregions (`"site"` by
#'  default), from sites to species clusters (`"species"`), or for `"both"`.
#'  The latter type of contribution is only available if a cluster has been
#'  assigned to the species in `bioregionalization` (see Details).
#' 
#' @param comat A co-occurrence `matrix` with sites as rows and species as
#' columns.
#' 
#' @param similarity The output object from [similarity()] or
#' [dissimilarity_to_similarity()].
#' 
#' @param index The name or number of the column to use as similarity. 
#' By default, the third column name of `similarity` is used.
#' 
#' @param verbose A `boolean` indicating whether to 
#' display progress messages. Set to `FALSE` to suppress these messages.
#' 
#' @return 
#' A `list` containing between one and six elements (listed below), depending on 
#' the selected indices (`bioregion_indices` and/or `bioregionalization_indices`) 
#' and the type of nodes (`site` and/or `species`).
#' \itemize{
#'   \item{**species_bioregions**: A `data.frame` containing the 
#'   species-to-bioregions indice(s) based on `comat`.}
#'   \item{**species_bioregionalization**: A `data.frame` containing the  
#'   species-to-bioregionalization indice(s) based on `comat`.}
#'   \item{**site_clusters**: A `data.frame` containing the  
#'   site-to-species clusters indice(s) based on `comat`.}
#'   \item{**site_clustering**: A `data.frame` containing the 
#'   site-to-species clustering indice(s) based on `comat`.}
#'   \item{**site_bioregions**: A `data.frame` containing the 
#'   site-to-bioregions indice(s) based on `similarity`.}
#'   \item{**site_bioregionalization**: A `data.frame` containing the 
#'   site-to-bioregionalization indice(s) based on `similarity`.}
#' }
#'
#' Note that if `bioregionalization` contains more than one partition 
#' (i.e., if `dim(bioregionalization$clusters) > 2`), a list of lists will be 
#' returned, with one sublist per partition.
#' 
#' @details
#' The **first type** of indices provided by this function is based on the 
#' contribution of each species to a given bioregion (`node_type = "site"` or 
#' `node_type = "both"`). This is calculated from the bioregion assigned to 
#' each site (as defined in `bioregionalization`) and 
#' an associated site–species co-occurrence matrix `comat` (a strict match 
#' between site IDs is required). 
#' 
#' Occurrence-based indices (`data_type = "occurrence"` or `data_type = "both"`) 
#' are derived from three core terms:
#' 
#' - **n_sb**: Number of sites belonging to bioregion **b** in which species 
#'   **s** is present.  
#' - **n_s**: Total number of sites in which species **s** is present.  
#' - **n_b**: Total number of sites belonging to bioregion **b**.
#' 
#' These **species-to-bioregion** indices include:  
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
#' **Species-to-bioregionalization** indices can also be computed, such as:  
#' 
#' - [**Participation**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#participation-occurrence), 
#'   as described in Denelle *et al.* (2020).  
#' 
#' Abundance-weighted versions of these indices (`data_type = "abundance"` or 
#' `data_type = "both"`) can also be derived using the following analogous core 
#' terms:  
#' 
#' - **w_sb**: Sum of abundances of species **s** in sites of bioregion **b**.  
#' - **w_s**: Total abundance of species **s**.  
#' - **w_b**: Total abundance of all species present in sites of bioregion **b**.  
#' 
#' These abundance-weighted terms correspond directly to their occurrence-based 
#' counterparts and allow computing abundance versions of the same 
#' indices. Detailed formulas and examples are provided in the package vignette.
#' 
#' The **second type** of indices provided by this function is based on the 
#' contribution of each site to a given *species cluster* (`node_type = "species"` 
#' or `node_type = "both"`, only when a 
#' cluster or bioregion has been assigned to species in `bioregionalization`). 
#' This is calculated from the cluster assigned to each species (as defined in 
#' `bioregionalization`) and an associated site–species co-occurrence matrix 
#' (a strict match between species IDs is required).
#' 
#' In this case, occurrence-based indices (`data_type = "occurrence"` or 
#' `data_type = "both"`) are derived from three core terms:
#' 
#' - **n_gc**: Number of species belonging to cluster **c** that are present 
#'   in site **g**.  
#' - **n_g**: Total number of species present in site **g**.  
#' - **n_c**: Total number of species belonging to cluster **c**.  
#' 
#' As for the first type, all indices (including their abundance-weighted 
#' versions when `data_type = "abundance"` or `data_type = "both"`) can be 
#' derived from these **site-to-clusters** relationships. Further 
#' details, mathematical definitions, and examples are provided in the package 
#' vignette.
#' 
#' The **third type** of indices provided by this function, not included by 
#' default, is based on the contribution of each site to a given bioregion. This 
#' is calculated from the bioregion assigned to each site (as defined in 
#' `bioregionalization`) and a site–site similarity metric (`similarity`) (a 
#' strict match between site IDs is required).
#' 
#' These **site-to-bioregion** indices include:  
#' 
#' - [**MeanSim**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#meansim): 
#'   The mean similarity of each site to the sites of every bioregion.  
#' - [**SdSim**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#sdsim): 
#'   The corresponding standard deviation of similarity values.  
#' 
#' **Site-to-bioregionalization** indices can also be computed, such as:  
#' 
#' - [**Silhouette**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#silhouette), 
#'   as described in Rousseeuw (1987) (similarity-based version).  
#'   
#' @note If `data_type = "auto"`, the choice between occurrence- or abundance-
#' based indices will be determined automatically from the input data, and a
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
#'                              bioregion_indices = "all",
#'                              bioregionalization_indices = "all",
#'                              data_type = "auto",
#'                              node_type = "site",
#'                              comat = fishmat,
#'                              similarity = fishsim,
#'                              index = 3,
#'                              verbose = TRUE)
#'                                      
#' @importFrom stats aggregate sd
#' 
#' @export
site_species_metrics <- function(bioregionalization,
                                 bioregion_indices = c("Specificity", 
                                                       "NSpecificity", 
                                                       "Fidelity", 
                                                       "IndVal", 
                                                       "NIndVal", 
                                                       "Rho", 
                                                       "CoreTerms"),
                                 bioregionalization_indices = "P",
                                 data_type = "auto",
                                 node_type = "site",
                                 comat,
                                 similarity = NULL,
                                 index = names(similarity)[3],
                                 verbose = TRUE){
  
  # Control verbose
  controls(args = verbose, data = NULL, type = "boolean")
  
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
  
  # Control indices
  
  # list of acronyms used:
  # s = species
  # g = site
  # b = bioregion
  # sb = species-to-bioregion and/or site-to-species cluster (gc) according to node_type
  # gb = site-to-bioregion 
  # sb_ind = indices based on species-to-bioregion relationships or 
  #          site-to-species cluster relationships (gc)
  # gb_ind = indices based on site-to-bioregion relationships
  # ind = bioregion indices
  # sb_agind = indices based on species-to-bioregionalization relationships or 
  #            site-to-species clustering relationships (gc)
  # gb_agind = indices based on site-to-bioregionalization relationships
  # agind = bioregionalization indices
  
  sb_ind <- c("Specificity", "NSpecificity", 
              "Fidelity", 
              "IndVal", "NIndVal", 
              "Rho", 
              "CoreTerms")
  gb_ind <- c("MeanSim", "SdSim")
  
  sb_agind <- c("P")
  gb_agind <- c("Silhouette")
  
  indnull <- FALSE 
  if(!is.null(bioregion_indices)){
    controls(args = bioregion_indices, data = NULL, type = "character_vector")
    if ("all" %in% bioregion_indices) {
      ind <- c(sb_ind, gb_ind)
    }else{
      ind <- bioregion_indices
    }
    if (length(intersect(c(sb_ind, gb_ind), ind)) !=
        length(ind)) {
      stop(paste0("One or several bioregion indices chosen are not", 
                  " available.\n",
                  "Please choose from the following:\n",
                  "Specificity, NSpecificity, Fidelity, IndVal, NIndVal, ",
                  "Rho and CoreTerms"),
           call. = FALSE)
    }
    # Check type of bioregion_indices (sb, gb or both)
    typeind <- "both"
    if(length(intersect(ind, gb_ind)) == 0){
      typeind <- "sb"
    }
    if(length(intersect(ind, sb_ind)) == 0){
      typeind <- "gb"
    }
  }else{
    indnull <- TRUE
    ind <- NULL
  }
  
  agindnull <- FALSE
  if(!is.null(bioregionalization_indices)){
    controls(args = bioregionalization_indices, 
             data = NULL, 
             type = "character_vector")
    if ("all" %in% bioregionalization_indices) {
      agind <- c(sb_agind, gb_agind)
    }else{
      agind <- bioregionalization_indices
    }
    if (length(intersect(c(sb_agind, gb_agind), agind)) !=
        length(agind)) {
      stop(paste0("One or several bioregionalization indices chosen are not", 
                  " available.\n",
                  "Please choose from the following:\n",
                  "P, Silhouette and Ratio"),
           call. = FALSE)
    }
    # Check type of bioregionalization_indices (sb, gb or both)
    typeagind <- "both"
    if(length(intersect(agind, gb_agind)) == 0){
      typeagind <- "sb"
    }
    if(length(intersect(agind, sb_agind)) == 0){
      typeagind <- "gb"
    }
  }else{
    agindnull <- TRUE
    agind <- NULL
  }
  
  if(indnull & agindnull){
    stop(paste0("At least one type of indices should be specified."),
         call. = FALSE)
  }
  
  # Check type of indices needed (sb, gb or both)
  type <- "both"
  if(indnull){
    type <- typeagind
  }
  if(agindnull){
    type <- typeind
  }
  if(!indnull & !agindnull){
    if(typeind == "sb" & typeagind == "sb"){
      type = "sb"
    }
    if(typeind == "gb" & typeagind == "gb"){
      type = "gb"
    }
  }
  
  # Control node_type
  if(type != "gb"){
    controls(args = node_type, data = NULL, type = "character")
    if (!(node_type %in% c("site", "species", "both"))) {
      stop(paste0("Please choose node_type from the following:\n",
                  "site, species or both."), 
           call. = FALSE)
    }
  }
  
  # Control conflicts
  if(type == "gb" | type == "both"){ # site-to-bioregions (gb) included
    if(b_node_type == "species"){
      stop(paste0("Site-to-bioregion metrics are not available ",
                  "when no bioregion are assigned to the site in ",
                  "bioregionalization."
      ), call. = FALSE)
    }
  }
  if(type == "sb" | type == "both"){ # species-to-bioregions (sb) included
    # and/or 
    # site-to-species_clusters (gc) included
    if(b_node_type == "site" & node_type != "site"){
      stop(paste0("Site-to-species_cluster metrics are not available ",
                  "when no clusters/bioregions are assigned to the species 
                  in bioregionalization."), 
           call. = FALSE)
    }
  }
  if(type == "sb"){ # species-to-bioregions (sb) only
    # and/or 
    # site-to-species_clusters (gc) only
    if(b_node_type == "species" & node_type != "species"){
      stop(paste0("Species-to-bioregion metrics are not available ",
                  "when no bioregion are assigned to the site in ",
                  "bioregionalization."), 
           call. = FALSE)
    }
  }
  
  # Control comat only if species-to-bioregion needed (or site-to-species cluster)
  if(type != "gb"){
    if(is.null(comat)){
      stop(paste0("comat is missing!"),
           call. = FALSE)
    }
    
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
    
    if(node_type != "species"){
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
    if(node_type != "site"){
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
  
  # Control data_type
  if(type != "gb"){
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
                           "based indices will be computed."))
          }
        }else{
          data_type <- "abundance"
          if(verbose){
            message(paste0("No data type detected in bioregionalization and ",
                           "comat is based on abundance data so abundance-",
                           "based indices will be computed."))
          }
        }
      }else{
        if(b_data_type == "occurrence"){
          data_type <- "occurrence"
          if(bin){
            if(verbose){
              message(paste0("The bioregionalization is based on occurence data ",
                             "and comat is based on occurence data so occurrence-",
                             "based indices will be computed."))
            }
          }else{
            if(verbose){
              message(paste0("The bioregionalization is based on occurence data ",
                             "but note that even if comat is based on abundance ", 
                             "data, occurrence-based indices will be computed. ",
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
                             "occurrence-based indices will be computed."))
            }
            
          }else{
            data_type <- "abundance"
            if(verbose){
              message(paste0("The bioregionalization is based on abundance data ",
                             "and comat is based on abundance data so ",
                             "abundance-based indices will be computed."))
            }
          }
        }
        
      }
    }
    #if(data_type == "abundance"  | data_type == "both"){
    #  if(bin){
    #    warning(paste0("comat is based on occurence data so abundance-based ",
    #                   "indices won't be computed!"))
    #    data_type = "occurrence"
    #  }
    #}
  }
  
  # Control similarity only if site-to-bioregion (gb) needed
  if(type != "sb"){
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
      if(node_type != "species"){
        sb <- sbgc(clusters = bioregionalization[attr(bioregionalization, 
                                                      "node_type") == "site",
                                                 (k+1)], 
                   comat = comat,
                   bioregion_indices = ind,
                   bioregionalization_indices = agind,
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
      if(node_type != "site"){
        gc <- sbgc(clusters = bioregionalization[attr(bioregionalization, 
                                                      "node_type") == "species",
                                                 (k+1)], 
                   comat = comat,
                   bioregion_indices = ind,
                   bioregionalization_indices = agind,
                   type = "gc",
                   data = data_type)
        
        if(!is.null(gc$bioregion)){
          output[[k]]$site_clusters <- gc$bioregion
        }
        if(!is.null(gc$bioregionalization)){
          output[[k]]$site_clustering <- gc$bioregionalization
        }
      }
    }
    
    #gb
    if(type != "sb"){
      gb <- gb(clusters = bioregionalization[attr(bioregionalization, 
                                                  "node_type") == "site",
                                             (k+1)],  
               similarity = similarity,
               bioregion_indices = ind,
               bioregionalization_indices = agind)
      
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
  
  class(output) <- append(class(output), "bioregion.site_species_metrics")
  return(output)
  
}

