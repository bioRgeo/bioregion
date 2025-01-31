#' Calculate contribution metrics for bioregions
#' 
#' This function calculates the number of sites per bioregion, as well as the
#' number of species these sites have, the number of endemic species, and the 
#' proportion of endemism.
#' 
#' @param bioregionalization A `bioregion.clusters` object.
#' 
#' @param comat A co-occurrence `matrix` with sites as rows and species as 
#' columns. 
#' 
#' @param map A spatial `sf data.frame` with sites and bioregions. It is the 
#' output of the function `map_bioregions`. `NULL` by default.
#'
#' @param col_bioregion An `integer` specifying the column position of the 
#' bioregion.
#' 
#' @return 
#' A `data.frame` with 5 columns, or 6 if spatial coherence is computed.
#'
#' @details 
#' Endemic species are species found only in the sites belonging to one 
#' bioregion.
#'
#' @seealso 
#' For more details illustrated with a practical example, 
#' see the vignette: 
#' \url{https://biorgeo.github.io/bioregion/articles/a5_3_summary_metrics.html}.
#' 
#' Associated functions: 
#' [site_species_metrics] [bioregionalization_metrics]
#'  
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) \cr
#' Boris Leroy (\email{leroy.boris@gmail.com}) \cr
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) 
#' 
#' @examples
#' comat <- matrix(sample(1000, 50), 5, 10)
#' rownames(comat) <- paste0("Site", 1:5)
#' colnames(comat) <- paste0("Species", 1:10)
#'
#' net <- similarity(comat, metric = "Simpson")
#' clust <- netclu_louvain(net)
#' 
#' bioregion_metrics(bioregionalization = clust, 
#'                   comat = comat) 
#' 
#' @export
bioregion_metrics <- function(bioregionalization, 
                              comat,
                              map = NULL, 
                              col_bioregion = NULL){
  
  # 1. Controls ---------------------------------------------------------------
  # input can be of format bioregion.clusters
  if (inherits(bioregionalization, "bioregion.clusters")) {
    if (inherits(bioregionalization$clusters, "data.frame")) {
      has.clusters <- TRUE
      clusters <- bioregionalization$clusters
      
      if(ncol(clusters) > 2) {
        stop(paste0("This function is designed to be applied on a single ",
                    "bioregionalization."), 
             call. = FALSE)
      }
      
    } else {
      if (bioregionalization$name == "hclu_hierarclust") {
        stop(paste0("No clusters have been generated for your hierarchical ",
                    "tree, please extract clusters from the tree before using ",
                    "bioregionalization_metrics().\n",
                    "See ?hclu_hierarclust or ?cut_tree."), 
             call. = FALSE)
      } else {
        stop(paste0("bioregionalization does not have the expected type of ",
                    "'clusters' slot."), 
             call. = FALSE)
      }
    }
  } else {
    stop(paste0("This function is designed to work on bioregion.clusters ",
                "objects and on a site x species matrix."), 
         call. = FALSE)
  }
  
  controls(args = NULL, data = comat, type = "input_matrix")
  
  if(!is.null(map)){
    if(!("sf" %in% class(map))){
      stop(paste0("map must be a 'sf' spatial data.frame with bioregions ",
                  "and sites."), 
           call. = FALSE)
    }
    if(!("data.frame" %in% class(map))){
      stop(paste0("map must be a 'sf' spatial data.frame with bioregions ",
                  "and sites."), 
           call. = FALSE)
    }
    if(ncol(map) < 3){
      stop(paste0("map must have at least 3 columns: sites, bioregions and ",
                  "geometry."), 
           call. = FALSE)
    }
    if(is.null(col_bioregion)){
      stop(paste0("col_bioregion must be defined ", 
                  "it is the column position ",
                  "of the bioregion."), 
           call. = FALSE)
    }
  }
  
  if(!is.null(col_bioregion)){
    if(is.null(map)){
      warning(paste0("col_bioregion is defined but is not considered since ",
                     "map is set to NULL."))
    }else{
      controls(args = col_bioregion, 
               data = NULL, 
               type = "strict_positive_integer")
      
      map_test <- map
      sf::st_geometry(map_test) <- NULL
      if(inherits(map_test[, col_bioregion], "logical")){
        stop("There is no bioregion in the Bioregion column.", 
             call. = FALSE)
      }
      rm(map_test)
    }
  }
  
  bioregion_df <- geometry <- NULL
  
  # 2. Function ---------------------------------------------------------------
  # If it is a bipartite object, we can directly count the species
  if(bioregionalization$inputs$bipartite == TRUE){
    clusters <- clusters[which(attributes(clusters)$node_type == "site"), ]
  }
  
  bioregion_df <- data.frame()
  
  # Loop over bioregions
  for(j in 1:bioregionalization$cluster_info$n_clust){
    focal_j <- unique(clusters[, 2])[j] # bioregion j
    
    # Sites belonging to bioregion j
    sites_j <- clusters[which(clusters[, 2] == focal_j), "ID"]
    
    # Subset site x species matrix with sites belonging to bioregion j
    comat_j <- comat[sites_j, colSums(comat[sites_j, ]) > 0, drop = FALSE]
    
    # Other sites
    comat_not_j <- comat[-which(rownames(comat) %in% sites_j), ]
    comat_not_j <- comat_not_j[, colSums(comat_not_j) > 0, drop = FALSE]
    
    # Number of endemics
    endemic_j <- sum(!(colnames(comat_j) %in% colnames(comat_not_j)))
    
    # Output
    bioregion_df <-
      rbind(bioregion_df,
            data.frame(Bioregion = focal_j,
                       Site_number = length(sites_j),
                       Species_number = ncol(comat_j),
                       Endemics = endemic_j,
                       Percentage_Endemic = 100*endemic_j/ncol(comat_j)))
  }
  
  # Test if all bioregions are in the output
  if(length(unique(bioregion_df$Bioregion)) !=
     bioregionalization$cluster_info$n_clust){
    warning("Not all bioregions are in the output.")
  }
  
  ## 2.2. Spatial coherence ---------------------------------------------------
  if(!is.null(map)){
    # Rename column with bioregion
    colnames(map)[col_bioregion] <- "Bioregion"
    
    # If one bioregion only: spatial coherence equals 100%
    if(dplyr::n_distinct(map$Bioregion) == 1){
      bioregion_df <- data.frame(bioregion_df,
                                 Coherence = 100)
    } else if(dplyr::n_distinct(map$Bioregion) > 1){
      spatial_coherence_df <- data.frame()
      for(j in 1:dplyr::n_distinct(map$Bioregion)){
        bioregion_j <-
          unique(map[which(map$Bioregion ==
                             unique(map$Bioregion)[j]), ]$Bioregion)
        
        # Merging sites belonging to bioregion j that touch each other
        map_j <- map[which(map$Bioregion == bioregion_j), ] 
        map_j <- dplyr::summarise(map_j, geometry = sf::st_union(geometry))
        map_j <- sf::st_cast(map_j, "POLYGON")
        map_j <- dplyr::mutate(map_j, ID = dplyr::row_number())
        
        # Adding area of each touching entity
        map_j$area <- as.numeric(sf::st_area(map_j))
        
        # Spatial coherence of bioregion i
        sp_coherence_j <- 100* max(map_j$area, na.rm = TRUE) /
          sum(map_j$area, na.rm = TRUE)
        
        # Storing spatial coherence
        spatial_coherence_df <- rbind(spatial_coherence_df,
                                      data.frame(Bioregion = bioregion_j,
                                                 Coherence = sp_coherence_j))
      }
      
      bioregion_df <- dplyr::left_join(bioregion_df, spatial_coherence_df,
                                       by = "Bioregion")
    } else{
      stop("There is no bioregion in the Bioregion column.", 
           call. = FALSE)
    }
  }
  
  return(bioregion_df)    
}
