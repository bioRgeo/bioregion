#' Calculate contribution metrics for bioregions
#' 
#' This function calculates the number of sites per bioregion, as well as the
#' the number of species these sites have, the number of endemic species and 
#' the proportion of endemism.
#' 
#' @param cluster_object a `bioregion.clusters` object or a `data.frame` or a 
#' list of `data.frame` containing multiple partitions. At least two partitions
#' are required. If a list of `data.frame` is provided, they should all have
#' the same number of rows (i.e., same items in the clustering for all
#' partitions). 
#' 
#' @param comat a co-occurrence `matrix` with sites as rows and species as
#' columns. 
#' 
#' @param map a spatial `sf data.frame` with sites and bioregions. It is the
#' output of the function `map_bioregions`. `NULL` by default.
#'
#' @param col_bioregion an `integer` specifying the column position of the
#' bioregion.
#'
#' @details 
#' Endemic species are species only found in the sites belonging to one
#' bioregion.
#'
#' @seealso [site_species_metrics]
#' @return 
#' A `data.frame` with 5 columns, or 6 if the spatial coherence is computed.
#'  
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}),
#' Boris Leroy (\email{leroy.boris@gmail.com}) and
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) 
#' 
#' @examples
#' comat_1 <- matrix(sample(0:1000, size = 10*12, replace = TRUE,
#' prob = 1/1:1001), 10, 12)
#' rownames(comat_1) <- paste0("Site", 1:10)
#' colnames(comat_1) <- paste0("Species", 1:12)
#' comat_1 <- cbind(comat_1,
#'                  matrix(0, 10, 8,
#'                         dimnames = list(paste0("Site", 1:10),
#'                                         paste0("Species", 13:20))))
#' 
#' comat_2 <- matrix(sample(0:1000, size = 10*12, replace = TRUE,
#'                          prob = 1/1:1001), 10, 12)
#' rownames(comat_2) <- paste0("Site", 11:20)
#' colnames(comat_2) <- paste0("Species", 9:20)
#' comat_2 <- cbind(matrix(0, 10, 8,
#'                         dimnames = list(paste0("Site", 11:20),
#'                                         paste0("Species", 1:8))),
#'                  comat_2)
#' 
#' comat <- rbind(comat_1, comat_2)
#' 
#' dissim <- dissimilarity(comat, metric = "Simpson")
#' clust1 <- nhclu_kmeans(dissim, n_clust = 3, index = "Simpson")
#' 
#' net <- similarity(comat, metric = "Simpson")
#' com <- netclu_greedy(net)
#' 
#' bioregion_metrics(cluster_object = clust1, comat = comat) 
#' 
#' # Spatial coherence
#' vegedissim <- dissimilarity(vegemat)
#' hclu <- nhclu_kmeans(dissimilarity = vegedissim, n_clust = 4)
#' vegemap <- map_bioregions(hclu, vegesf, write_clusters = TRUE, plot = FALSE)
#' 
#' bioregion_metrics(cluster_object = hclu, comat = vegemat, map = vegemap,
#' col_bioregion = 2) 
#' 
#' @export

bioregion_metrics <- function(cluster_object, comat,
                              map = NULL, col_bioregion = NULL){
  # 1. Controls ---------------------------------------------------------------
  # input can be of format bioregion.clusters
  if (inherits(cluster_object, "bioregion.clusters")) {
    if (inherits(cluster_object$clusters, "data.frame")) {
      has.clusters <- TRUE
      clusters <- cluster_object$clusters
      
      if(ncol(clusters) > 2) {
        stop("This function is designed to be applied on a single partition.",
             "Your cluster_object has multiple partitions (select only one).")
      }
      
    } else {
      if (cluster_object$name == "hierarchical_clustering") {
        stop("No clusters have been generated for your hierarchical tree,
        please extract clusters from the tree before using partition_metrics()
        See ?hclu_hierarclust or ?cut_tree")
      } else {
        stop(
          "cluster_object does not have the expected type of 'clusters' slot")
      }
    }
  } else {
    stop("This function is designed to work on bioregion.clusters objects and
         on a site x species matrix.")
  }
  
  controls(args = NULL, data = comat, type = "input_matrix")
  
  if(!is.null(map)){
    if(!("sf" %in% class(map))){
      stop("map must be a 'sf' spatial data.frame with bioregions and sites.")
    }
    if(!("data.frame" %in% class(map))){
      stop("map must be a 'sf' spatial data.frame with bioregions and sites.")
    }
    if(ncol(map) < 3){
      stop("map must have at least 3 columns: sites, bioregions and geometry.")
    }
    if(is.null(col_bioregion)){
      stop("col_bioregion must be defined, it is the column position of the
           bioregion.")
    }
  }
  
  if(!is.null(col_bioregion)){
    if(is.null(map)){
      warning("col_bioregion is defined but is not considered since map is set
              to NULL.")
    }
    
    controls(args = col_bioregion, data = NULL, type = "positive_integer")
    
    map_test <- map
    sf::st_geometry(map_test) <- NULL
    if(inherits(map_test[, col_bioregion], "logical")){
      stop("There is no bioregion in the Bioregion column.")
    }
    rm(map_test)
  }
  
  bioregion_df <- geometry <- NULL
  
  # 2. Function ---------------------------------------------------------------
  # If it is a bipartite object, we can directly count the species
  if(cluster_object$inputs$bipartite == TRUE){
    clusters <- clusters[which(attributes(clusters)$node_type == "site"), ]
  }
  
  bioregion_df <- data.frame()
  
  # Loop over bioregions
  for(j in 1:cluster_object$cluster_info$n_clust){
    focal_j <- unique(clusters[, 2])[j] # bioregion j
    
    # Sites belonging to bioregion j
    sites_j <- clusters[which(clusters[, 2] == focal_j), "ID"]
    
    # Subset site x species matrix with sites belonging to bioregion j
    comat_j <- comat[sites_j, colSums(comat[sites_j, ]) > 0, drop = FALSE]
    
    # Other sites
    comat_not_j <- comat[-which(rownames(comat) %in% sites_j), ]
    comat_not_j <- comat_not_j[, colSums(comat_not_j) > 0]
    
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
     cluster_object$cluster_info$n_clust){
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
      stop("There is no bioregion in the Bioregion column.")
    }
  }
  
  return(bioregion_df)    
}
