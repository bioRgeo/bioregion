#' Calculate metrics for bioregions
#' 
#' This function calculates the number of sites, the number of species, 
#' the number of endemic species and the proportion of endemism per bioregion. The 
#' spatial coherence can be optionally computed if a spatial object is provided. 
#' 
#' @param bioregionalization A `bioregion.clusters` object.
#' 
#' @param comat A site-species `matrix` with sites as rows and species as
#' columns.
#' 
#' @param map A spatial object that can be handled by `sf` or `terra`. 
#' The first attribute or layer should correspond to the sites' ID 
#' (see Details). Needed only for the spatial coherence (`NULL` by default).
#'
#' @param col_bioregion Deprecated.
#' 
#' @return 
#' A `data.frame` with 5 columns (Bioregion ID and metrics described below) or 
#' 7 if spatial coherence is computed.
#' 
#' - **NbSites**: Number of sites per bioregion
#' - **Richness**: Number of distinct species per bioregion.
#' - **Rich_Endemics**: Number of species found only in the bioregion.
#' - **Prop_Endemics**: Fraction of endemics species.
#' - [**SC_size**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#bioregion-metrics-spatial-coherence): 
#' Spatial coherence based on size, fraction of the number of site contained in the bioregion's largest contiguous patch.
#' - [**SC_area**](https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html#bioregion-metrics-spatial-coherence): 
#' Spatial coherence based on area, fraction of the bioregion area contained in its largest contiguous patch.
#' 
#' Note that if `bioregionalization` contains multiple partitions 
#' (i.e., if `dim(bioregionalization$clusters) > 2`), a list will be 
#' returned.
#'
#' @details
#' `map` should be the output of 
#' `map_bioregions(bioregionalization,   
#'                 geometry,
#'                 write_clusters = TRUE)`
#' 
#' @seealso 
#' For more details illustrated with a practical example, 
#' see the vignette: 
#' \url{https://biorgeo.github.io/bioregion/articles/a5_2_summary_metrics.html}.
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
  
  # Control col_bioregion
  if (!is.null(col_bioregion)) {
    warning("col_bioregion is deprecated.", 
            call. = FALSE)
  }
  
  # Control bioregionalization
  controls(args = NULL, 
           data = bioregionalization, 
           type ="input_bioregionalization")
  
  # Extract sites
  if(bioregionalization$inputs$node_type == "species"){
    stop(paste0("No bioregions are assigned to the sites in ",
                "bioregionalization."),
         call. = FALSE)
  }else{
    clusters <- as.data.frame(bioregionalization$clusters[attr(bioregionalization$clusters, 
                                                               "node_type") == "site",])
    attr(clusters, "node_type") <- NULL
    b_site <- clusters[,1]
  }
  
  # Control comat
  controls(args = NULL, data = comat, type = "input_matrix")
  minco <- min(comat)
  if (minco < 0) {
    stop("Negative value(s) detected in comat!", 
         call. = FALSE)
  }
  comat_site <- rownames(comat)
  
  # Check comat_site & b_site
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
  
  # Control map
  if(!is.null(map)){
    map <- map_bioregions(bioregionalization,
                          map = map,
                          partition_index = NULL,
                          map_as_output = TRUE,
                          plot = FALSE)
  }
  
  # Loop over partitions
  bioregionalization <- bioregionalization$clusters
  nb_partitions <- dim(bioregionalization)[2] - 1
  
  output <- list()
  length(output) <- nb_partitions
  names(output) <- colnames(bioregionalization)[-1]
  
  for(k in 1:nb_partitions){
    
    # partition
    partition <- clusters[,(k+1)]
    
    # comat_bin
    comat_bin <- comat
    comat_bin[comat_bin > 0] <- 1
    
    # nj
    tab <- stats::aggregate(comat_bin[,1], list(partition), length)
    colnames(tab) <- c("Bioregion", "NbSites")
    
    # nij
    temp <- stats::aggregate(comat_bin, list(partition), sum)
    nij <- t(as.matrix(temp[,-1]))
    nij <- (nij > 0)*1
    rownames(nij) <- colnames(temp)[-1]
    colnames(nij) <- temp[,1]

    # Richness
    tab$Richness <- apply(nij, 2, sum)
    
    # Rich_Endemics
    endemics <- (apply(nij, 1, sum) == 1)
    endemics <- nij*endemics
    tab$Rich_Endemics <- apply(endemics, 2, sum)
    
    # Prop_Endemics
    tab$Prop_Endemics <- tab$Rich_Endemics / tab$Richness
    
    # Spatial coherence
    if(!is.null(map)){
      
      SC <- NULL
      SCs <- NULL
      
      # Loop over the bioregions
      for (i in 1:dim(tab)[1]) {
        
        if(tab$NbSites[i] == 1){
          SC <- c(SC, 1)
          SCs <- c(SCs, 1)
        }else{
          
          # Subset of map
          mapki <- map[, (k+1)]
          mapki <- mapki[mapki[, 1, drop = TRUE] == tab$Bioregion[i],]
          
          # Neighbors
          nb <- sf::st_touches(mapki, sparse = TRUE)
          
          # Graph
          g <- igraph::graph_from_adj_list(nb)
          
          # Patches
          patch <- igraph::components(g)$membership
          
          # Area
          areas <- as.numeric(sf::st_area(mapki))
          
          # Patch size or area
          patch_size <- aggregate(areas, list(patch), length)[,2]
          patch_area <- aggregate(areas, list(patch), sum)[,2]
          
          # Spatial coherence
          SC <- c(SC, max(patch_size) / sum(patch_size))
          SCs <- c(SCs, max(patch_area) / sum(patch_area))
          
        }
      }
       
      tab$SC_Size <- SC
      tab$SC_Area <- SCs
    }
    
    # Update output
    output[[k]] <- tab
             
  }
  
  # Return output
  if(nb_partitions == 1){
    output <- output[[1]]
  }
  
  return(output)

}
