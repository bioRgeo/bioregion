#' Create a map of bioregions
#'
#' This plot function can be used to visualize bioregions based on a 
#' `bioregion.clusters` object combined with a geometry (`sf` objects). 
#'
#' @param clusters An object of class `bioregion.clusters` or a `data.frame`. 
#' If a `data.frame` is used, the first column should represent the sites' 
#' ID, and the subsequent column(s) should represent the clusters.
#'  
#' @param geometry A spatial object that can be handled by the `sf` package. 
#' The first attribute should correspond to the sites' ID (see Details).
#' 
#' @param bioregionalization An `integer`, `character`, or `NULL` specifying 
#' which bioregionalization(s) to plot. If `NULL` (default), all 
#' bioregionalizations are plotted. If an `integer` or vector of `integers`, 
#' bioregionalization(s) are selected by column number(s) in the `clusters` 
#' data.frame (starting from 1 after the ID column). If a `character` or 
#' vector of `characters`, bioregionalization(s) are selected by name(s) 
#' matching column names in `clusters`.
#' 
#' @param write_clusters A `boolean` indicating if the `clusters` should be 
#' added to the `geometry`.
#' 
#' @param plot A `boolean` indicating if the plot should be drawn.
#' 
#' @param ... Further arguments to be passed to `sf::plot()`.
#' 
#' @return One or several maps of bioregions if `plot = TRUE` and the 
#' geometry with additional clusters' attributes if `write_clusters = TRUE`.
#' 
#' @details
#' The `clusters` and `geometry` site IDs should correspond. They should 
#' have the same type (i.e., `character` if `clusters` is a 
#' `bioregion.clusters` object) and the sites of `clusters` should be 
#' included in the sites of `geometry`.
#' 
#' **Bipartite networks**: If the `clusters` object is from a bipartite network
#' (containing both sites and species), only site nodes will be mapped. The 
#' function automatically filters to site nodes using the `node_type` attribute.
#' 
#' **Colors**: If the `clusters` object contains colors (added via 
#' `bioregion_colors()`), these colors will be automatically used for plotting. 
#' Otherwise, the default `sf` color scheme will be applied.
#' 
#' @author
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) \cr
#' Boris Leroy (\email{leroy.boris@gmail.com}) \cr
#' Pierre Denelle (\email{pierre.denelle@gmail.com})
#'
#' @examples
#' data(fishmat)
#' data(fishdf) # (data.frame version of fishmat)
#' data(fishsf)
#' 
#' net <- similarity(fishmat, metric = "Simpson")
#' clu <- netclu_greedy(net)
#' map <- map_bioregions(clu, fishsf, write_clusters = TRUE, plot = FALSE)
#' 
#' # With colors
#' clu_colored <- bioregion_colors(clu)
#' map_bioregions(clu_colored, fishsf, plot = TRUE)
#' 
#' # With bipartite network (sites and species)
#' clu_bip <- netclu_greedy(fishdf, bipartite = TRUE)
#' clu_bip_colored <- bioregion_colors(clu_bip)
#' map_bioregions(clu_bip_colored, fishsf, plot = TRUE)
#' 
#' # With multiple bioregionalizations, plot only specific ones
#' dissim <- dissimilarity(fishmat, metric = "Simpson")
#' clu_multi <- hclu_hierarclust(dissim, n_clust = 2:5)
#' map_bioregions(clu_multi, fishsf, bioregionalization = c(1, 3), plot = TRUE)  # By index
#' map_bioregions(clu_multi, fishsf, bioregionalization = c("K_2", "K_4"), plot = TRUE)  # By name
#' 
#' @importFrom sf st_geometry
#' 
#' @export

map_bioregions <- function(clusters, 
                           geometry,
                           bioregionalization = NULL,
                           write_clusters = FALSE,
                           plot = TRUE, 
                           ...) {

  controls(args = write_clusters, data = NULL, type = "boolean")
  controls(args = plot, data = NULL, type = "boolean")
  
  # Control clusters 
  has_colors <- FALSE
  color_list <- NULL
  
  if (inherits(clusters, "bioregion.clusters")) {
    clu <- TRUE
    df <- clusters$clusters
    
    # Check if colors are available
    if (!is.null(clusters$colors) && is.list(clusters$colors)) {
      has_colors <- TRUE
      # Extract colors for each partition
      # colors is a list with one data.frame per partition
      color_list <- lapply(clusters$colors, function(color_df) {
        # Create a named vector: cluster ID -> color
        setNames(color_df$color, color_df$cluster)
      })
      names(color_list) <- names(clusters$colors)
    }
  }else{  
    # data.frame
    if (!is.data.frame(clusters)) {
      stop(
        "If not a bioregion.clusters's object, clusters must be a data.frame.",
        call. = FALSE)
    }
    # at least two columns
    if (dim(clusters)[2] < 2) {
      stop("clusters must be a data.frame with at least two columns.",
           call. = FALSE)
    }
    # no duplictaed ID
    if (sum(duplicated(clusters[,1])) > 0) {
      message("Duplicated site ID detected!")
    }
    # no NAs
    nbna <- sum(is.na(clusters))
    if (nbna > 0) {
      stop("NA(s) detected in the data.frame!", 
           call. = FALSE)
    }
    clu <- FALSE
    df <- clusters
  }
  
  # Control geometry
  if(class(geometry)[1] != "sf"){
    stop("It seems that the geometry used is not an sf object.",
         call. = FALSE)
  }
  
  # Handle bipartite networks - filter to sites only
  node_type <- attr(df, "node_type")
  if (!is.null(node_type)) {
    # This is a bipartite network with both sites and species
    # Keep only sites for mapping
    df <- df[node_type == "site", , drop = FALSE]
  }
  
  # Handle bioregionalization selection
  partition_names <- colnames(df)[-1]  # Exclude ID column
  n_partitions <- length(partition_names)
  
  if (!is.null(bioregionalization)) {
    # Validate bioregionalization parameter
    if (is.numeric(bioregionalization)) {
      # Integer indices (1-based, after ID column)
      if (any(!(bioregionalization %% 1 == 0))) {
        stop("bioregionalization must be an integer or a vector of integers.", 
             call. = FALSE)
      }
      if (any(bioregionalization < 1) || any(bioregionalization > n_partitions)) {
        stop(paste0("bioregionalization indices must be between 1 and ", n_partitions, 
                    " (number of available bioregionalizations)."), 
             call. = FALSE)
      }
      # Select bioregionalizations by index
      selected_indices <- bioregionalization
    } else if (is.character(bioregionalization)) {
      # Bioregionalization names
      controls(args = bioregionalization, data = NULL, type = "character_vector")
      
      # Check that all names exist
      invalid_names <- setdiff(bioregionalization, partition_names)
      if (length(invalid_names) > 0) {
        stop(paste0("bioregionalization name(s) not found in clusters: ",
                    paste(invalid_names, collapse = ", "), 
                    "\nAvailable bioregionalizations: ",
                    paste(partition_names, collapse = ", ")), 
             call. = FALSE)
      }
      # Select bioregionalizations by name
      selected_indices <- match(bioregionalization, partition_names)
    } else {
      stop("bioregionalization must be NULL, an integer (or vector of integers), or a character (or vector of characters).",
           call. = FALSE)
    }
    
    # Filter df to keep only ID and selected bioregionalization columns
    df <- df[, c(1, selected_indices + 1), drop = FALSE]
    
    # Update partition_names and color_list if colors are present
    partition_names <- colnames(df)[-1]
    if (has_colors && !is.null(color_list)) {
      # Filter color_list to keep only selected partitions
      color_list <- color_list[names(color_list) %in% partition_names]
    }
  }
  
  # Control that cluster IDs are included in geometry IDs
  idc <- df[,1]
  idg <- geometry[, 1, drop = TRUE]
  
  # Check for missing sites
  missing_sites <- setdiff(idc, idg)
  if(length(missing_sites) > 0){
    stop(paste0("Some cluster sites are not found in the geometry:\n",
                "  Missing sites: ", paste(head(missing_sites, 10), collapse = ", "),
                if(length(missing_sites) > 10) paste0(" ... (", length(missing_sites) - 10, " more)") else "",
                "\n  Please ensure that all sites in 'clusters' have corresponding entries in 'geometry'."),
         call. = FALSE)
  }
  
  # Control parameters
  controls(args = write_clusters, data = NULL, type = "boolean")
  controls(args = plot, data = NULL, type = "boolean")
  
  # Prepare geometry
  sp <- geometry[match(idc, idg), ]
  nbsp <- dim(sp)[2]
  nbdf <- dim(df)[2]
  sp <- cbind(sp, df[, -1])
  colnames(sp)[nbsp:(nbsp+nbdf-2)] <- colnames(df)[-1]
  
  # Plot
  if(plot){
    
    geomsp <- sf::st_geometry(sp)
    plotsp <- sp[, -(1:(nbsp-1))]
    nbplotsp <- dim(plotsp)[2]-1
    
    # Get partition names (column names of cluster assignments)
    partition_names <- colnames(df)[-1]
    
    if(nbplotsp == 1){ 
      # Single partition plot
      partition_name <- partition_names[1]
      
      if(has_colors && partition_name %in% names(color_list)) {
        # Use colors from bioregion_colors
        cluster_values <- plotsp[[partition_name]]
        colors_map <- color_list[[partition_name]]
        
        # sf::plot expects colors in order of sorted unique values
        unique_vals <- sort(unique(cluster_values))
        pal <- colors_map[as.character(unique_vals)]
        
        plot(plotsp, pal = pal, ...)
      } else {
        # Use default coloring
        plot(plotsp, ...)
      }
    }else{ 
      # Multiple partitions
      # Note: sf::plot() does not support different palettes for each attribute
      # when plotting multiple attributes at once. We need to plot each partition
      # individually to preserve custom colors.
      
      if(has_colors) {
        # Plot each partition individually to preserve custom colors
        # Set up grid layout based on number of partitions
        if(nbplotsp <= 4) {
          # Single row layout for 2-4 partitions
          ncols <- nbplotsp
          nrows <- 1
        } else {
          # Multiple rows: 2 columns per row
          ncols <- 2
          nrows <- ceiling(nbplotsp / 2)
        }
        
        # Set up the plotting layout
        old_par <- par(mfrow = c(nrows, ncols))
        on.exit(par(old_par), add = TRUE)
        
        # Plot each partition with its colors
        for(i in seq_len(nbplotsp)) {
          partition_name <- partition_names[i]
          
          if(partition_name %in% names(color_list)) {
            # Get colors for this partition
            cluster_values <- plotsp[[partition_name]]
            colors_map <- color_list[[partition_name]]
            # sf::plot expects colors in order of sorted unique values
            unique_vals <- sort(unique(cluster_values))
            pal <- colors_map[as.character(unique_vals)]
            
            # Plot with custom colors
            # Per sf documentation: when using par(mfrow), must set key.pos=NULL and reset=FALSE
            plot(plotsp[i], pal = pal, key.pos = NULL, reset = FALSE, ...)
          } else {
            # No colors for this partition
            plot(plotsp[i], key.pos = NULL, reset = FALSE, ...)
          }
        }
      } else {
        # No colors - use standard multi-panel layout
        mod4q <- floor(nbplotsp/4)
        mod4r <- nbplotsp-mod4q*4
        
        if(mod4q == 0){
          # Plot all partitions at once
          plot(plotsp, ...)
        }else{
          # Multiple panels of 4
          for(k in 1:mod4q){
            grDevices::dev.new()
            start_idx <- (k-1)*4+1
            end_idx <- (k-1)*4+4
            plot(plotsp[start_idx:end_idx], ...)
          }
          if(mod4r>0){
            grDevices::dev.new()
            start_idx <- (mod4q)*4+1
            end_idx <- (mod4q)*4+mod4r
            plot(plotsp[start_idx:end_idx], key.pos=NULL, ...)
          }
        }
      }
    }
  }
  
  # Write
  if(write_clusters){
    return(sp)
  }
}
