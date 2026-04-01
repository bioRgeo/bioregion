#' Create a map of bioregions
#'
#' This plot function can be used to visualize bioregions based on a 
#' `bioregion.clusters` object combined with a spatial object (`sf` or `terra`). 
#'
#' @param bioregionalization A `bioregion.clusters` object.
#'  
#' @param map A spatial object that can be handled by `sf` or `terra`. 
#' The first attribute or layer should correspond to the sites' ID 
#' (see Details).
#' 
#' @param partition_index An `integer`, `character`, or `NULL` specifying 
#' which `bioregionalization`'s partition(s) to plot. By default (`NULL`), all 
#' partitions are plotted. If an `integer` or vector of `integers` is provided, 
#' partition(s) are selected by column number(s) in the `bioregionalization` 
#' data.frame (starting from 1 after the ID column). If a `character` or 
#' vector of `characters`, partition(s) are selected by name(s) 
#' matching column names in `bioregionalization`.
#' 
#' @param map_as_output A `boolean` indicating if the `sf` `data.frame` object 
#' used for the plot should be returned.
#' 
#' @param plot A `boolean` indicating if the plot should be drawn.
#' 
#' @param ... Further arguments to be passed to `sf::plot()`.
#' 
#' @param clusters Deprecated. Use `bioregionalization` instead. The former
#' `bioregionalization` has been replaced by `partition_index`.
#' 
#' @param geometry Deprecated. Use `map` instead.
#' 
#' @param write_clusters Deprecated. Use `map_as_output` instead.
#' 
#' @return One or several maps of bioregions if `plot = TRUE` and the 
#' `sf` `data.frame` object used for the plot if `map_as_output = TRUE`.
#' 
#' @details
#' The site IDs in `bioregionalization` and `map` should correspond. They must 
#' have the same type (i.e., `character` if `bioregionalization` is a 
#' `bioregion.clusters` object), and the sites in `bioregionalization` should be 
#' included among the sites in `map`. If `map` is an `sf` or a `SpatVector` 
#' (`terra`) object, it should contain an attribute table with the IDs in the 
#' first column. If `map` is a `SpatRaster` (`terra`) object, it should contain 
#' the IDs in the first layer.
#' 
#' If the `bioregionalization` object contains both types of nodes (sites and 
#' species), only site  will be mapped. The function automatically filters to 
#' site nodes using the `node_type` attribute.
#' 
#' **Colors**: If the `bioregionalization` object contains colors (added via 
#' `bioregion_colors()`), these colors will be automatically used for plotting. 
#' Otherwise, the default `sf` color scheme will be applied.
#' 
#' @seealso 
#' For more details illustrated with a practical example, 
#' see the vignette: 
#' \url{https://biorgeo.github.io/bioregion/articles/a5_1_visualization.html}.
#' 
#' Associated functions: 
#' [bioregion_colors]
#' 
#' @author
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) \cr
#' Boris Leroy (\email{leroy.boris@gmail.com}) \cr
#' Pierre Denelle (\email{pierre.denelle@gmail.com})
#'
#' @examples
#' data(fishmat)
#' data(fishsf)
#' 
#' net <- similarity(fishmat, metric = "Simpson")
#' clu <- netclu_greedy(net)
#' mapclu <- map_bioregions(clu, 
#'                          map = fishsf, 
#'                          map_as_output = TRUE, 
#'                          plot = FALSE)
#' 
#' # With colors
#' clu_colored <- bioregion_colors(clu)
#' mapclu <- map_bioregions(clu_colored, 
#'                          map = fishsf, 
#'                          map_as_output = TRUE, 
#'                          plot = FALSE)
#'            
#' @export
map_bioregions <- function(bioregionalization,
                           map,
                           partition_index = NULL,
                           map_as_output = FALSE,
                           plot = TRUE, 
                           clusters = NULL, 
                           geometry = NULL,
                           write_clusters = NULL,
                           ...) {
  
  # Deprecated arguments
  if (!is.null(clusters)) {
     warning(paste0("clusters is deprecated. ", 
                    "It has been replaced by bioregionalization.",
                   " Old bioregionalization argument has been replaced by",
                   " partition_index."), 
            call. = FALSE)
  }
  if (!is.null(geometry)) {
    stop("geometry is deprecated. It has been replaced by map.", 
         call. = FALSE)
  }
  if (!is.null(write_clusters)) {
    stop("write_clusters is deprecated. It has been replaced by map_as_output.", 
         call. = FALSE)
  }
  
  # Control write_clusters & plot
  controls(args = map_as_output, data = NULL, type = "boolean")
  controls(args = plot, data = NULL, type = "boolean")
  if(!map_as_output & !plot){
    stop(paste0("At least one argument among map_as_output and plot should ",
                "be set to TRUE."),
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
    partition_names <- colnames(clusters)[-1]
  }

  # Check if colors are available
  has_colors <- FALSE
  color_list <- NULL
  if (!is.null(bioregionalization$colors) && is.list(bioregionalization$colors)) {
    has_colors <- TRUE
    # Extract colors for each partition
    # colors is a list with one data.frame per partition
    color_list <- lapply(bioregionalization$colors, function(color_df) {
      # Create a named vector: cluster ID -> color
      stats::setNames(color_df$color, color_df$cluster)
    })
    names(color_list) <- names(bioregionalization$colors)
  }
  
  # Control partition_index 
  if (!is.null(partition_index)) {
    
    controls(args = partition_index, 
             data = clusters, 
             type = "input_partition_index")
    
    if(is.character(partition_index)){
      clusters <- clusters[, c(1, match(partition_index, colnames(clusters)))]
    }else{
      clusters <- clusters[, c(1, partition_index)]
    }
    partition_names <- colnames(clusters)[-1]
    
    # Update color_list if colors are present
    if (has_colors && !is.null(color_list)) {
      # Filter color_list to keep only selected partitions
      color_list <- color_list[names(color_list) %in% partition_names]
    }
    
  }
  
  # Control map
  controls(args = NULL, data = map, type = "input_map")
  
  # Convert SpatVector in sf data.frame 
  if(inherits(map, "SpatVector")){
    map <- sf::st_as_sf(map)
    map[[1]] <- as.character(map[[1]])
  }
  
  # Convert SpatRaster in sf data.frame
  if(inherits(map, "SpatRaster")){
    if (!requireNamespace("terra", quietly = TRUE)) {
      stop(paste0("The 'terra' package is required to use this function.\n", 
                  "Please install it with install.packages('terra')."),
           .call = FALSE)
    }
    pols <- terra::as.polygons(map, dissolve=FALSE, na.rm=TRUE)
    map <- sf::st_as_sf(pols)
    map[[1]] <- as.character(map[[1]])
  }
  
  # Clean map
  map <- map[, c(1, which(colnames(map) == attr(map, "sf_column")))]
  colnames(map)[1] <- "Site"
  map_site <- map$Site
  
  # Check map_site & b_site 
  missing_sites <- setdiff(b_site, map_site)
  if(length(missing_sites) > 0){
    stop(paste0("Some sites are not found in map:\n",
                "  Missing sites: ", paste(utils::head(missing_sites, 10), collapse = ", "),
                if(length(missing_sites) > 10) paste0(" ... (", length(missing_sites) - 10, " more)") else "",
                "\n  Please ensure that all sites in 'bioregionalization' have corresponding entries in 'map'."),
         call. = FALSE)
  }
  
  # Prepare geometry
  map <- map[match(b_site, map_site), ]
  map <- cbind(map, clusters[, -1])
  colnames(map)[2:(dim(clusters)[2])] <- colnames(clusters)[-1]
  
  #print(map[1:10,])
  #if(has_colors){
  #  print(color_list)
  #}
  
  # Plot
  if(plot){

    geomsp <- sf::st_geometry(map)
    plotsp <- map[, -1]
    nbplotsp <- dim(plotsp)[2]-1

    # Single partition plot
    if(nbplotsp == 1){
      
      partition_name <- partition_names[1]

      if(has_colors) {
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
        old_par <- graphics::par(mfrow = c(nrows, ncols))
        on.exit(graphics::par(old_par), add = TRUE)

        # Plot each partition with its colors
        for(i in seq_len(nbplotsp)) {
          partition_name <- partition_names[i]
          
          # Get colors for this partition
          cluster_values <- plotsp[[partition_name]]
          colors_map <- color_list[[partition_name]]
          # sf::plot expects colors in order of sorted unique values
          unique_vals <- sort(unique(cluster_values))
          pal <- colors_map[as.character(unique_vals)]
          
          # Plot with custom colors
          # Per sf documentation: when using par(mfrow), must set key.pos=NULL and reset=FALSE
          plot(plotsp[i], pal = pal, key.pos = NULL, reset = FALSE, ...)

        }
      } else {
        # No colors - use standard multi-panel layout
        mod4q <- floor(nbplotsp/4)
        mod4r <- nbplotsp - mod4q*4

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
            plot(plotsp[start_idx:end_idx], ...)
          }
        }
      }
    }
  }
  
  # Return map with added partition(s)
  if(map_as_output){
    return(map)
  }
}
