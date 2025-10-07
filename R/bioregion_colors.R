#' Add color palettes to bioregion cluster objects
#'
#' This function assigns colors to clusters in a `bioregion.clusters` object
#' using color palettes from the `rcartocolor` package. It  
#' handles large numbers of clusters by assigning vivid colors to the most
#' important clusters (based on size), grey shades to less important clusters,
#' and optionally black to insignificant clusters.
#'
#' @param clusters An object of class `bioregion.clusters`, typically output
#' from clustering functions such as [netclu_infomap()], [hclu_hierarclust()],
#' or [nhclu_pam()].
#' 
#' @param palette A `character` string indicating which color palette from
#' `rcartocolor` to use. Default is `"Vivid"`. Other qualitative palettes 
#' include `"Bold"`, `"Prism"`, `"Safe"`, `"Antique"`, and `"Pastel"`.
#' 
#' @param cluster_ordering A `character` string indicating the criterion for
#' ranking clusters to determine color assignment priority. Options are:
#' \itemize{
#'   \item `"n_sites"` (default): Rank by number of sites in each cluster
#'   \item `"n_species"`: Rank by number of species (bipartite networks only)
#'   \item `"n_both"`: Rank by combined sites + species (bipartite networks only)
#' }
#' Larger clusters (by the chosen criterion) receive vivid colors first.
#' 
#' @param cutoff_insignificant A `numeric` value or `NULL` (default). When
#' specified, clusters with values at or below this threshold (based on the
#' `cluster_ordering` criterion) are considered insignificant and colored
#' black, reducing visual clutter on maps. If `NULL`, all clusters receive
#' distinct colors.
#' 
#' @return 
#' A modified `bioregion.clusters` object with two additional elements:
#' \itemize{
#'   \item `colors`: A `list` where each element corresponds to a partition 
#'     (bioregionalization). Each list element is a `data.frame` with two columns:
#'     \itemize{
#'       \item `cluster` (`character`): Cluster identifier for that partition
#'       \item `color` (`character`): Hex color code (e.g., "#FF5733")
#'     }
#'   \item `clusters_colors`: A `data.frame` with the same structure as the
#'     `clusters` element, but with cluster IDs replaced by their corresponding
#'     hex color codes for direct use in plotting functions.
#' }
#' 
#' @details
#' The function uses a two-step algorithm to assign colors:
#' 
#' **Step 1: Identify insignificant clusters** (if `cutoff_insignificant` is specified)
#' 
#' Insignificant clusters are those with a marginal size compared to others. 
#' This is a subjective threshold set by the user. All such clusters are assigned
#' the color black (#000000) to minimize their visual impact.  
#' Clusters with values at or below the threshold are assigned black (#000000).
#' 
#' **Step 2: Assign colors to significant clusters**
#' 
#' Remaining clusters are ranked by the `cluster_ordering` criterion:
#' \itemize{
#'   \item **Top clusters** (up to 12): Receive distinct colors from the chosen
#'     palette. This limit is because above 12 the human eye struggles to 
#'     distinguish between colors.
#'   \item **Remaining clusters** (beyond top 12): Receive shades of grey from
#'     light (#CCCCCC) to dark (#404040), maintaining visual distinction but
#'     with less prominence.
#' }
#' 
#' **Multiple partitions**: If the cluster object contains multiple partitions
#' (e.g., from hierarchical clustering with different k values), colors are
#' assigned independently for each partition. Each partition gets its own color
#' scale optimized for the number of clusters in that partition.
#' 
#' @note
#' The colored cluster object can be directly used with [map_bioregions()],
#' which will automatically detect and apply the color scheme when present.
#' 
#' @seealso 
#' [map_bioregions()] for visualizing colored clusters on maps
#' 
#' 
#' @references
#' Color palettes from the `rcartocolor` package:
#' Nowosad J (2018). "CARTOColors: color palettes inspired by CARTO."
#' \url{https://github.com/Nowosad/rcartocolor}
#' 
#' @author
#' Boris Leroy (\email{leroy.boris@@gmail.com}) \cr
#' Maxime Lenormand (\email{maxime.lenormand@@inrae.fr})
#' Pierre Denelle (\email{pierre.denelle@gmail.com})
#' 
#' @examples
#' data(fishmat)
#' data(fishsf)
#' 
#' # Basic example with few clusters
#' sim <- similarity(fishmat, metric = "Simpson")
#' clust <- netclu_greedy(sim)
#' clust_colored <- bioregion_colors(clust)
#' print(clust_colored)
#' 
#' \dontrun{
#' # Map with automatic colors
#' map_bioregions(clust_colored, fishsf)
#' 
#' # Example with many clusters and cutoff
#' dissim <- similarity_to_dissimilarity(sim)
#' clust <- hclu_hierarclust(dissim,
#'                           optimal_tree_method = "best",
#'                           n_clust = 15)
#' clust_colored2 <- bioregion_colors(clust, 
#'                                    cluster_ordering = "n_sites",
#'                                    cutoff_insignificant = 1)
#' map_bioregions(clust_colored2, fishsf)
#' 
#' # Example with different palette
#' clust_colored3 <- bioregion_colors(clust, palette = "Bold")
#' map_bioregions(clust_colored3, fishsf)
#' 
#' 
#' # Example with bipartite network
#' clust_bip <- netclu_greedy(fishdf, bipartite = TRUE)
#' clust_bip_colored <- bioregion_colors(clust_bip, 
#'                                       cluster_ordering = "n_both")
#' map_bioregions(clust_bip_colored, fishsf)
#'                                       
#' }
#' 
#' @export
bioregion_colors <- function(clusters,
                             palette = "Vivid",
                             cluster_ordering = "n_sites",
                             cutoff_insignificant = NULL) {
  
  # 1. Input validation --------------------------------------------------------
  
  # Check that clusters is a bioregion.clusters object
  if (!inherits(clusters, "bioregion.clusters")) {
    stop("clusters must be a bioregion.clusters object.", call. = FALSE)
  }
  
  # Check that clusters has a clusters element
  if (is.null(clusters$clusters) || !inherits(clusters$clusters, "data.frame")) {
    stop("clusters object must contain a clusters data.frame.", call. = FALSE)
  }
  
  # Validate palette parameter
  controls(args = palette, data = NULL, type = "character")
  
  # Test if palette is valid by trying to get colors
  tryCatch({
    rcartocolor::carto_pal(n = 3, name = palette)
  }, error = function(e) {
    stop(paste0("Invalid palette name '", palette, "'. ",
                "See rcartocolor::carto_pal() for valid palette names."),
         call. = FALSE)
  })
  
  # Validate cluster_ordering parameter
  controls(args = cluster_ordering, data = NULL, type = "character")
  if (!cluster_ordering %in% c("n_sites", "n_species", "n_both")) {
    stop(paste0("cluster_ordering must be one of: 'n_sites', 'n_species', ",
                "or 'n_both'."), call. = FALSE)
  }
  
  # Check bipartite requirements for n_species and n_both
  if (cluster_ordering %in% c("n_species", "n_both")) {
    if (is.null(clusters$inputs$bipartite) || !clusters$inputs$bipartite) {
      stop(paste0("cluster_ordering '", cluster_ordering, 
                  "' can only be used with bipartite clustering."),
           call. = FALSE)
    }
  }
  
  # Validate cutoff_insignificant parameter
  if (!is.null(cutoff_insignificant)) {
    controls(args = cutoff_insignificant, data = NULL, 
             type = "positive_numeric")
    if (length(cutoff_insignificant) != 1) {
      stop("cutoff_insignificant must be a single numeric value.", 
           call. = FALSE)
    }
  }
  
  # 2. Extract cluster information ---------------------------------------------
  
  # Get cluster columns (exclude first column which is ID)
  cluster_cols <- clusters$clusters[, 2:ncol(clusters$clusters), drop = FALSE]
  nb_partitions <- ncol(cluster_cols)
  partition_names <- names(cluster_cols)
  
  # Initialize list to store colors for each partition
  colors_list <- list()
  clusters_colors <- clusters$clusters
  
  # Process each partition separately
  for (part_idx in seq_len(nb_partitions)) {
    partition_name <- partition_names[part_idx]
    partition_col <- cluster_cols[[part_idx]]
    
    # Get unique cluster IDs for this partition
    # CRITICAL: Convert everything to character
    partition_cluster_ids <- unique(as.character(partition_col))
    partition_cluster_ids <- partition_cluster_ids[!is.na(partition_cluster_ids)]
    
    # 3. Calculate cluster sizes based on cluster_ordering --------------------
    
    if (cluster_ordering == "n_sites") {
      # Count number of sites in each cluster for this partition
      cluster_sizes <- sapply(partition_cluster_ids, function(cid) {
        sum(as.character(partition_col) == cid, na.rm = TRUE)
      })
      names(cluster_sizes) <- partition_cluster_ids
      
    } else if (cluster_ordering == "n_species") {
      # For bipartite networks, count species in each cluster
      # Need to access node_type attribute
      node_type <- attr(clusters$clusters, "node_type")
      if (is.null(node_type)) {
        warning(paste0("node_type attribute not found in clusters object. ",
                       "Falling back to n_sites ordering."))
        # Recalculate with n_sites
        cluster_sizes <- sapply(partition_cluster_ids, function(cid) {
          sum(as.character(partition_col) == cid, na.rm = TRUE)
        })
        names(cluster_sizes) <- partition_cluster_ids
      } else {
        # Count only species nodes
        cluster_sizes <- sapply(partition_cluster_ids, function(cid) {
          matches <- as.character(partition_col) == cid
          sum(matches & node_type == "species", na.rm = TRUE)
        })
        names(cluster_sizes) <- partition_cluster_ids
      }
      
    } else if (cluster_ordering == "n_both") {
      # For bipartite networks, count both sites and species
      cluster_sizes <- sapply(partition_cluster_ids, function(cid) {
        sum(as.character(partition_col) == cid, na.rm = TRUE)
      })
      names(cluster_sizes) <- partition_cluster_ids
    }
    
    # 4. Separate insignificant clusters ---------------------------------------
    
    insignificant_clusters <- character(0)
    significant_clusters <- partition_cluster_ids
    
    if (!is.null(cutoff_insignificant)) {
      insignificant_clusters <- partition_cluster_ids[cluster_sizes <= cutoff_insignificant]
      significant_clusters <- partition_cluster_ids[cluster_sizes > cutoff_insignificant]
    }
    
    # 5. Rank significant clusters ---------------------------------------------
    
    if (length(significant_clusters) > 0) {
      # Sort by size descending
      ranked_clusters <- significant_clusters[order(cluster_sizes[significant_clusters], 
                                                     decreasing = TRUE)]
    } else {
      ranked_clusters <- character(0)
    }
    
    # 6. Assign colors ---------------------------------------------------------
    
    max_vivid_colors <- 12
    nb_significant <- length(ranked_clusters)
    
    # Initialize colors data frame for this partition
    partition_colors_df <- data.frame(
      cluster = character(0),
      color = character(0),
      stringsAsFactors = FALSE
    )
    
    # Assign black to insignificant clusters
    if (length(insignificant_clusters) > 0) {
      insignificant_colors <- data.frame(
        cluster = insignificant_clusters,
        color = "#000000",
        stringsAsFactors = FALSE
      )
      partition_colors_df <- rbind(partition_colors_df, insignificant_colors)
    }
    
    # Assign colors to significant clusters
    if (nb_significant > 0) {
      if (nb_significant <= max_vivid_colors) {
        # All significant clusters get vivid colors
        # rcartocolor requires at least 3 colors, so we handle small cases
        if (nb_significant < 3) {
          # For 1-2 clusters, get 3 colors and use the first ones
          all_vivid_colors <- rcartocolor::carto_pal(n = 3, name = palette)
          vivid_colors <- all_vivid_colors[1:nb_significant]
        } else {
          vivid_colors <- rcartocolor::carto_pal(n = nb_significant, name = palette)
        }
        significant_colors <- data.frame(
          cluster = ranked_clusters,
          color = vivid_colors,
          stringsAsFactors = FALSE
        )
        partition_colors_df <- rbind(partition_colors_df, significant_colors)
        
      } else {
        # Top 12 get vivid colors, rest get grey shades
        vivid_colors <- rcartocolor::carto_pal(n = max_vivid_colors, name = palette)
        top_colors <- data.frame(
          cluster = ranked_clusters[1:max_vivid_colors],
          color = vivid_colors,
          stringsAsFactors = FALSE
        )
        
        # Generate grey shades for remaining clusters
        nb_grey <- nb_significant - max_vivid_colors
        # Create grey scale from light to dark
        grey_values <- seq(from = 204, to = 64, length.out = nb_grey)
        grey_colors <- grDevices::rgb(grey_values, grey_values, grey_values, 
                                      maxColorValue = 255)
        
        remaining_colors <- data.frame(
          cluster = ranked_clusters[(max_vivid_colors + 1):nb_significant],
          color = grey_colors,
          stringsAsFactors = FALSE
        )
        
        partition_colors_df <- rbind(partition_colors_df, top_colors, remaining_colors)
      }
    }
    
    # Ensure cluster column is character
    partition_colors_df$cluster <- as.character(partition_colors_df$cluster)
    
    # Store colors for this partition
    colors_list[[partition_name]] <- partition_colors_df
    
    # 7. Create clusters_colors column for this partition ---------------------
    
    # Convert cluster IDs to character
    cluster_ids_char <- as.character(partition_col)
    
    # Replace each cluster ID with its color
    color_vec <- character(length(cluster_ids_char))
    for (i in seq_along(cluster_ids_char)) {
      if (!is.na(cluster_ids_char[i])) {
        # Find matching color
        color_match <- partition_colors_df$color[partition_colors_df$cluster == cluster_ids_char[i]]
        if (length(color_match) > 0) {
          color_vec[i] <- color_match[1]
        } else {
          # Fallback to grey if no match found (shouldn't happen)
          color_vec[i] <- "#808080"
        }
      } else {
        color_vec[i] <- NA_character_
      }
    }
    
    # Update clusters_colors for this partition
    clusters_colors[[part_idx + 1]] <- color_vec  # +1 to skip ID column
  }
  
  # 8. Add new elements to clusters object -------------------------------------
  
  clusters$colors <- colors_list
  clusters$clusters_colors <- clusters_colors
  
  # 9. Return enhanced object --------------------------------------------------
  
  return(clusters)
}
