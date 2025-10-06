#' Export a network to GDF format for Gephi visualization
#'
#' This function exports a network (unipartite or bipartite) from a 
#' `data.frame` to the GDF (Graph Data Format) file format, which can be 
#' directly imported into Gephi visualization software. The function handles 
#' edge data, node attributes, and color specifications.
#'
#' @param df A two- or three-column `data.frame` where each row represents 
#' an edge (interaction) between two nodes. The first two columns contain the 
#' node identifiers, and an optional third column can contain edge weights.
#'
#' @param col1 A `character` string specifying the name of the first column 
#' in `df` containing node identifiers. Defaults to `"Node1"`.
#'
#' @param col2 A `character` string specifying the name of the second column 
#' in `df` containing node identifiers. Defaults to `"Node2"`.
#'
#' @param weight A `character` string specifying the name of the column in 
#' `df` containing edge weights. If `NULL` (default), edges are unweighted.
#'
#' @param bioregions An optional `bioregion.clusters` object (typically 
#' from clustering functions like [netclu_greedy()]) or a `data.frame` 
#' containing bioregionalization results. When a `bioregion.clusters` object with 
#' colors (from [bioregion_colors()]) is provided, colors and bioregion 
#' assignments are automatically extracted and used for visualization. 
#' Alternatively, a `data.frame` with bioregionalization data can be provided, where 
#' each row represents a node with one column containing node identifiers that 
#' match those in `df`.
#'
#' @param bioregionalization A `character` string or a positive `integer` with 
#' two different uses depending on the type of `bioregions`: 
#' \itemize{
#'   \item When `bioregions` is a `bioregion.clusters` object with multiple 
#'   partitions: specifies which partition to use. Can be either a character 
#'   string with the partition name (e.g., "K_3", "K_5") or a positive integer 
#'   indicating the partition index (e.g., 1 for first partition, 2 for second). 
#'   If `NULL` (default), the first partition is used.
#'   \item When `bioregions` is a `data.frame`: specifies the name of the column 
#'   containing node identifiers that match those in `df`. Must be a character 
#'   string. Defaults to the first column name if not specified.
#' }
#'
#' @param color_column A `character` string specifying the name of a column 
#' in `bioregions` containing color information in hexadecimal format (e.g., 
#' "#FF5733"). If specified, colors will be converted to RGB format for Gephi. 
#' If `NULL` (default), colors are automatically extracted when `bioregions` 
#' is a `bioregion.clusters` object with colors. When `bioregions` is a plain 
#' `data.frame`, this parameter must be specified to include colors.
#'
#' @param file A `character` string specifying the output file path. 
#' Defaults to `"output.gdf"`.
#'
#' @details
#' The GDF format is a simple text-based format used by Gephi to define graph 
#' structure. This function creates a GDF file with two main sections:
#' \itemize{
#'   \item \strong{nodedef}: Defines nodes and their attributes (name, label, 
#'   and any additional bioregionalization information from `bioregions`)
#'   \item \strong{edgedef}: Defines edges between nodes, optionally with 
#'   weights
#' }
#' 
#' If `color_column` is specified, hexadecimal color codes are automatically 
#' converted to RGB format (e.g., "#FF5733" becomes "255,87,51") as required 
#' by Gephi's color specification.
#' 
#' Attributes are automatically typed as VARCHAR (text), DOUBLE (numeric), 
#' or color (for color attributes).
#' 
#' \strong{Important note on zero-weight edges}: Gephi does not handle edges 
#' with weight = 0 properly. If a weight column is specified and edges with 
#' weight = 0 are detected, they will be automatically removed from the exported 
#' network, and a warning will be issued.
#'
#' @return
#' The function writes a GDF file to the specified path and returns nothing 
#' (`NULL` invisibly). The file can be directly opened in Gephi for network 
#' visualization and analysis.
#'
#' @author
#' Boris Leroy (\email{leroy.boris@gmail.com}) \cr
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) \cr
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr})
#'
#' @examples
#' # Create a simple network
#' net <- data.frame(
#'   Node1 = c("A", "A", "B", "C"),
#'   Node2 = c("B", "C", "C", "D"),
#'   Weight = c(1.5, 2.0, 1.0, 3.5)
#' )
#' 
#' # Export network with weights
#' \dontrun{
#' exportGDF(net, weight = "Weight", file = "my_network.gdf")
#' }
#' 
#' # Create bioregionalization data with colors (as data.frame)
#' bioregion_data <- data.frame(
#'   node_id = c("A", "B", "C", "D"),
#'   cluster = c("1", "2", "3", "4"),
#'   node_color = c("#FF5733", "#33FF57", "#3357FF", "#FF33F5")
#' )
#' 
#' # Export network with bioregionalization and colors
#' \dontrun{
#' exportGDF(net, 
#'           weight = "Weight",
#'           bioregions = bioregion_data,
#'           bioregionalization = "node_id",
#'           color_column = "node_color",
#'           file = "my_network_with_bioregions.gdf")
#' }
#' 
#' # Using bioregion.clusters object with colors (recommended)
#' \dontrun{
#' data(fishmat)
#' net <- similarity(fishmat, metric = "Simpson")
#' clust <- netclu_greedy(net)
#' clust_colored <- bioregion_colors(clust)
#' 
#' # Convert to network format
#' net_df <- mat_to_net(fishmat, weight = TRUE)
#' 
#' # Export with automatic colors from clustering - very simple!
#' exportGDF(net_df, 
#'           weight = "weight",
#'           bioregions = clust_colored,
#'           file = "my_network_colored.gdf")
#' 
#' # With multiple partitions, specify which one to use
#' dissim <- similarity_to_dissimilarity(similarity(fishmat, metric = "Simpson"))
#' clust_hier <- hclu_hierarclust(dissim, n_clust = c(3, 5, 8))
#' clust_hier_colored <- bioregion_colors(clust_hier)
#' 
#' # Using partition name
#' exportGDF(net_df,
#'           weight = "weight",
#'           bioregions = clust_hier_colored,
#'           bioregionalization = "K_5",
#'           file = "my_network_K5.gdf")
#' 
#' # Or using partition index (2 = second partition)
#' exportGDF(net_df,
#'           weight = "weight",
#'           bioregions = clust_hier_colored,
#'           bioregionalization = 2,
#'           file = "my_network_partition2.gdf")
#' }
#'
#' @export
exportGDF <- function(df,
                      col1 = "Node1",
                      col2 = "Node2",
                      weight = NULL,
                      bioregions = NULL,
                      bioregionalization = NULL,
                      color_column = NULL,
                      file = "output.gdf") {
  
  # Control df (network data.frame)
  controls(args = NULL, data = df, type = "input_net")
  
  # Control col1 and col2 (column names for nodes)
  controls(args = col1, data = NULL, type = "character")
  controls(args = col2, data = NULL, type = "character")
  
  if (!(col1 %in% colnames(df))) {
    stop(paste0("col1 ('", col1, "') is not a column name in df."),
         call. = FALSE)
  }
  
  if (!(col2 %in% colnames(df))) {
    stop(paste0("col2 ('", col2, "') is not a column name in df."),
         call. = FALSE)
  }
  
  # Control weight (optional column name for edge weights)
  if (!is.null(weight)) {
    controls(args = weight, data = NULL, type = "character")
    
    if (!(weight %in% colnames(df))) {
      stop(paste0("weight ('", weight, "') is not a column name in df."),
           call. = FALSE)
    }
    
    if (!is.numeric(df[[weight]])) {
      stop("The weight column must be numeric.",
           call. = FALSE)
    }
    
    if (sum(is.na(df[[weight]])) > 0) {
      stop("NA(s) detected in the weight column.",
           call. = FALSE)
    }
  }
  
  # Handle bioregion.clusters objects
  is_bioregion_clusters <- inherits(bioregions, "bioregion.clusters")
  
  if (is_bioregion_clusters) {
    # Extract information from bioregion.clusters object
    if (is.null(bioregions$clusters)) {
      stop("bioregion.clusters object does not contain clusters data.",
           call. = FALSE)
    }
    
    # Determine which partition to use
    partition_names <- names(bioregions$clusters)[-1]  # Exclude ID column
    
    if (is.null(bioregionalization)) {
      # Use first partition by default
      bioregionalization_partition <- partition_names[1]
      message(paste0("Using partition '", bioregionalization_partition, 
                    "' for cluster assignments and colors."))
    } else {
      # Control bioregionalization (character or positive integer)
      controls(args = bioregionalization, data = NULL, 
               type = "character_or_positive_integer")
      
      # Handle integer input (partition index)
      if (is.numeric(bioregionalization)) {
        if (bioregionalization > length(partition_names)) {
          stop(paste0("bioregionalization index ", bioregionalization, 
                     " is out of range. Available partitions: 1 to ",
                     length(partition_names), " (",
                     paste(partition_names, collapse = ", "), ")."),
               call. = FALSE)
        }
        bioregionalization_partition <- partition_names[bioregionalization]
        message(paste0("Using partition ", bioregionalization, " ('", 
                      bioregionalization_partition, 
                      "') for cluster assignments and colors."))
      } else {
        # Handle character input (partition name)
        if (!bioregionalization %in% partition_names) {
          stop(paste0("bioregionalization '", bioregionalization, 
                     "' not found in available partitions: ",
                     paste(partition_names, collapse = ", ")),
               call. = FALSE)
        }
        bioregionalization_partition <- bioregionalization
      }
    }
    
    # Create bioregions data.frame from clusters
    bioregions_from_clusters <- bioregions$clusters[, c(1, which(names(bioregions$clusters) == bioregionalization_partition))]
    names(bioregions_from_clusters) <- c("name", "cluster")
    
    # Add colors if available
    if (!is.null(bioregions$clusters_colors)) {
      colors_col <- bioregions$clusters_colors[[bioregionalization_partition]]
      bioregions_from_clusters$node_color <- colors_col
      
      if (is.null(color_column)) {
        color_column <- "node_color"
      }
    }
    
    # Replace bioregions with extracted data
    bioregions <- bioregions_from_clusters
    bioregionalization <- "name"
  }
  
  # Control bioregions (optional data.frame with bioregionalization data)
  if (!is.null(bioregions)) {
    if (!is_bioregion_clusters) {
      controls(args = NULL, data = bioregions, type = "input_data_frame")
    }
    
    # Set default for bioregionalization if not specified
    if (is.null(bioregionalization)) {
      bioregionalization <- colnames(bioregions)[1]
    }
    
    # Control bioregionalization (must be character for data.frame)
    if (!is_bioregion_clusters) {
      controls(args = bioregionalization, data = NULL, type = "character")
    }
    
    if (!(bioregionalization %in% colnames(bioregions))) {
      stop(paste0("bioregionalization ('", bioregionalization, 
                  "') is not a column name in bioregions."),
           call. = FALSE)
    }
  }
  
  # Control color_column (optional column name for node colors)
  if (!is.null(color_column)) {
    controls(args = color_column, data = NULL, type = "character")
  }
  
  # Control file (output file path)
  controls(args = file, data = NULL, type = "character")
  
  # Function body starts here
  # Filter out edges with weight = 0 if weight column is specified
  if (!is.null(weight)) {
    n_zero <- sum(df[[weight]] == 0)
    if (n_zero > 0) {
      warning(paste0(n_zero, " edge(s) with weight = 0 detected and removed. ",
                    "Gephi does not handle zero-weight edges properly."),
              call. = FALSE)
      df <- df[df[[weight]] != 0, ]
    }
  }
  
  nodes <- unique(c(df[[col1]], df[[col2]]))

  node_df <- data.frame(name = nodes, label = nodes, stringsAsFactors = FALSE)

  if (!is.null(bioregions)) {
    colnames(bioregions)[which(colnames(bioregions) ==
                                 bioregionalization)] <- "name"
    if (!('name' %in% colnames(bioregions))) {
      stop("The bioregions data.frame must contain a 'name' column.",
           call. = FALSE)
    }
    node_df <- merge(node_df, bioregions, by = 'name', all.x = TRUE)
  }
  
  if (!is.null(color_column)) {
    if (!color_column %in% names(node_df)) {
      stop(paste0("Color column '", color_column, 
                  "' not found in node characteristics."),
           call. = FALSE)
    }
    
    # Convert the color column from hex to RGB
    hex_to_rgb <- function(hex) {
      # Handle NA values - assign a neutral grey for nodes without colors
      # (e.g., species in bipartite networks when only sites were clustered)
      if (is.na(hex)) {
        return("200,200,200")  # Light grey for uncolored nodes
      }
      # Remove '#' if present
      hex <- gsub("#", "", hex)
      # Split into R, G, B components
      rgb_values <- sapply(seq(1, nchar(hex), 2), function(i) substr(hex, i, i+1))
      # Convert hex to decimal
      rgb_values <- as.numeric(sapply(rgb_values, function(x) strtoi(x, 16L)))
      # Return as a string "R,G,B"
      paste(rgb_values, collapse=",")
    }
    
    node_df$color <- sapply(node_df[[color_column]], hex_to_rgb)
    # Remove the original color column if not needed
    node_df[[color_column]] <- NULL
  }
  
  node_attributes <- setdiff(names(node_df), c("name", "label"))
  attribute_definitions <- sapply(node_attributes, function(attr) {
    if (attr == "color") {
      "color VARCHAR"
    } else if (is.numeric(node_df[[attr]])) {
      paste0(attr, " DOUBLE")
    } else {
      paste0(attr, " VARCHAR")
    }
  })
  
  node_header <- paste("nodedef>name VARCHAR,label VARCHAR", paste(attribute_definitions, collapse=","), sep=",")
  
  # Build node lines with proper quoting for fields containing commas
  node_lines <- apply(node_df, 1, function(x) {
    # Quote fields that contain commas (like RGB colors) with single quotes for Gephi
    quoted_fields <- sapply(x, function(field) {
      if (grepl(",", field)) {
        paste0("'", field, "'")
      } else {
        field
      }
    })
    paste(quoted_fields, collapse=",")
  })
  
  if (!is.null(weight)) {
    edge_header <- "edgedef>node1 VARCHAR,node2 VARCHAR,weight DOUBLE"
    edge_lines <- apply(df[, c(col1, col2, weight)], 1, paste, collapse=",")
  } else {
    edge_header <- "edgedef>node1 VARCHAR,node2 VARCHAR"
    edge_lines <- apply(df[, c(col1, col2)], 1, paste, collapse=",")
  }
  
  all_lines <- c(node_header, node_lines, edge_header, edge_lines)
  writeLines(all_lines, con = file)
}
