# Helper Functions for site_species_metrics Validation Tests
#
# This script provides utility functions to create mock bioregion.clusters 
# objects for testing purposes. These mock objects mimic the structure of 
# actual clustering outputs (e.g., from netclu_infomap, netclu_greedy) 
# without requiring full clustering algorithms to run.
create_mock_clusters <- function(comat, 
                                 clusters_df, 
                                 bipartite = FALSE,
                                 species_clusters = NULL,
                                 algorithm = "mock_clustering",
                                 weight = FALSE,
                                 weight_index = 3) {
  
  # Validate inputs
  if (!is.matrix(comat) && !is.data.frame(comat)) {
    stop("comat must be a matrix or data.frame")
  }
  
  if (!is.data.frame(clusters_df)) {
    stop("clusters_df must be a data.frame")
  }
  
  if (bipartite && is.null(species_clusters)) {
    stop("species_clusters must be provided when bipartite = TRUE")
  }
  
  # Check that cluster IDs match matrix rownames
  if (!all(clusters_df[, 1] %in% rownames(comat))) {
    stop("First column of clusters_df must contain site IDs matching rownames(comat)")
  }
  
  # For bipartite, check species clusters too
  if (bipartite && !is.null(species_clusters)) {
    if (!all(species_clusters[, 1] %in% colnames(comat))) {
      stop("First column of species_clusters must contain species IDs matching colnames(comat)")
    }
  }
  
  # Determine data_type based on weight parameter
  data_type <- if (weight) "abundance" else "occurrence"
  
  # Create minimal bioregion.clusters structure
  # This mimics the output of netclu_infomap()
  mock_clusters <- list(
    name = algorithm,
    args = list(
      weight = weight,
      index = weight_index
    ), 
    inputs = list(
      bipartite = bipartite,
      weight = weight,
      data_type = data_type,  
      pairwise = !bipartite,  
      pairwise_metric = if (!bipartite) "Simpson" else NA,
      dissimilarity = !bipartite,
      nb_sites = nrow(comat),
      hierarchical = FALSE
    ),
    algorithm = list(),  
    clusters = clusters_df,
    cluster_info = data.frame(
      partition_name = colnames(clusters_df)[-1],
      n_clust = apply(clusters_df[, -1, drop = FALSE], 2, 
                     function(x) length(unique(x[!is.na(x)]))),
      stringsAsFactors = FALSE
    )
  )
  
  # Add species clusters for bipartite case
  if (bipartite && !is.null(species_clusters)) {
    
    # Check that column names match between site and species clusters
    if (!identical(colnames(clusters_df), colnames(species_clusters))) {
      stop("clusters_df and species_clusters must have identical column names")
    }
    
    # Combine site and species clusters into single data.frame
    combined_clusters <- rbind(clusters_df, species_clusters)
    
    # Add node_type attribute to distinguish sites from species
    node_types <- c(
      rep("site", nrow(clusters_df)),
      rep("species", nrow(species_clusters))
    )
    attr(combined_clusters, "node_type") <- node_types
    
    mock_clusters$clusters <- combined_clusters
    
    # Update cluster_info to reflect combined clustering
    mock_clusters$cluster_info <- data.frame(
      partition_name = colnames(combined_clusters)[-1],
      n_clust = apply(combined_clusters[, -1, drop = FALSE], 2, 
                     function(x) length(unique(x[!is.na(x)]))),
      stringsAsFactors = FALSE
    )
  }
  
  class(mock_clusters) <- c("bioregion.clusters", "list")
  
  return(mock_clusters)
}
