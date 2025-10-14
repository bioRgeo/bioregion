#' Calculate contribution metrics of sites and species
#' 
#' This function calculates metrics to assess the contribution of a given
#' species or site to its bioregion.
#' 
#' @param bioregionalization A `bioregion.clusters` object.
#' 
#' @param comat Optional (if `net` is provided). A co-occurrence `matrix` with 
#' sites as rows and species as columns. Either `comat` or `net` must be 
#' provided. If both are provided, `comat` takes precedence.
#' 
#' @param indices A `character` specifying the contribution metric to compute. 
#' Available options are `rho`, `affinity`, `fidelity`, `indicator_value` and
#' `Cz`.
#' 
#' @param net Optional (if `comat` is provided). A `data.frame` where each row 
#' represents an interaction between two nodes (site and species) and an 
#' optional third column indicating the interaction's weight. Required for `Cz` 
#' indices in bipartite networks. Either `comat` or `net` must be provided.
#' 
#' @param net_site_col A `character` (column name) or `integer` (column index) 
#' indicating the position of the column containing sites in `net`. Default is 
#' `1`. If not provided and available in `bioregionalization$args$site_col`, 
#' that value will be used.
#' 
#' @param net_species_col A `character` (column name) or `integer` (column index) 
#' indicating the position of the column containing species in `net`. Default 
#' is `2`. If not provided and available in `bioregionalization$args$species_col`, 
#' that value will be used.
#' 
#' @param net_weight A `logical` or `NULL`. When `net` is provided without `comat`, 
#' determines whether weights should be used when converting to matrix format. 
#' If `NULL` (default), weight usage is automatically determined based on the 
#' clustering algorithm inputs. Set to `TRUE` to force using weights 
#' or `FALSE` to force ignoring them.
#' 
#' @param net_weight_col A `character` (column name) or `integer` (column index) 
#' specifying which column in `net` contains the weights (i.e., abundance,
#' density or any other quantitative measure of your taxa). If `NULL` (default), 
#' the third column is used when `net_weight = TRUE`. 
#' 
#' @param verbose A `logical` indicating whether to display progress messages
#' during computation. Default is `TRUE`.
#' 
#' @return 
#' A named list of class `bioregion.site.species.metrics` containing one element 
#' per bioregionalization. Each element is a list with:
#' \itemize{
#' \item{`name`: Character, the bioregionalization name}
#' \item{`metrics`: A `data.frame` with columns `Bioregion`, `Species`, and 
#' the requested contribution metrics (e.g., `rho`, `affinity`, `fidelity`, 
#' `indicator_value`)}
#' \item{`indices`: Character vector of computed indices}
#' \item{`cz_metrics`: (If `Cz` requested) A `data.frame` with Cz metrics for 
#' nodes (sites and species), including columns `Node`, `Bioregion`, `Category`, 
#' `C` (participation coefficient), and `z` (within-bioregion degree z-score)}
#' \item{`species_clusters`: (If `Cz` requested for unipartite networks) A 
#' `data.frame` showing species assignments to bioregions based on indicator 
#' values, with columns `Species` and `Bioregion`}
#' \item{`tied_species`: (If applicable for unipartite networks) A `data.frame` 
#' with species names that had equal indicator values across multiple 
#' bioregions, including columns `Species`, `Tied_Bioregions`, and `Tie_Reason`}
#' \item{`args`: List of input arguments}
#' }
#' 
#' For a single bioregionalization, the list contains one element accessible via 
#' `result[[1]]` or `result$bioregionalization_name`.
#' 
#' For multiple bioregionalizations, access individual results via 
#' `result$bioregionalization_name` or `result[[i]]`.
#' 
#' @details
#' When multiple bioregionalizations are provided (i.e., the `clusters` arg has more 
#' than 2 columns), the function computes metrics for each bioregionalization 
#' separately. The output is always returned as a named list of class 
#' `bioregion.site.species.metrics`, regardless of the number of 
#' bioregionalizations. 
#' 
#' 
#' The \eqn{\rho} metric is derived from Lenormand et al. (2019) with the 
#' following formula:
#' 
#' \eqn{\rho_{ij} = \frac{n_{ij} - \frac{n_i n_j}{n}}{\sqrt{\left(\frac{n - n_j}{
#' n-1}\right) \left(1-\frac{n_j}{n}\right) \frac{n_i n_j}{n}}}}
#' 
#' where \eqn{n} is the number of sites, \eqn{n_i} is the number of sites in 
#' which species \eqn{i} is present, \eqn{n_j} is the number of sites in 
#' bioregion \eqn{j}, and \eqn{n_{ij}} is the number of occurrences of species 
#' \eqn{i} in sites of bioregion \eqn{j}.
#' 
#' Affinity \eqn{A}, fidelity \eqn{F}, and individual contributions
#' \eqn{IndVal} describe how species are linked to their bioregions. These
#' metrics are described in Bernardo-Madrid et al. (2019):
#' 
#' - Affinity of species to their region: 
#'   \eqn{A_i = \frac{R_i}{Z}}, where \eqn{R_i} is the occurrence/range size 
#'   of species \eqn{i} in its associated bioregion, and \eqn{Z} is the total 
#'   size (number of sites) of the bioregion. High affinity indicates that the 
#'   species occupies most sites in its bioregion.
#' 
#' - Fidelity of species to their region: 
#'   \eqn{F_i = \frac{R_i}{D_i}}, where \eqn{R_i} is the occurrence/range size 
#'   of species \eqn{i} in its bioregion, and \eqn{D_i} is its total range size. 
#'   High fidelity indicates that the species is not present in other regions.
#' 
#' - Indicator Value of species: 
#'   \eqn{IndVal = F_i \cdot A_i}.
#' 
#' `Cz` metrics are derived from Guimerà & Amaral (2005), and to be computed
#' they require both SITES and SPECIES to be assigned to bioregions.
#' 
#' For **bipartite networks**, both sites and species are already assigned to 
#' bioregions, and Cz metrics are computed directly from these assignments.
#' 
#' For **unipartite networks**, only sites are assigned to bioregions. 
#' Hence, Cz metrics can normally not be computed on unipartite networks.
#' We propose here a solution to be able to do that, but it is an experimental
#' feature which has not been extensively tested yet. To compute Cz metrics,
#' species must first be assigned to bioregions. 
#' To do that, we assign species to the bioregion where they have the highest
#' indicator value (IndVal). In case of ties:
#' \enumerate{
#'   \item Species are assigned to the bioregion where they have the highest 
#'   number of occurrences
#'   \item If still tied, species are assigned to the first bioregion in 
#'   sorted order
#' }
#' Species assignments and tie information are returned in the output for 
#' verification and reproducibility.
#' 
#' - Participation coefficient: 
#'   \eqn{C_i = 1 - \sum_{s=1}^{N_M}{\left(\frac{k_{is}}{k_i}\right)^2}}, where 
#'   \eqn{k_{is}} is the number of links of node \eqn{i} to nodes in bioregion 
#'   \eqn{s}, and \eqn{k_i} is the total degree of node \eqn{i}. A high value 
#'   means links are uniformly distributed; a low value means links are within 
#'   the node's bioregion.
#' 
#' - Within-bioregion degree z-score: 
#'   \eqn{z_i = \frac{k_i - \overline{k_{si}}}{\sigma_{k_{si}}}}, where 
#'   \eqn{k_i} is the number of links of node \eqn{i} to nodes in its bioregion 
#'   \eqn{s_i}, \eqn{\overline{k_{si}}} is the average degree of nodes in 
#'   \eqn{s_i}, and \eqn{\sigma_{k_{si}}} is the standard deviation of degrees 
#'   in \eqn{s_i}.
#' 
#' @references
#' Bernardo-Madrid R, Calatayud J, González‐Suárez M, Rosvall M, Lucas P, 
#' Antonelli A & Revilla E (2019) Human activity is altering the world’s 
#' zoogeographical regions. \emph{Ecology Letters} 22, 1297--1305.
#' 
#' Guimerà R & Amaral LAN (2005) Functional cartography of complex metabolic 
#' networks. \emph{Nature} 433, 895--900.
#' 
#' Lenormand M, Papuga G, Argagnon O, Soubeyrand M, Alleaume S & Luque S (2019)
#' Biogeographical network analysis of plant species distribution in the 
#' Mediterranean region. \emph{Ecology and Evolution} 9, 237--250.
#'
#' @seealso 
#' For more details illustrated with a practical example, 
#' see the vignette: 
#' \url{https://biorgeo.github.io/bioregion/articles/a5_3_summary_metrics.html}.
#' 
#' Associated functions: 
#' [bioregion_metrics] [bioregionalization_metrics]
#'  
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) \cr
#' Boris Leroy (\email{leroy.boris@gmail.com}) \cr
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) 
#' 
#' @examples
#' comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
#'                 20, 25)
#' rownames(comat) <- paste0("Site",1:20)
#' colnames(comat) <- paste0("Species",1:25)
#' 
#' # Single bioregionalization
#' dissim <- dissimilarity(comat, metric = "Simpson")
#' clust1 <- nhclu_kmeans(dissim, n_clust = 3, index = "Simpson")
#' 
#' metrics <- site_species_metrics(bioregionalization = clust1, 
#'                                 comat = comat,
#'                                 indices = "rho")
#' 
#' # Print shows summary
#' print(metrics)
#' 
#' # Access the data
#' metrics[[1]]$metrics  # The contribution metrics data.frame
#' head(metrics[[1]]$metrics)
#' 
#' # Multiple contribution metrics
#' net <- similarity(comat, metric = "Simpson")
#' com <- netclu_greedy(net)
#' metrics_all <- site_species_metrics(bioregionalization = com, 
#'                                     comat = comat,
#'                                     indices = c("rho", "affinity", 
#'                                                "fidelity", "indicator_value"))
#' head(metrics_all[[1]]$metrics)
#' 
#' # Multiple bioregionalizations
#' # Create two different clustering solutions
#' clust_k3 <- nhclu_kmeans(dissim, n_clust = 3, index = "Simpson")
#' clust_k4 <- nhclu_kmeans(dissim, n_clust = 4, index = "Simpson")
#' 
#' # Combine them into one bioregion.clusters object
#' clust_multi <- clust_k3
#' clust_multi$clusters <- cbind(clust_k3$clusters, K_4 = clust_k4$clusters[, 2])
#' 
#' metrics_multi <- site_species_metrics(bioregionalization = clust_multi, 
#'                                       comat = comat,
#'                                       indices = "rho")
#' 
#' # Pretty print of results
#' print(metrics_multi)
#' 
#' # Access individual bioregionalization by name
#' metrics_multi$K_3$metrics
#' # Or by index
#' metrics_multi[[1]]$metrics
#' 
#' # See structure
#' str(metrics_multi)
#' 
#' # Cz indices (bipartite networks)
#' net_bip <- mat_to_net(comat, weight = TRUE)
#' clust_bip <- netclu_greedy(net_bip, bipartite = TRUE)
#' metrics_cz <- site_species_metrics(bioregionalization = clust_bip, 
#'                                    comat = comat, 
#'                                    net = net_bip, 
#'                                    indices = "Cz")
#' metrics_cz[[1]]$cz_metrics
#' 
#' # Using net instead of comat 
#' net <- mat_to_net(comat, weight = TRUE)
#' metrics_from_net <- site_species_metrics(bioregionalization = clust1,
#'                                          net = net,
#'                                          indices = "rho")
#' 
#' # Results should be equivalent
#' metrics_from_net[[1]]$metrics
#' 
#' # Cz indices for unipartite networks
#' dissim_uni <- dissimilarity(comat, metric = "Simpson")
#' clust_uni <- nhclu_kmeans(dissim_uni, n_clust = 3)
#' net_uni <- mat_to_net(comat, weight = TRUE)
#' 
#' metrics_cz_uni <- site_species_metrics(
#'   bioregionalization = clust_uni,
#'   comat = comat,
#'   net = net_uni,
#'   indices = "Cz"
#' )
#' 
#' # View Cz metrics for sites and species
#' head(metrics_cz_uni[[1]]$cz_metrics)
#' 
#' # View species assignments to bioregions
#' head(metrics_cz_uni[[1]]$species_clusters)
#' 
#' # Check if any species had tied indicator values
#' if("tied_species" %in% names(metrics_cz_uni[[1]])) {
#'   print(metrics_cz_uni[[1]]$tied_species)
#' }
#' 
#' @export

site_species_metrics <- function(bioregionalization,
                                 comat = NULL,
                                 indices = c("rho"),
                                 net = NULL,
                                 net_site_col = 1,
                                 net_species_col = 2,
                                 net_weight = NULL,
                                 net_weight_col = NULL,
                                 verbose = FALSE){
  
  # 1. Controls ---------------------------------------------------------------
  if(verbose) {
    message("Starting site_species_metrics computation...")
  }
  
  # input can be of format bioregion.clusters
  if (inherits(bioregionalization, "bioregion.clusters")) {
    if (inherits(bioregionalization$clusters, "data.frame")) {
      has.clusters <- TRUE
      clusters <- bioregionalization$clusters
      
      # Message if multiple bioregionalizations detected
      if(ncol(clusters) > 2) {
        if(verbose) {
          message("Multiple bioregionalizations detected (", ncol(clusters) - 1, 
                  " columns). Computing metrics for each...")
        }
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
  
  # 2. Input validation and setup ----------------------------------------------
  if(verbose) {
    message("Validating input data...")
  }
  
  # Map renamed arguments to internal variable names for consistency
  site_col <- net_site_col
  species_col <- net_species_col
  weight <- net_weight
  
  # Extract site_col and species_col from bioregionalization if not provided
  if(is.null(net_site_col) && !is.null(bioregionalization$args$site_col)) {
    site_col <- bioregionalization$args$site_col
    if(verbose) {
      message("Using site_col from bioregionalization object.")
    }
  }
  if(is.null(net_species_col) && !is.null(bioregionalization$args$species_col)) {
    species_col <- bioregionalization$args$species_col
    if(verbose) {
      message("Using net_species_col from bioregionalization object.")
    }
  }
  
  # Get index/weight column name if available from bioregionalization
  user_index <- net_weight_col
  
  # For bipartite networks, use args$index (the weight column name from original net)
  # For unipartite networks, this should not be used as it refers to the similarity metric
  if(!is.null(bioregionalization$args$index) && bioregionalization$inputs$bipartite) {
    if(is.null(net_weight_col)) {
      user_index <- bioregionalization$args$index
    }
  }
  
  # Validate that at least one of comat or net is provided
  comat_provided <- !is.null(comat)
  net_provided <- !is.null(net)
  
  if (!comat_provided && !net_provided) {
    stop("Either 'comat' or 'net' must be provided.", call. = FALSE)
  }
  
  # 3. Ensure both comat AND net exist -----------------------------------------
  
  if (comat_provided && !net_provided) {
    # Scenario 1: Only comat provided - create net from comat
    if (verbose) {
      message("Only 'comat' provided. Creating 'net' from 'comat'...")
    }
    
    # Determine weight usage using priority system
    weight_info <- determine_weight_usage(
      bioregionalization = bioregionalization,
      user_weight = weight,
      user_index = user_index,
      net = NULL,
      verbose = verbose
    )
    
    # Create net from comat
    net <- mat_to_net(comat, weight = weight_info$use_weight, remove_zeroes = TRUE)
    colnames(net)[1:2] <- c("Site", "Species")
    site_col <- 1
    species_col <- 2
    
  } else if (!comat_provided && net_provided) {
    # Scenario 2: Only net provided - create comat from net
    if (verbose) {
      message("Only 'net' provided. Creating 'comat' from 'net'...")
    }
    
    # Validate net format
    if (!is.data.frame(net)) {
      stop(paste0("net should be a data.frame with at least two columns, ",
                  "corresponding to the sites and species. By default, sites are ",
                  "considered to be in the first column, and species in the second. ",
                  "This can be changed with the arguments 'net_site_col' and 'net_species_col'."),
           call. = FALSE)
    }
    
    # Validate and process site_col and species_col
    controls(args = site_col, data = NULL, type = "character_or_positive_integer")
    controls(args = species_col, data = NULL, type = "character_or_positive_integer")
    
    # Convert character column names to indices if needed
    if (is.character(site_col)) {
      if (!(site_col %in% colnames(net))) {
        stop(paste0("site_col '", site_col, "' not found in net column names."), 
             call. = FALSE)
      }
      site_col <- which(colnames(net) == site_col)
    }
    
    if (is.character(species_col)) {
      if (!(species_col %in% colnames(net))) {
        stop(paste0("species_col '", species_col, "' not found in net column names."), 
             call. = FALSE)
      }
      species_col <- which(colnames(net) == species_col)
    }
    
    if (site_col > ncol(net)) {
      stop("The site column ('net_site_col') is incorrect.", call. = FALSE)
    }
    
    if (species_col > ncol(net)) {
      stop("The species column ('net_species_col') is incorrect.", call. = FALSE)
    }
    
    # Determine weight usage using priority system
    weight_info <- determine_weight_usage(
      bioregionalization = bioregionalization,
      user_weight = weight,
      user_index = user_index,
      net = net,
      verbose = verbose
    )
    
    # Validate weight requirements
    if (weight_info$use_weight && ncol(net) < 3) {
      stop(paste0("weight is TRUE but 'net' has only ", ncol(net), 
                  " columns. A third column with weights is required."),
           call. = FALSE)
    }
    
    # Convert net to comat
    comat <- net_to_mat(
      net = net,
      weight = weight_info$use_weight,
      squared = FALSE,
      symmetrical = FALSE,
      missing_value = 0
    )
    
    # Ensure comat has sites as rows and species as columns
    if (site_col > species_col) {
      comat <- t(comat)
    }
    
  } else {
    # Scenario 3: Both comat and net provided - keep both
    if (verbose) {
      message("Both 'comat' and 'net' provided. Using both.")
    }
    
    # TODO: Add validation for consistency between comat and net
    # For now, just determine weight usage
    weight_info <- determine_weight_usage(
      bioregionalization = bioregionalization,
      user_weight = weight,
      user_index = user_index,
      net = net,
      verbose = verbose
    )
  }
  
  # Update weight variable for later use
  weight <- weight_info$use_weight
  
  # Now validate comat
  controls(args = NULL, data = comat, type = "input_matrix")
  
  # Check if rownames of comat are present in clusters
  if(!all(rownames(comat) %in% clusters[, 1])){
    stop(paste0("Site names of comat cannot be found in clusters of object ",
                "bioregionalization.\n",
                "Did you invert column/row names? Sites should be in rows and ",
                "species in columns in comat."),
         call. = FALSE)
  }
  
  controls(args = indices, data = NULL, type = "character_vector")
  
  if(!isTRUE(unique(indices %in% c("rho", "Cz", "affinity", "fidelity",
                                   "indicator_value")))){
    stop(paste0("Please choose indices from the following:\n",
                "rho, affinity, fidelity, indicator_value or Cz."),
         call. = FALSE)
  }
  
  # 4. Validate indices and prepare for computation ----------------------------
  
  controls(args = indices, data = NULL, type = "character_vector")
  
  if(!isTRUE(unique(indices %in% c("rho", "Cz", "affinity", "fidelity",
                                   "indicator_value")))){
    stop(paste0("Please choose indices from the following:\n",
                "rho, affinity, fidelity, indicator_value or Cz."),
         call. = FALSE)
  }
  
  # At this point, both comat and net are guaranteed to exist and be properly formatted
  # Validate net structure (basic checks)
  if(!is.data.frame(net)){
    stop("net should be a data.frame with at least two columns,
         corresponding to the sites and species. By default, sites are
         considered to be in the first column, and species in the second.
         This can be changed with the arguments 'net_site_col' and
         'net_species_col'.")
  }
  
  # 5. Main computation ---------------------------------------------------------
  if(verbose) {
    message("Input validation completed. Starting main computation...")
  }
  
  # Store bipartite status
  is_bipartite <- bioregionalization$inputs$bipartite
  
  # Inform user if computing Cz for unipartite network
  if("Cz" %in% indices && !is_bipartite) {
    if(verbose) {
      message("Computing Cz indices for UNIPARTITE network.")
      message("Species will be assigned to bioregions based on indicator values.")
      message("Note: This feature assigns each species to the bioregion where it has the highest indicator value.")
      message("This is an experimental feature which has not been extensively tested.")
    }
  }
  
  # Handle multiple bioregionalizations
  if(ncol(clusters) > 2) {
    # Multiple bioregionalizations
    bioregion_names <- colnames(clusters)[-1]
    
    results_list <- lapply(2:ncol(clusters), function(col_idx) {
      
      if(verbose) {
        message("Computing metrics for bioregionalization ", col_idx - 1, 
                " of ", ncol(clusters) - 1, ": ", bioregion_names[col_idx - 1])
      }
      
      # Extract single bioregionalization
      single_clusters <- clusters[, c(1, col_idx)]
      colnames(single_clusters) <- c(colnames(clusters)[1], "cluster")
      
      # Preserve node_type attribute for bipartite cases
      if(is_bipartite) {
        attr(single_clusters, "node_type") <- attr(clusters, "node_type")
      }
      
      # Make a copy of net for this bioregionalization (if needed)
      net_copy <- if(!is.null(net)) net else NULL
      
      # Compute metrics for this bioregionalization
      compute_single_bioregionalization_metrics(
        single_clusters = single_clusters,
        bioregion_name = bioregion_names[col_idx - 1],
        weight = bioregionalization$args$weight,
        weight_index = bioregionalization$args$index,
        comat = comat,
        indices = indices,
        net = net_copy,
        site_col = site_col,
        species_col = species_col,
        is_bipartite = is_bipartite,
        verbose = verbose
      )
    })
    
    names(results_list) <- bioregion_names
    
  } else {
    # Single bioregionalization - also return list for consistency
    bioregion_name <- colnames(clusters)[2]
    
    if(verbose) {
      message("Computing metrics for single bioregionalization: ", bioregion_name)
    }
    
    results_list <- list(
      compute_single_bioregionalization_metrics(
        single_clusters = clusters,
        bioregion_name = bioregion_name,
        weight = bioregionalization$args$weight,
        weight_index = bioregionalization$args$index,
        comat = comat,
        indices = indices,
        net = net,
        site_col = site_col,
        species_col = species_col,
        is_bipartite = is_bipartite,
        verbose = verbose
      )
    )
    names(results_list) <- bioregion_name
  }
  
  # Add class and attributes for S3 method dispatch
  class(results_list) <- c("bioregion.site.species.metrics", "list")
  attr(results_list, "indices") <- indices
  attr(results_list, "n_bioregionalizations") <- length(results_list)
  
  if(verbose) {
    message("Computation completed successfully!")
  }
  
  return(results_list)
}


# Helper function to compute metrics for a single bioregionalization
#' @noRd
compute_single_bioregionalization_metrics <- function(single_clusters,
                                                      bioregion_name,
                                                      weight,
                                                      weight_index,
                                                      comat,
                                                      indices,
                                                      net = NULL,
                                                      site_col = 1,
                                                      species_col = 2,
                                                      is_bipartite = FALSE,
                                                      verbose = TRUE) {
  
  rho_df <- affinity_df <- fidelity_df <- indval_df <- bipartite_df <- NULL
  

  weight <- ifelse(is.null(weight), FALSE, weight)
  
  # 1. Cz indices ------------------------------------------
  if("Cz" %in% indices) {
    if(verbose) {
      message("  Computing Cz indices (participation coefficient and z-score)...")
    }
    
    # Store original column name for weight if it exists
    # weight_index can be either a character (column name) or numeric (column position)
    original_weight_colname <- NULL
    if(weight && !is.null(weight_index)) {
      if(is.character(weight_index)) {
        # weight_index is a column name
        original_weight_colname <- weight_index
      } else if(is.numeric(weight_index) && weight_index <= ncol(net)) {
        # weight_index is a column position
        original_weight_colname <- colnames(net)[weight_index]
      }
    }
    
    # Auto-detect weight from net structure if not explicitly specified
    # (happens when user provides net directly without going through a clustering function)
    if(is.null(original_weight_colname) && ncol(net) >= 3) {
      # Net has 3+ columns, assume it's weighted
      if(!weight) {
        # Update weight to TRUE if net has weights but weight was FALSE
        weight <- TRUE
        if(verbose) {
          message("  Auto-detected weighted net (3+ columns). Setting weight = TRUE for Cz calculations.")
        }
      }
      # Assume third column is weights
      original_weight_colname <- colnames(net)[3]
      if(verbose) {
        message("  Using column '", original_weight_colname, "' as weights for Cz calculations.")
      }
    }
    
    # If net has only 2 columns, it's unweighted
    if(ncol(net) < 3 && weight) {
      if(verbose) {
        message("  Note: weight = TRUE but net has only ", ncol(net), " columns. Using unweighted calculations.")
      }
      weight <- FALSE
    }
    
    # Rename site and species columns as Sites and Species
    colnames(net)[site_col] <- "Site"
    colnames(net)[species_col] <- "Species"
    
    # Handle bipartite vs unipartite differently
    if(is_bipartite) {
      bipartite_df <- single_clusters
      # Add a column category (site or species) to bipartite_df
      bipartite_df$cat <- attributes(single_clusters)$node_type
      colnames(bipartite_df) <- c("Node", "Bioregion", "Category")
      
    } else {
  
      # Get bioregions
      bioregions <- sort(unique(single_clusters[, 2]))
      
      # Assign species using helper function
      species_assignment <- assign_species_to_bioregions(
        comat = comat,
        single_clusters = single_clusters,
        bioregions = bioregions
      )
      
      # Create bipartite_df structure
      # Sites
      site_df <- data.frame(
        Node = single_clusters[, 1],
        Bioregion = single_clusters[, 2],
        Category = "site",
        stringsAsFactors = FALSE
      )
      
      # Species
      species_df <- data.frame(
        Node = species_assignment$species_clusters$Species,
        Bioregion = species_assignment$species_clusters$Bioregion,
        Category = "species",
        stringsAsFactors = FALSE
      )
      
      # Combine
      bipartite_df <- rbind(site_df, species_df)
      
      # Store tied species info for output
      attr(bipartite_df, "tied_species") <- species_assignment$tied_species
      
      if(nrow(species_assignment$tied_species) > 0 && verbose) {
        message("    ", nrow(species_assignment$tied_species), 
                " species had tied indicator values (see $tied_species for details).")
      }
    }
    
    # Add bioregions of the sites to the bipartite data.frame
    net$Site <- as.character(net$Site)
    bipartite_df$Node <- as.character(bipartite_df$Node)
    
    site_bioregions <- bipartite_df[, c("Node", "Bioregion")]
    colnames(site_bioregions) <- c("Site", "Bioregion_site")
    net <- merge(net, site_bioregions, by = "Site", all.x = TRUE)
    
    # Add bioregions of the species to the bipartite data.frame
    net$Species <- as.character(net$Species)
    species_bioregions <- bipartite_df[, c("Node", "Bioregion")]
    colnames(species_bioregions) <- c("Species", "Bioregion_species")
    net <- merge(net, species_bioregions, by = "Species", all.x = TRUE)
    
    if(weight && !is.null(original_weight_colname)) {
      weight_col_idx <- which(colnames(net) == original_weight_colname)
      if(length(weight_col_idx) == 0) {
        warning("Weight column '", original_weight_colname, "' not found after merge. Using unweighted calculations.")
        weight <- FALSE
      } else {
        weight_col_idx <- weight_col_idx[1]
      }
    }
    
    # Compute coefficient of participation C
    if(verbose) {
      message("    Computing participation coefficient C for sites and species...")
    }

    # Vectorized computation for sites
    # For sites: C measures distribution of (weighted) links across species bioregions
    dat_com <- bipartite_df[which(bipartite_df$Category == "site"), ]
    
    if(weight) {
      # Sum weights by site and species bioregion
      C_site <- by(net, net$Site, function(site_net) {
        weight_by_bioregion <- tapply(site_net[[weight_col_idx]], 
                                      site_net$Bioregion_species, 
                                      sum, na.rm = TRUE)
        total_weight <- sum(weight_by_bioregion, na.rm = TRUE)
        if(total_weight == 0) return(NA)
        1 - sum((weight_by_bioregion/total_weight)^2)
      })
    } else {
      # Count occurrences by site and species bioregion
      C_site <- tapply(net$Bioregion_species, net$Site, function(x) {
        tmp <- table(x)
        1 - sum((tmp/sum(tmp))^2)
      })
    }
    
    C_site <- data.frame(
      Node = names(C_site),
      C = as.numeric(C_site),
      Category = "site",
      stringsAsFactors = FALSE
    )
    
    # Vectorized computation for species
    # For species: C measures distribution of (weighted) links across site bioregions
    dat_sp <- bipartite_df[which(bipartite_df$Category == "species"), ]
    
    if(weight) {
      # Sum weights by species and site bioregion
      C_sp <- by(net, net$Species, function(species_net) {
        weight_by_bioregion <- tapply(species_net[[weight_col_idx]], 
                                      species_net$Bioregion_site, 
                                      sum, na.rm = TRUE)
        total_weight <- sum(weight_by_bioregion, na.rm = TRUE)
        if(total_weight == 0) return(NA)
        1 - sum((weight_by_bioregion/total_weight)^2)
      })
    } else {
      # Count occurrences by species and site bioregion
      C_sp <- tapply(net$Bioregion_site, net$Species, function(x) {
        tmp <- table(x)
        1 - sum((tmp/sum(tmp))^2)
      })
    }
    
    C_sp <- data.frame(
      Node = names(C_sp),
      C = as.numeric(C_sp),
      Category = "species",
      stringsAsFactors = FALSE
    )
    
    C_dat <- rbind(C_site, C_sp)
    
    # Merge results with bipartite_df using base R
    bipartite_df <- merge(bipartite_df, C_dat, by = c("Node", "Category"), 
                         all.x = TRUE)
    
    # Compute z
    if(verbose) {
      message("    Computing within-bioregion degree z-score...")
    }
    
    # For sites: (weighted) sum of links to species that are in the SITE'S bioregion
    # For species: (weighted) sum of links to sites that are in the SPECIES' bioregion
    # In code terms: count/sum links where the bioregion of the node
    # matches the bioregion of its partner

    if(weight) {
      site_links <- sapply(dat_com$Node, function(site_node) {
        site_bioregion <- dat_com$Bioregion[dat_com$Node == site_node]
        site_net <- net[net$Site == site_node, ]
        # Sum weights where species is in same bioregion
        sum(site_net[[weight_col_idx]][site_net$Bioregion_species == site_bioregion], 
            na.rm = TRUE)
      })
    } else {
      site_links <- sapply(dat_com$Node, function(site_node) {
        site_bioregion <- dat_com$Bioregion[dat_com$Node == site_node]
        site_net <- net[net$Site == site_node, ]
        sum(site_net$Bioregion_species == site_bioregion)
      })
    }
    
    site_links_df <- data.frame(
      Node = dat_com$Node,
      n_link_bioregion = as.numeric(site_links),
      stringsAsFactors = FALSE
    )
    
    # Count/sum within-bioregion links for species
    # (links to sites in the same bioregion as the species)
    if(weight) {
      species_links <- sapply(dat_sp$Node, function(species_node) {
        species_bioregion <- dat_sp$Bioregion[dat_sp$Node == species_node]
        species_net <- net[net$Species == species_node, ]
        # Sum weights where species is in same bioregion
        sum(species_net[[weight_col_idx]][species_net$Bioregion_site == species_bioregion], 
            na.rm = TRUE)
      })
    } else {
      species_links <- sapply(dat_sp$Node, function(species_node) {
        species_bioregion <- dat_sp$Bioregion[dat_sp$Node == species_node]
        species_net <- net[net$Species == species_node, ]
        sum(species_net$Bioregion_site == species_bioregion)
      })
    }
    
    species_links_df <- data.frame(
      Node = dat_sp$Node,
      n_link_bioregion = as.numeric(species_links),
      stringsAsFactors = FALSE
    )
    
    # Merge with bipartite_df
    all_links <- rbind(site_links_df, species_links_df)
    bipartite_df <- merge(bipartite_df, all_links, by = "Node", all.x = TRUE)
    
    # Calculate mean and SD separately for sites and species within each bioregion
    # Create a grouping key
    bipartite_df$group_key <- paste(bipartite_df$Bioregion, 
                                     bipartite_df$Category, sep = "_")
    
    # Compute mean by group
    mean_by_group <- tapply(bipartite_df$n_link_bioregion,
                           bipartite_df$group_key,
                           mean, na.rm = TRUE)
    mean_df <- data.frame(
      group_key = names(mean_by_group),
      mean_link_bioregion = as.numeric(mean_by_group),
      stringsAsFactors = FALSE
    )
    
    # Compute SD by group
    sd_by_group <- tapply(bipartite_df$n_link_bioregion,
                         bipartite_df$group_key,
                         stats::sd, na.rm = TRUE)
    sd_df <- data.frame(
      group_key = names(sd_by_group),
      sd_link_bioregion = as.numeric(sd_by_group),
      stringsAsFactors = FALSE
    )
    
    # Merge back to bipartite_df
    bipartite_df <- merge(bipartite_df, mean_df, by = "group_key", all.x = TRUE)
    bipartite_df <- merge(bipartite_df, sd_df, by = "group_key", all.x = TRUE)
    
    # Calculate z-score
    bipartite_df$z <- (bipartite_df$n_link_bioregion -
                         bipartite_df$mean_link_bioregion) /
      bipartite_df$sd_link_bioregion
    
    # Replace NaN with 0 when SD is 0 (all nodes have same degree)
    bipartite_df$z[is.nan(bipartite_df$z)] <- 0
    
    # Remove intermediate columns
    bipartite_df <- bipartite_df[, c("Node", "Bioregion", "Category", "C", "z")]
  }
  
  # 2. Contribution metrics (rho, affinity, fidelity, indicator_value) ------
  if(any(c("rho", "affinity", "fidelity", "indicator_value") %in% indices)){
    if(verbose) {
      selected_indices <- intersect(indices, 
                                   c("rho", "affinity", "fidelity", "indicator_value"))
      message("  Computing contribution metrics: ", 
              paste(selected_indices, collapse = ", "))
    }
    # Binary site-species matrix
    comat_bin <- comat
    comat_bin[comat_bin > 0] <- 1
    
    # If it is a bipartite object, split site and species into two objects
     if(is_bipartite){
       all_clusters <- single_clusters
       single_clusters <- single_clusters[
         which(attributes(single_clusters)$node_type == "site"), ]
     }
    
    # Data.frames with output
    rho_df <- data.frame(Bioregion = character(),
                         Species = character(),
                         rho = character())
    affinity_df <- data.frame(Bioregion = character(),
                              Species = character(),
                              affinity = character())
    fidelity_df <- data.frame(Bioregion = character(),
                              Species = character(),
                              fidelity = character())
    indval_df <- data.frame(Bioregion = character(),
                            Species = character(),
                            indval = character())
    
    # Formula components
    n <- nrow(comat) # number of sites
    rich <- rowSums(comat_bin) # species richness per site
    n_i <- colSums(comat_bin) # number of occurrences per species (previously n_i)
    n_species <- ncol(comat)
    
    # Get unique bioregions
    bioregions <- unique(single_clusters[, 2])
    n_bioregions <- length(bioregions)
    
    if(verbose) {
      message("    Processing ", n_bioregions, " bioregions...")
    }
    
    # Vectorized computation over bioregions
    # Pre-allocate storage lists for efficiency
    rho_list <- affinity_list <- fidelity_list <- indval_list <- vector("list", n_bioregions)
    
    # Create a mapping matrix: sites x bioregions (TRUE if site belongs to bioregion)
    # This avoids repeated subsetting operations
    site_bioregion_map <- sapply(bioregions, function(br) {
      single_clusters[, 2] == br
    })
    rownames(site_bioregion_map) <- single_clusters[, 1]
    
    # Count sites per bioregion
    n_j_vec <- colSums(site_bioregion_map)
    
    # Compute n_ij matrix: species (rows) x bioregions (columns)
    # n_ij[i, j] = number of occurrences of species i in sites of bioregion j
    # Matrix multiplication: t(comat_bin) %*% site_bioregion_map gives us this directly
    n_ij_mat <- t(comat_bin) %*% site_bioregion_map
    
    # Now vectorize calculations across all bioregions at once
    if("rho" %in% indices){
      # Compute rho for all species x bioregion combinations at once
      # Original formula: p_ij <- (n_ij - ((n_i*n_j)/n))/(sqrt(((n - n_j)/(n-1))*(1-(n_j/n))*((n_i*n_j)/n)))
      # where:
      #   n_ij = occurrences of species i in bioregion j (matrix: species x bioregions)
      #   n_i = total occurrences of species i across all sites (vector: per species)
      #   n_j = number of sites in bioregion j (vector: per bioregion)
      #   n = total number of sites (scalar)
      
      # Replicate vectors to match n_ij_mat dimensions (species x bioregions)
      # n_i_mat: replicate n_i (per species) across bioregions (columns)
      n_i_mat <- matrix(n_i, nrow = n_species, ncol = n_bioregions, byrow = FALSE)
      # n_j_mat: replicate n_j_vec (per bioregion) across species (rows)
      n_j_mat <- matrix(n_j_vec, nrow = n_species, ncol = n_bioregions, byrow = TRUE)
      
      # Calculate rho using the original formula, vectorized
      numerator <- n_ij_mat - ((n_i_mat * n_j_mat) / n)
      denominator <- sqrt(((n - n_j_mat) / (n - 1)) * (1 - (n_j_mat / n)) * ((n_i_mat * n_j_mat) / n))
      p_ij_mat <- numerator / denominator
      
      # Convert matrix to long format data.frame
      rho_df <- data.frame(
        Bioregion = rep(as.character(bioregions), each = n_species),
        Species = rep(colnames(comat), times = n_bioregions),
        rho = as.vector(p_ij_mat)
      )
    }
    
    if("affinity" %in% indices) {
      # Affinity: n_ij / n_j for each bioregion
      # Divide each column of n_ij_mat by corresponding n_j
      affinity_mat <- sweep(n_ij_mat, 2, n_j_vec, FUN = "/")
      
      affinity_df <- data.frame(
        Bioregion = rep(as.character(bioregions), each = n_species),
        Species = rep(colnames(comat), times = n_bioregions),
        affinity = as.vector(affinity_mat)
      )
    }
    
    if("fidelity" %in% indices) {
      # Fidelity: n_ij / n_i for each species
      # Divide each row of n_ij_mat by corresponding n_i
      fidelity_mat <- sweep(n_ij_mat, 1, n_i, FUN = "/")
      
      fidelity_df <- data.frame(
        Bioregion = rep(as.character(bioregions), each = n_species),
        Species = rep(colnames(comat), times = n_bioregions),
        fidelity = as.vector(fidelity_mat)
      )
    }
    
    if("indicator_value" %in% indices) {
      # Indicator value: (n_ij/n_j) * (n_ij/n_i)
      affinity_mat <- sweep(n_ij_mat, 2, n_j_vec, FUN = "/")
      fidelity_mat <- sweep(n_ij_mat, 1, n_i, FUN = "/")
      indval_mat <- affinity_mat * fidelity_mat
      
      indval_df <- data.frame(
        Bioregion = rep(as.character(bioregions), each = n_species),
        Species = rep(colnames(comat), times = n_bioregions),
        indval = as.vector(indval_mat)
      )
    }
    
    # Merge all outputs together into a single data.frame
    if(verbose) {
      message("    Merging results...")
    }
    res_df <- dplyr::full_join(rho_df, affinity_df, 
                               by = c("Bioregion", "Species"))
    res_df <- dplyr::full_join(res_df, fidelity_df, 
                               by = c("Bioregion", "Species"))
    res_df <- dplyr::full_join(res_df, indval_df, 
                               by = c("Bioregion", "Species"))
    
    # Remove columns full of NAs (indices not selected)
    res_df <- res_df[, colSums(is.na(res_df)) != nrow(res_df)]
    
    # Controls on the output
    # test if all bioregions are there
    if(length(unique(res_df$Bioregion)) != n_bioregions){
      warning("Not all bioregions are in the output.")
    }
    
    # test if all species are there
    if(length(unique(res_df$Species)) != ncol(comat)){
      warning("Not all species are in the output.")
    }
    
    # test if all species are there X times (X = nb of bioregions)
    if(length(unique(table(res_df$Species))) != 1 ||
       unique(table(res_df$Species)) != n_bioregions){
      warning("Not all species x bioregions combinations are in the output.")
    }
  } else {
    res_df <- NULL
  }
  
  # Return structured list
  output <- list(
    name = bioregion_name,
    metrics = res_df,
    indices = indices,
    args = list(
      site_col = site_col,
      species_col = species_col
    )
  )
  
  # If Cz indices were computed, add them
  if("Cz" %in% indices && !is.null(bipartite_df)) {
    output$cz_metrics <- bipartite_df
    
    # For unipartite Cz calculations, add species assignments and tie information
    if(!is_bipartite) {
      # Extract species clusters
      species_clusters <- bipartite_df[bipartite_df$Category == "species", 
                                        c("Node", "Bioregion")]
      colnames(species_clusters) <- c("Species", "Bioregion")
      output$species_clusters <- species_clusters
      
      # Add tied species info if it exists
      tied_species <- attr(bipartite_df, "tied_species")
      if(!is.null(tied_species) && nrow(tied_species) > 0) {
        output$tied_species <- tied_species
      }
    }
  }
  
  return(output)
}




#' Assign species to bioregions based on indicator values
#' 
#' For unipartite clusterings, assign each species to the bioregion where it 
#' has the highest indicator value. Ties are broken by highest occurrence, 
#' then by first bioregion in sorted order.
#' 
#' @param comat Site x species matrix
#' @param single_clusters Data frame with site assignments (ID and cluster)
#' @param bioregions Vector of unique bioregion IDs
#' @param tie_method Character specifying tie-breaking method: 
#'   "occurrence" (default) uses highest occurrence,
#'   "first" uses first bioregion in sorted order
#' 
#' @return List containing:
#'   - species_clusters: data.frame with Species and Bioregion columns
#'   - tied_species: data.frame with species names and tied bioregions
#' 
#' @noRd
assign_species_to_bioregions <- function(comat, single_clusters, bioregions,
                                          tie_method = "occurrence") {
  # TODO: vectorize loop for efficiency
  # Binary matrix
  comat_bin <- comat
  comat_bin[comat_bin > 0] <- 1
  
  n_species <- ncol(comat)
  species_names <- colnames(comat)
  n_bioregions <- length(bioregions)
  
  # Calculate indicator values for all species x bioregion combinations
  
  # Create site-bioregion mapping matrix
  site_bioregion_map <- sapply(bioregions, function(br) {
    single_clusters[, 2] == br
  })
  rownames(site_bioregion_map) <- single_clusters[, 1]
  
  # Count sites per bioregion
  n_j_vec <- colSums(site_bioregion_map)
  
  # Count species occurrences per bioregion
  n_ij_mat <- t(comat_bin) %*% site_bioregion_map
  
  # Total occurrences per species
  n_i <- colSums(comat_bin)
  
  # Replicate for matrix operations
  n_i_mat <- matrix(n_i, nrow = n_species, ncol = n_bioregions, byrow = FALSE)
  n_j_mat <- matrix(n_j_vec, nrow = n_species, ncol = n_bioregions, byrow = TRUE)
  
  # Calculate indicator value: affinity * fidelity
  affinity_mat <- n_ij_mat / n_j_mat
  fidelity_mat <- n_ij_mat / n_i_mat
  indval_mat <- affinity_mat * fidelity_mat
  
  # Replace NaN with 0 (occurs when species has 0 occurrences)
  indval_mat[is.nan(indval_mat)] <- 0
  
  # For each species, find bioregion with highest indval
  max_indval_idx <- apply(indval_mat, 1, which.max)
  max_indval <- apply(indval_mat, 1, max)
  
  # Identify ties: species where multiple bioregions have max indval
  tied_species_list <- list()
  assigned_bioregions <- bioregions[max_indval_idx]
  
  for(i in 1:n_species) {
    # Find all bioregions with max indval for this species
    max_val <- max_indval[i]
    bioregions_with_max <- which(abs(indval_mat[i, ] - max_val) < 1e-10)
    
    if(length(bioregions_with_max) > 1) {
      # Tie exists - resolve by specified method
      if(tie_method == "occurrence") {
        # Resolve by highest occurrence
        occurrences_in_tied <- n_ij_mat[i, bioregions_with_max]
        max_occurrence <- max(occurrences_in_tied)
        bioregions_with_max_occ <- bioregions_with_max[
          occurrences_in_tied == max_occurrence
        ]
        
        if(length(bioregions_with_max_occ) > 1) {
          # Still tied - use first bioregion in sorted order
          tied_bioregions <- bioregions[bioregions_with_max]
          tied_species_list[[species_names[i]]] <- data.frame(
            Species = species_names[i],
            Tied_Bioregions = paste(tied_bioregions, collapse = ", "),
            Tie_Reason = "Equal indval and occurrence",
            stringsAsFactors = FALSE
          )
          assigned_bioregions[i] <- bioregions[bioregions_with_max_occ[1]]
        } else {
          # Resolved by occurrence
          tied_bioregions <- bioregions[bioregions_with_max]
          tied_species_list[[species_names[i]]] <- data.frame(
            Species = species_names[i],
            Tied_Bioregions = paste(tied_bioregions, collapse = ", "),
            Tie_Reason = "Equal indval (resolved by occurrence)",
            stringsAsFactors = FALSE
          )
          assigned_bioregions[i] <- bioregions[bioregions_with_max_occ]
        }
      } else {
        # Use first bioregion in sorted order
        tied_bioregions <- bioregions[bioregions_with_max]
        tied_species_list[[species_names[i]]] <- data.frame(
          Species = species_names[i],
          Tied_Bioregions = paste(tied_bioregions, collapse = ", "),
          Tie_Reason = "Equal indval",
          stringsAsFactors = FALSE
        )
        assigned_bioregions[i] <- bioregions[bioregions_with_max[1]]
      }
    }
  }
  
  # Create output data.frame
  species_clusters <- data.frame(
    Species = species_names,
    Bioregion = assigned_bioregions,
    stringsAsFactors = FALSE
  )
  
  # Combine tied species info
  tied_species_df <- if(length(tied_species_list) > 0) {
    do.call(rbind, tied_species_list)
  } else {
    data.frame(
      Species = character(),
      Tied_Bioregions = character(),
      Tie_Reason = character(),
      stringsAsFactors = FALSE
    )
  }
  
  return(list(
    species_clusters = species_clusters,
    tied_species = tied_species_df
  ))
}

