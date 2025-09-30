#' Calculate contribution metrics of sites and species
#' 
#' This function calculates metrics to assess the contribution of a given
#' species or site to its bioregion.
#' 
#' @param bioregionalization A `bioregion.clusters` object.
#' 
#' @param comat A co-occurrence `matrix` with sites as rows and species as 
#' columns. 
#' 
#' @param indices A `character` specifying the contribution metric to compute. 
#' Available options are `rho`, `affinity`, `fidelity`, `indicator_value` and
#' `Cz`.
#' 
#' @param net `NULL` by default. Required for `Cz` indices. A 
#' `data.frame` where each row represents an interaction between two nodes 
#' and an optional third column indicating the interaction's weight.  
#' 
#' @param site_col A number indicating the position of the column containing
#' the sites in `net`. 1 by default.
#' 
#' @param species_col A number indicating the position of the column
#' containing the species in `net`. 2 by default.
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
#' \item{`cz_metrics`: (If `Cz` requested) A `data.frame` with Cz metrics for nodes}
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
#' `Cz` metrics are derived from Guimerà & Amaral (2005):
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
#' # Cz indices (bipartite only)
#' net_bip <- mat_to_net(comat, weight = TRUE)
#' clust_bip <- netclu_greedy(net_bip, bipartite = TRUE)
#' metrics_cz <- site_species_metrics(bioregionalization = clust_bip, 
#'                                    comat = comat, 
#'                                    net = net_bip, 
#'                                    indices = "Cz")
#' metrics_cz[[1]]$cz_metrics
#' 
#' @export


site_species_metrics <- function(bioregionalization,
                                 comat,
                                 indices = c("rho"),
                                 net = NULL,
                                 site_col = 1,
                                 species_col = 2){
  
  # 1. Controls ---------------------------------------------------------------
  # input can be of format bioregion.clusters
  if (inherits(bioregionalization, "bioregion.clusters")) {
    if (inherits(bioregionalization$clusters, "data.frame")) {
      has.clusters <- TRUE
      clusters <- bioregionalization$clusters
      
      # Message if multiple bioregionalizations detected
      if(ncol(clusters) > 2) {
        message("Multiple bioregionalizations detected (", ncol(clusters) - 1, 
                " columns). Computing metrics for each...")
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
  
  controls(args = indices, data = NULL, type = "character_vector")
  
  if(!isTRUE(unique(indices %in% c("rho", "Cz", "affinity", "fidelity",
                                   "indicator_value")))){
    stop(paste0("Please choose indices from the following:\n",
                "rho, affinity, fidelity, indicator_value or Cz."),
         call. = FALSE)
  }
  
  if("Cz" %in% indices && is.null(net)){
    stop("net is needed to compute Cz indices.",
         call. = FALSE)
  }
  
  if("Cz" %in% indices && bioregionalization$inputs$bipartite == FALSE){
    stop(paste0("Cz metrics can only be computed for a bipartite ",
                "bioregionalization ",
                "(where both sites and species are assigned to a bioregion."),
         call. = FALSE)
  }
  
  if(!is.null(net)){
    if(!is.data.frame(net)){
      stop("net should be a data.frame with at least two columns,
           corresponding to the sites and species. By default, sites are
           considered to be in the first column, and species in the second.
           This can be changed with the arguments 'site_col' and
           'species_col'.")
    }
    controls(args = site_col, data = NULL, type = "strict_positive_numeric")
    controls(args = species_col, data = NULL, type = "strict_positive_numeric")
    
    if(site_col > ncol(net)){
      stop("The site column ('site_col') is incorrect.")
    }
    
    if(species_col > ncol(net)){
      stop("The species column ('species_col') is incorrect.")
    }
  }
  
  # 2. Main computation -------------------------------------------------------
  # Store bipartite status
  is_bipartite <- bioregionalization$inputs$bipartite
  
  # Handle multiple bioregionalizations
  if(ncol(clusters) > 2) {
    # Multiple bioregionalizations
    bioregion_names <- colnames(clusters)[-1]
    
    results_list <- lapply(2:ncol(clusters), function(col_idx) {
      
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
        comat = comat,
        indices = indices,
        net = net_copy,
        site_col = site_col,
        species_col = species_col,
        is_bipartite = is_bipartite
      )
    })
    
    names(results_list) <- bioregion_names
    
  } else {
    # Single bioregionalization - also return list for consistency
    bioregion_name <- colnames(clusters)[2]
    
    results_list <- list(
      compute_single_bioregionalization_metrics(
        single_clusters = clusters,
        bioregion_name = bioregion_name,
        comat = comat,
        indices = indices,
        net = net,
        site_col = site_col,
        species_col = species_col,
        is_bipartite = is_bipartite
      )
    )
    names(results_list) <- bioregion_name
  }
  
  # Add class and attributes for S3 method dispatch
  class(results_list) <- c("bioregion.site.species.metrics", "list")
  attr(results_list, "indices") <- indices
  attr(results_list, "n_bioregionalizations") <- length(results_list)
  
  return(results_list)
}


# Helper function to compute metrics for a single bioregionalization
#' @noRd
compute_single_bioregionalization_metrics <- function(single_clusters,
                                                      bioregion_name,
                                                      comat,
                                                      indices,
                                                      net = NULL,
                                                      site_col = 1,
                                                      species_col = 2,
                                                      is_bipartite = FALSE) {
  
  rho_df <- affinity_df <- fidelity_df <- indval_df <- bipartite_df <- NULL
  
  # 1. Cz indices (bipartite only) ------------------------------------------
  if("Cz" %in% indices && is_bipartite) {
    # Rename site and species columns as Sites and Species
    colnames(net)[site_col] <- "Site"
    colnames(net)[species_col] <- "Species"
    
    bipartite_df <- single_clusters
    # Add a column category (site or species) to bipartite_df
    bipartite_df$cat <- attributes(single_clusters)$node_type
    colnames(bipartite_df) <- c("Node", "Bioregion", "Category")
    
    # Add bioregions of the sites to the bipartite data.frame
    net$Site <- as.character(net$Site)
    bipartite_df$Node <- as.character(bipartite_df$Node)
    net <- dplyr::left_join(net,
                            bipartite_df[, c("Node", "Bioregion")],
                            by = c("Site" = "Node"))
    colnames(net)[colnames(net) == "Bioregion"] <- "Bioregion_site"
    
    # Add bioregions of the species to the bipartite data.frame
    net$Species <- as.character(net$Species)
    bipartite_df$Node <- as.character(bipartite_df$Node)
    net <- dplyr::left_join(net,
                            bipartite_df[, c("Node", "Bioregion")],
                            by = c("Species" = "Node"))
    colnames(net)[colnames(net) == "Bioregion"] <- "Bioregion_species"
    
    # Compute coefficient of participation C
    C_site <- data.frame()
    dat_com <- bipartite_df[which(bipartite_df$Category == "site"), ]
    for(i in 1:nrow(dat_com)){
      tmp <- table(net[which(net$Site == dat_com[i, "Node"]), "Bioregion_species"])
      C_site <- rbind(C_site,
                      data.frame(Node = dat_com[i, "Node"],
                                 C = 1 - sum((tmp/sum(tmp))^2),
                                 Category = "site"))
    }
    
    C_sp <- data.frame()
    dat_sp <- bipartite_df[which(bipartite_df$Category == "species"), ]
    for(i in 1:nrow(dat_sp)){
      tmp <- table(net[which(net$Species == dat_sp[i, "Node"]), "Bioregion_site"])
      C_sp <- rbind(C_sp,
                    data.frame(Node = dat_sp[i, "Node"],
                               C = 1 - sum((tmp/sum(tmp))^2),
                               Category = "species"))
    }
    C_dat <- rbind(C_site, C_sp)
    
    # Merge results with bipartite_df
    bipartite_df <- dplyr::left_join(bipartite_df, C_dat,
                                     by = c("Node",  "Category"))
    
    # Compute z
    bipartite_df$n_link_bioregion <- NA
    for(i in 1:nrow(bipartite_df)){
      if(bipartite_df[i, "Category"] == "site"){
        tmp <- net[which(net$Site == bipartite_df[i, "Node"]), ]
        bipartite_df[i, "n_link_bioregion"] <-
          nrow(tmp[which(tmp$Bioregion_site == tmp$Bioregion_species), ])
      } else{
        tmp <- net[which(net$Species == bipartite_df[i, "Node"]), ]
        bipartite_df[i, "n_link_bioregion"] <-
          nrow(tmp[which(tmp$Bioregion_site == tmp$Bioregion_species), ])
      }
    }
    
    # Average number of links within a bioregion
    mean_link_bioregion <- tapply(bipartite_df$n_link_bioregion,
                                  bipartite_df$Bioregion,
                                  mean)
    mean_link_bioregion <-
      data.frame(Bioregion = names(mean_link_bioregion),
                 mean_link_bioregion = as.numeric(mean_link_bioregion))
    bipartite_df <- dplyr::left_join(bipartite_df, mean_link_bioregion,
                                     by = "Bioregion")
    
    # Standard deviation of the number of links within a bioregion
    sd_link_bioregion <- tapply(bipartite_df$n_link_bioregion,
                                bipartite_df$Bioregion,
                                stats::sd)
    sd_link_bioregion <-
      data.frame(Bioregion = names(sd_link_bioregion),
                 sd_link_bioregion = as.numeric(sd_link_bioregion))
    bipartite_df <- dplyr::left_join(bipartite_df, sd_link_bioregion,
                                     by = "Bioregion")
    
    # z
    bipartite_df$z <- (bipartite_df$n_link_bioregion -
                         bipartite_df$mean_link_bioregion) /
      bipartite_df$sd_link_bioregion
    
    # Remove intermediate columns
    bipartite_df <- bipartite_df[, c("Node", "Bioregion", "Category", "C", "z")]
  }
  
  # 2. Contribution metrics (rho, affinity, fidelity, indicator_value) ------
  if(any(c("rho", "affinity", "fidelity", "indicator_value") %in% indices)){
    # Binary site-species matrix
    comat_bin <- comat
    comat_bin[comat_bin > 0] <- 1
    
    # If it is a bipartite object, we just consider the sites
    if(is_bipartite){
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
    n_i <- colSums(comat_bin) # number of occurrences per species
    n_species <- ncol(comat)
    
    # Get unique bioregions
    bioregions <- unique(single_clusters[, 2])
    n_bioregions <- length(bioregions)
    n_clust <- n_bioregions
    
    # Loop over bioregions
    for(j in 1:n_bioregions){
      focal_j <- as.character(bioregions[j]) # bioregion j (ensure character for consistent joins)
      
      # Sites belonging to bioregion j
      focal_sites <- single_clusters[which(single_clusters[, 2] == bioregions[j]), 1]
      # Number of sites belonging to bioregion j
      n_j <- sum(single_clusters[, 2] == focal_j)
      
      # Occurrences per species in each of these sites to get n_ij
      n_ij <- colSums(comat_bin[focal_sites, , drop = FALSE])
      
      if("rho" %in% indices){
        # Contribution of species i to bioregion j
        p_ij <- (n_ij - ((n_i*n_j)/n))/(sqrt(((n - n_j)/(n-1))*(1-(n_j/n))*((n_i*n_j)/n)))
        
        rho_df <- rbind(rho_df,
                        data.frame(Bioregion = focal_j,
                                   Species = names(p_ij),
                                   rho = as.numeric(p_ij)))
      }
      if("affinity" %in% indices){
        # Affinity of species i to bioregion j
        affinity_df <- rbind(affinity_df,
                             data.frame(Bioregion = focal_j,
                                        Species = names(n_ij),
                                        affinity = n_ij/n_j))
      }
      
      if("fidelity" %in% indices){
        # Fidelity of species i to bioregion j
        fidelity_df <- rbind(fidelity_df,
                             data.frame(Bioregion = focal_j,
                                        Species = names(n_ij),
                                        fidelity = n_ij/n_i))
      }
      
      if("indicator_value" %in% indices){
        # Individual contribution of species i to bioregion j
        indval_df <- rbind(indval_df,
                           data.frame(Bioregion = focal_j,
                                      Species = names(n_ij),
                                      indval = (n_ij/n_j)*(n_ij/n_i)))
      }
    }
    
    # Merge all outputs together into a single data.frame
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
    if(length(unique(res_df$Bioregion)) != n_clust){
      warning("Not all bioregions are in the output.")
    }
    
    # test if all species are there
    if(length(unique(res_df$Species)) != ncol(comat)){
      warning("Not all species are in the output.")
    }
    
    # test if all species are there X times (X = nb of bioregions)
    if(length(unique(table(res_df$Species))) != 1 ||
       unique(table(res_df$Species)) != n_clust){
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
  }
  
  return(output)
}
