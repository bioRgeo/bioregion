#' Calculate contribution metrics of sites and species
#' 
#' This function calculates metrics that assess the contribution of a given
#' species or site to its bioregion.
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
#' @param bipartite_link `NULL` by default. Needed for `Cz` indices. A
#' `data.frame` where each row represents the interaction between two nodes
#'  and an optional third column indicating the weight of the interaction.  
#' 
#' @param indices a `character` specifying the contribution metric to compute.
#' Available options are `rho` and `Cz`.
#' 
#' @details 
#' The \eqn{\rho} metric is derived from \insertRef{Lenormand2019}{bioregion}.
#' Its formula is the following:
#' \eqn{\rho_{ij} = (n_ij - ((n_i n_j)/n))/(sqrt(((n - n_j)/(n-1)) (1-(n_j/n)) ((n_i n_j)/n)))}
#' 
#' with \eqn{n} the number of sites, \eqn{n_i} the number of sites in which
#' species \eqn{i} is present, \eqn{n_j} the number of sites belonging to the
#' bioregion \eqn{j}, \eqn{n_ij} the number of occurrences of species \eqn{i}
#' in sites belonging to the bioregion \eqn{j}.
#' 
#' Affinity, fidelity and individual contributions describe how species are
#' linked to their bioregions. These metrics are described in
#' \insertRef{Bernardo2019}{bioregion}
#' Affinity of species to their region, \eqn{A_i = R_i / Z}
#' where \eqn{R_i} is the occurrence/range size of species \eqn{i} in its
#' associated bioregion, and \eqn{Z} the total size (number of sites) of the
#' bioregion.
#' 
#' A high affinity means that the species is occupying most sites of its
#' associated bioregion.
#' 
#' Fidelity of species to their region, \eqn{F_i = R_i / D_i}
#' where \eqn{R_i} is the occurrence/range size of species \eqn{i} in its
#' associated bioregion, and \eqn{D_i} is its total occurrence/range size.
#' 
#' A high fidelity means that the species is not present in other bioregions
#' than their associated one.
#' 
#' Indicator Value of species, \eqn{IndVal = F_i * A_i}
#' 
#' `Cz` metrics are derived from \insertRef{Guimera2005}{bioregion}.
#' Their respective formula are:
#' \eqn{C_i = 1 - \sum_{s=1}^{N_M}{{(\frac{k_is}{k_i}})^2}}
#' 
#' where \eqn{k_{is}} is the number of links of node (species or site) \eqn{i}
#' to nodes in bioregion \eqn{s}, and \eqn{k_i} is the total degree of node
#' \eqn{i}. The participation coefficient of a node is therefore close to 1 if
#' its links are uniformly distributed among all the bioregions and 0 if all
#' its links are within its own bioregion.
#' 
#' And:
#' \eqn{z_i = \frac{k_i - \overline{k_{si}}}{\sigma_{k_{si}}}}
#' 
#' where \eqn{k_i} is the number of links of node (species or site) \eqn{i} to
#' other nodes in its bioregion \eqn{s_i}, \eqn{\overline{k_{si}}} is the
#' average of \eqn{k} over all the nodes in \eqn{s_i}, and
#' \eqn{\sigma_{k_{si}}} is the standard deviation of \eqn{k} in \eqn{s_i}.
#' The within-bioregion degree z-score measures how well-connected node \eqn{i}
#' is to other nodes in the bioregion.
#'
#' @seealso [bioregionalization_metrics]
#' @return 
#' A `data.frame` with the columns Bioregion, Species, and the desired summary
#' statistics, or a list of `data.frames` if `Cz` and other indices are
#' selected.
#'  
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}),
#' Boris Leroy (\email{leroy.boris@gmail.com}) and
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) 
#' 
#' @examples
#' comat <- matrix(sample(1000, 50), 5, 10)
#' rownames(comat) <- paste0("Site", 1:5)
#' colnames(comat) <- paste0("Species", 1:10)
#' 
#' comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
#'                 20, 25)
#' rownames(comat) <- paste0("Site",1:20)
#' colnames(comat) <- paste0("Species",1:25)
#' 
#' dissim <- dissimilarity(comat, metric = "Simpson")
#' clust1 <- nhclu_kmeans(dissim, n_clust = 3, index = "Simpson")
#' 
#' net <- similarity(comat, metric = "Simpson")
#' com <- netclu_greedy(net)
#' 
#' site_species_metrics(cluster_object = clust1, comat = comat,
#' indices = "rho")
#' 
#' # Contribution metrics
#' site_species_metrics(cluster_object = com, comat = comat,
#' indices = c("rho", "affinity", "fidelity", "indicator_value"))
#' 
#' # Cz indices
#' net_bip <- mat_to_net(comat, weight = TRUE)
#' clust_bip <- netclu_greedy(net_bip, bipartite = TRUE)
#' site_species_metrics(cluster_object = clust_bip, comat = comat, 
#' bipartite_link = net_bip, indices = "Cz")
#' 
#' @export

site_species_metrics <- function(cluster_object,
                                 comat,
                                 indices = c("rho"),
                                 bipartite_link = NULL){
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
        please extract clusters from the tree before using bioregionalization_metrics()
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
  
  # controls(args = indices, data = NULL, type = "character")
  if (!is.character(indices)) {
    stop(paste0(deparse(substitute(indices)),
                " must be a character or a vector of characters."),
         call. = FALSE
    )
  }
  
  if(!isTRUE(unique(indices %in% c("rho", "Cz", "affinity", "fidelity",
                                   "indicator_value")))){
    stop("Please choose algorithm among the followings values:
    rho, affinity, fidelity, indicator_value or Cz.", call. = FALSE)
  }
  
  if("Cz" %in% indices && is.null(bipartite_link)){
    stop("bipartite_link is needed to compute Cz indices.")
  }
  
  if("Cz" %in% indices && cluster_object$inputs$bipartite == FALSE){
    stop("Cz metrics can only be computed for a bipartite partition (where
         both sites and species are assigned to a bioregion.")
  }
  
  # Add controls for bipartite_link
  
  rho_df <- affinity_df <- fidelity_df <- indval_df <- NULL
  
  # 2. Function ---------------------------------------------------------------
  ## 2.1. Cz ------------------------------------------------------------------
  # only for bipartite cases; not implemented yet
  if("Cz" %in% indices){
    # cluster_object <- clust2 # remove
    # bipartite_link <- net_bip # remove
    
    # Needs two data frames as inputs:
    # bipartite_link is a data.frame of the links between nodes containing four
    # columns: site, species, bioregion_site and bioregion_species 
    # 
    # bipartite_df contains three columns: node, bioregion and cat which
    # respectively stand for the name of the node, its bioregion and its
    # bipartite category
    
    # Rename columns Node1 as Sites and Node2 as Species, and throw warning if
    # this is the case
    if("Node1" %in% colnames(bipartite_link)){
      colnames(bipartite_link)[colnames(bipartite_link) == "Node1"] <- "Site"
      warning("Column 'Node1' has been renamed 'Sites'. If this column does not
            correspond to sites, rename it before running the function.")
    }
    if("Node2" %in% colnames(bipartite_link)){
      colnames(bipartite_link)[colnames(bipartite_link) == "Node2"] <- "Species"
      warning("Column 'Node2' has been renamed 'Species'. If this column does not
            correspond to species, rename it before running the function.")
    }
    
    bipartite_df <- cluster_object$clusters
    # Add a column category (site or species) to bipartite_df
    bipartite_df$cat <- attributes(bipartite_df)$node_type
    colnames(bipartite_df) <- c("Node", "Bioregion", "Category")
    
    # Add bioregions of the sites to the bipartite data.frame
    bipartite_link$Site <- as.character(bipartite_link$Site)
    bipartite_df$Node <- as.character(bipartite_df$Node)
    bipartite_link <- dplyr::left_join(bipartite_link,
                                       bipartite_df[, c("Node", "Bioregion")],
                                       by = c("Site" = "Node"))
    colnames(bipartite_link)[colnames(bipartite_link) == "Bioregion"] <-
      "Bioregion_site"
    
    # Add bioregions of the species to the bipartite data.frame
    bipartite_link$Species <- as.character(bipartite_link$Species)
    bipartite_df$Node <- as.character(bipartite_df$Node)
    bipartite_link <- dplyr::left_join(bipartite_link,
                                       bipartite_df[, c("Node", "Bioregion")],
                                       by = c("Species" = "Node"))
    colnames(bipartite_link)[colnames(bipartite_link) == "Bioregion"] <-
      "Bioregion_species"
    
    
    # Compute coefficient of participation C
    C_site <- data.frame()
    dat_com <- bipartite_df[which(bipartite_df$Category == "site"), ]
    for(i in 1:nrow(dat_com)){
      tmp <-
        table(bipartite_link[which(bipartite_link$Site == dat_com[i, "Node"]),
                             "Bioregion_species"])
      C_site <- rbind(C_site,
                      data.frame(Node = dat_com[i, "Node"],
                                 C = 1 - sum((tmp/sum(tmp))^2),
                                 Category = "site"))
    }
    
    C_sp <- data.frame()
    dat_sp <- bipartite_df[which(bipartite_df$Category == "species"), ]
    for(i in 1:nrow(dat_sp)){
      tmp <-
        table(bipartite_link[which(bipartite_link$Species == dat_sp[i, "Node"]),
                             "Bioregion_site"])
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
        tmp <- bipartite_link[which(bipartite_link$Site == 
                                      bipartite_df[i, "Node"]), ]
        bipartite_df[i, "n_link_bioregion"] <-
          nrow(tmp[which(tmp$Bioregion_site == tmp$Bioregion_species), ])
      } else{
        tmp <- bipartite_link[which(bipartite_link$Species ==
                                      bipartite_df[i, "Node"]), ]
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
    bipartite_df <-
      bipartite_df[, c("Node", "Bioregion", "Category", "C", "z")]
  }
  
  ## 2.2. Metrics -------------------------------------------------------------
  if(any(c("rho", "affinity", "fidelity", "indicator_value") %in% indices)){
    # Binary site-species matrix
    comat_bin <- comat
    comat_bin[comat_bin > 0] <- 1
    
    # If it is a bipartite object, we just consider the sites
    if(cluster_object$inputs$bipartite == TRUE){
      cluster_object$clusters <-
        cluster_object$clusters[
          which(attributes(cluster_object$clusters)$node_type == "site"), ]
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
    fidelity_df <- data.frame(Bioregion = character(),
                              Species = character(),
                              fidelity = character())
    indval_df <- data.frame(Bioregion = character(),
                            Species = character(),
                            indval = character())
    
    # Formula
    n <- nrow(comat) # number of sites
    n_i <- colSums(comat_bin) # number of occurrences per species
    n_j <- table(cluster_object$clusters[, 2]) # number of sites per bioregion
    
    # Loop over bioregions
    for(j in 1:cluster_object$cluster_info$n_clust){
      focal_j <- unique(cluster_object$clusters[, 2])[j] # bioregion j
      
      # Sites belonging to bioregion j
      focal_sites <- cluster_object$clusters[which(
        cluster_object$clusters[, 2] == focal_j), 1]
      # Number of sites belonging to bioregion j
      n_j <- table(cluster_object$clusters[, 2])[[focal_j]]
      
      # Occurrences per species in each of these sites to get n_ij
      n_ij <- colSums(comat[focal_sites, , drop = FALSE])
      
      if("rho" %in% indices){
        # Contribution of species i to bioregion j
        p_ij <-
          (n_ij - ((n_i*n_j)/n))/(sqrt(((n - n_j)/
                                          (n-1))*(1-(n_j/n))*((n_i*n_j)/n)))
        
        rho_df <- rbind(rho_df,
                        data.frame(Bioregion = focal_j,
                                   Species = names(p_ij),
                                   rho = as.numeric(p_ij)))
      }
      if("affinity" %in% indices){
        # Affinity of species i to bioregion j
        affinity_df <- rbind(affinity_df,
                             data.frame(Bioregion = focal_j,
                                        Species = names(p_ij),
                                        affinity = n_ij/n_j))
      }
      
      if("fidelity" %in% indices){
        # Fidelity of species i to bioregion j
        fidelity_df <- rbind(fidelity_df,
                             data.frame(Bioregion = focal_j,
                                        Species = names(p_ij),
                                        fidelity = n_ij/n_i))
      }
      
      if("indicator_value" %in% indices){
        # Individual contribution of species i to bioregion j
        indval_df <- rbind(indval_df,
                           data.frame(Bioregion = focal_j,
                                      Species = names(p_ij),
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
    if(length(unique(res_df$Bioregion)) !=
       cluster_object$cluster_info$n_clust){
      warning("Not all bioregions are in the output.")
    }
    
    # test if all species are there
    if(length(unique(res_df$Species)) != ncol(comat)){
      warning("Not all species are in the output.")
    }
    
    # test if all species are there X times (X = nb of bioregions)
    if(length(unique(table(res_df$Species))) != 1 ||
       unique(table(res_df$Species)) != cluster_object$cluster_info$n_clust){
      warning("Not all species x bioregions combinations are in the output.")
    }
  }
  
  if("Cz" %in% indices){
    if(length(indices) == 1){
      return(bipartite_df)
    }else{
      return(list(res = res_df,
                  cz = bipartite_df))
    }
  }else{
    return(res_df)    
  }
}
