# clusterMetrics <- function(cluster_object,
#                            net = NULL,
#                            site_col = 1,
#                            species_col = 2,
#                            site.area = NULL, # A data.frame containing at least two columns with site name ("name") and site area ("area")
#                            partition = "all"
#                            
#                            
#                            network = NULL,
#                            site.field = colnames(db)[1],
#                            species.field = colnames(db)[2],
#                            level = "lvl1",
#                            site.metrics = TRUE,
#                            species.metrics = TRUE,
#                            stats.per.cluster = FALSE,
#                            ...)
# {
#   tree1$clusters
# 
#   data(fishmat)
#   dissim <- dissimilarity(fishmat)
#   tree1 <- hclu_hierarclust(dissim,
#                             n_clust = c(4, 8))
#                     
#   
#   cluster_object <- tree1
#   clusters <- cluster_object$clusters
#   
#   net <- mat_to_net(fishmat)
#   site_col = "Node1"
#   species_col = "Node2"
#   
#   # CHECK HERE that all nodes from cluster_object are in net
#   
#   # If all nodes from cluster_object are in net, then the next line is ok
#   nodetype <- data.frame(node_id = clusters$ID,
#                          node_type = ifelse(clusters$ID %in% net[, site_col],
#                                             "site",
#                                             "species"))
#   
#   clusters$nodetype <- NA
#   net$nodetype[network$Name %in% net[, site_col]] <- "site"
#   net$nodetype[network$Name %in% net[, species_col]] <- "species"
#   
#   
#   
#   
#   
#   region.stats <- data.frame(cluster = unique(network[, level]),
#                              nb.sites = sapply(unique(network[, level]),
#                                                function(x, net, database)
#                                                {
#                                                  length(unique(database[which(database[, site.field] %in%
#                                                                                 net$Name[which(net[, level] == x)]), site.field]))
#                                                },
#                                                net = network,
#                                                database = db))
#   
#   
#   net <- data.frame(
#     net, 
#     cluster_object$clusters[data.table::chmatch(net[, site_col],
#                                                 cluster_object$clusters$ID),
#                             -1])
#   
#   # Visible binding for global variable
#   N <- endemism <- end_richness <- pc_endemism <- NULL 
#   
#   
#   endemism_results <- lapply(
#     cluster_object$cluster_info$partition_name,
#     function(x, network., species_col.) {
#       # Create species per cluster network in data.table format for faster
#       # calculations
#       species_cluster <- data.table::as.data.table(
#         unique(network.[, c(species_col., x)]))
#       # Calculate richness per cluster
#       rich_clusters <- species_cluster[, .N, by = x]
#       # Calculate species occurrence across clusters
#       occ_sp <- species_cluster[, .N, by = species_col.]
#       # Add new column with endemism status
#       occ_sp[, "endemism" := N == 1] 
#       # Then add endemism status to species_cluster table
#       species_cluster$endemism <-
#         occ_sp$endemism[data.table::chmatch(species_cluster[[species_col.]],
#                                             occ_sp[[species_col.]])] 
#       # Calculate richness of endemics per cluster
#       end_clusters <- species_cluster[endemism == TRUE, .N, by = x]
#       # Merge total & endemism richness tables
#       rich_clusters[, end_richness := end_clusters[
#         data.table::chmatch(rich_clusters[[x]],
#                             end_clusters[[x]]), N]]
#       # Replace NAs (i.e. no endemics) by zeros
#       rich_clusters[is.na(rich_clusters)] <- 0 
#       # Calculate percentage of endemism
#       rich_clusters[, pc_endemism := end_richness/N] 
#       return(rich_clusters)
#     }, network. = net, species_col. = species_col
#   )
#   names(endemism_results) <- cluster_object$cluster_info$partition_name
#   
#   
#   # if(!("nodetype" %in% colnames(network)))
#   # {
#   #   network$nodetype <- NA
#   #     }
#   # 
#   # network
#   # 
#   
#   # 
#   # dots <- list(...)
#   # cat("Network data.frame: ", deparse(substitute(network)), "\n")
#   # cat("Database data.frame: ", deparse(substitute(db)), "\n")
#   # 
# # 
# #   
# #   max.lvl <- 0
# #   loop <- TRUE
# #   while(loop)
# #   {
# #     max.lvl <- max.lvl + 1
# #     res <- try(network[[paste0("lvl", max.lvl)]])
# #     if(is.null(res))
# #     {
# #       loop <- FALSE
# #       max.lvl <- max.lvl - 1
# #     }
# #   }
# #   
# 
#   
# 
#   
#   if(!is.null(site.area))
#   {
#     region.stats <- data.frame(region.stats,
#                                area = sapply(unique(network[, level]),
#                                              function(x, net, surf){
#                                                sum(surf$area[which(surf$name %in%
#                                                                      net$Name[which(net[, level] %in% x)])])
#                                              },
#                                              net = network,
#                                              surf = site.area))
#   }
#   
#   
#   db[, "sp.cluster"] <- network[match(db[, species.field],
#                                       network$Name), level]
#   db[, "site.cluster"] <- network[match(db[, site.field],
#                                         network$Name), level]
#   
#   
#   if("manual.site.correction" %in% names(dots))
#   {
#     db[, "site.cluster"] <- dots$manual.site.correction[match(db[, site.field],
#                                                               dots$manual.site.correction[, 1]), 2]
#   }
#   
#   # db[which(is.na(db[, "site.cluster"])), ]
#   
#   # Contigency table for clusters
#   cluster.contingency <- as.data.frame.matrix(table(db[, species.field],
#                                                     db[, "site.cluster"]))
#   cluster.contingency[cluster.contingency > 1] <- 1
#   
#   # Occurrence in clusters
#   cluster.occurrences <- rowSums(cluster.contingency)
#   
#   res$network.stats <- data.frame(
#     nb.nodes = nrow(network),
#     nb.species.network = length(which(network$nodetype == "species")),
#     nb.sites.network = length(which(network$nodetype == "site")),
#     nb.species.db = length(unique(db[, species.field])),
#     nb.sites.db = length(unique(db[, site.field])),
#     nb.links = nrow(db),
#     nb.endemics = length(which(cluster.occurrences == 1)),
#     max.level = max.lvl
#   )
#   
#   res$nb.regions.per.sp <- cluster.occurrences
#   
#   cat("  - ",
#       scales::comma(res$network.stats$nb.species.network),
#       " species and ",
#       scales::comma(res$network.stats$nb.sites.network), " sites in the network\n")
#   cat("  - ",
#       scales::comma(res$network.stats$nb.species.db),
#       " species and ",
#       scales::comma(res$network.stats$nb.sites.db), " sites in the database\n")
#   cat("  - ",
#       scales::comma(res$network.stats$nb.links),
#       " links in the database\n")
#   cat("  - Maximum hierarchical level: ",
#       res$network.stats$max.level, "\n")
#   if(any(is.na(db[, "site.cluster"])))
#   {
#     cat("  - Note that ",
#         length(unique(db[which(is.na(db[, "site.cluster"])), site.field])),
#         " sites do not have a cluster, i.e. they exist in db",
#         " but do not exist in network",
#         ". Some metrics will not be calculated for these sites.", sep = "")
#   }
#   if(any(is.na(db[, "sp.cluster"])))
#   {
#     cat("  - Note that ",
#         length(unique(db[which(is.na(db[, "sp.cluster"])), species.field])),
#         " species do not have a cluster, i.e. they exist in db",
#         " but do not exist in network",
#         ". Some metrics will not be calculated for these species.\n", sep = "")
#   }
#   
#   cat("Cluster level under evaluation: ", level, "\n")
#   cat("  - Number of endemics: ",
#       scales::comma(res$network.stats$nb.endemics),
#       " species\n")
#   
#   
#   region.stats[, c("richness", "char.richness", "end.richness", "nested.levels")] <- NA
#   for (i in 1:nrow(region.stats))
#   {
#     reg <- region.stats$cluster[i]
#     
#     # Characteristic sites of current cluster
#     characteristic.sites <- network$Name[which(network[, level] %in% reg &
#                                                  network$nodetype == "site")]
#     
#     # Characteristic species of current region
#     characteristic.species <- network$Name[which(network[, level] %in% reg &
#                                                    network$nodetype == "species")]
#     region.stats$char.richness[i] <- length(characteristic.species)
#     
#     # Contingency matrix site x  species of current region
#     cur.reg.contin <- as.data.frame.matrix(table(db[which(db[, site.field] %in% characteristic.sites), species.field],
#                                                  db[which(db[, site.field] %in% characteristic.sites), site.field]))
#     
#     # Occurrence of species in current region
#     reg.occ <- rowSums(cur.reg.contin)
#     reg.occ[reg.occ > 1] <- 1
#     
#     # Names species of current region
#     sp.in.cur.reg <- names(reg.occ)[reg.occ > 0]
#     
#     # Species richness of current region
#     region.stats$richness[i] <- length(sp.in.cur.reg)
#     
#     # Total cluster occurrence of species in current region
#     cluster.occurrences.cur.reg <- cluster.occurrences[
#       names(cluster.occurrences) %in% sp.in.cur.reg]
#     region.stats$end.richness[i] <- length(which(cluster.occurrences.cur.reg == 1))
#     
#     
#     # Counting the number of hierarchical levels
#     tmp <- network[which(network[[level]] %in% reg), ]
#     
#     
#     region.stats$nested.levels[i] <-
#       max(which(!is.na(tmp[, paste0("lvl", 1:max.lvl)]), arr.ind = T)[, 2])
#   }
#   res$region.stats <- region.stats
#   
#   
#   
#   if(species.metrics) {
#     
#     cat("Computing species indices...\n")
#     sp.stats <- data.frame(species = unique(db[, species.field]))
#     
#     # Cluster of our species
#     sp.stats[, "cluster"] <- network[match(sp.stats$species,
#                                            network$Name), level]
#     
#     
#     sp.stats[, "Endemism"] <- (res$nb.regions.per.sp == 1)[match(sp.stats$species,
#                                                                  names(res$nb.regions.per.sp))]
#     rownames(sp.stats) <- sp.stats$species
#     for (sp in sp.stats$species)
#     {
#       # All sites where species occurs
#       site.occurrence.total <- unique(db[which(db[, species.field] == sp), site.field])
#       # Sites of the native region where the species occurs
#       site.occurrence.region <- site.occurrence.total[which(site.occurrence.total %in%
#                                                               network$Name[which(network[, level] ==
#                                                                                    sp.stats[sp, "cluster"])])]
#       
#       
#       
#       # Raw occurrence of the species in the region to which it was assigned
#       sp.stats[sp, "Occ.Ri"] <- length(site.occurrence.region)
#       # Total raw occurrence of the species
#       sp.stats[sp, "Occ.Di"] <- length(site.occurrence.total)
#       
#       
#       # Occurrence-based affinity
#       sp.stats[sp, "Occ.Ai"] <- sp.stats[sp, "Occ.Ri"] / region.stats$nb.sites[
#         which(region.stats$cluster == sp.stats[sp, "cluster"])]
#       # Occurrence-based Fidelity
#       sp.stats[sp, "Occ.Fi"] <- sp.stats[sp, "Occ.Ri"] / sp.stats[sp, "Occ.Di"]
#       
#       # Occurrence-based IndVal
#       sp.stats[sp, "Occ.IndVal"] <- sp.stats[sp, "Occ.Ai"] * sp.stats[sp, "Occ.Fi"]
#       # Occurrence-based DilVal
#       sp.stats[sp, "Occ.DilVal"] <- sp.stats[sp, "Occ.Ai"] * (1 - sp.stats[sp, "Occ.Fi"])
#       
#       if(!is.null(site.area))
#       {
#         # area of the distribution range of the species in the region to which it was assigned
#         sp.stats[sp, "Ri"] <- sum(site.area$area[which(site.area$name %in% site.occurrence.region)])
#         # Total area the distribution range of the species
#         sp.stats[sp, "Di"] <- sum(site.area$area[which(site.area$name %in% site.occurrence.total)])
#         
#         # area-based affinity
#         sp.stats[sp, "Ai"] <- sp.stats[sp, "Ri"] / region.stats$area[
#           which(region.stats$cluster == sp.stats[sp, "cluster"])]
#         # area-based Fidelity
#         sp.stats[sp, "Fi"] <- sp.stats[sp, "Ri"] / sp.stats[sp, "Di"]
#         
#         # area-based IndVal
#         sp.stats[sp, "IndVal"] <- sp.stats[sp, "Ai"] * sp.stats[sp, "Fi"]
#         # area-based DilVal
#         sp.stats[sp, "DilVal"] <- sp.stats[sp, "Ai"] * (1 - sp.stats[sp, "Fi"])
#       }
#     }
#     
#     res$species.stats <- sp.stats
#   }
#   
#   if(site.metrics) {
#     
#     cat("Computing site indices...\n")
#     site.stats <- data.frame(site = unique(db[, site.field]))
#     rownames(site.stats) <- site.stats$site
#     
#     # Cluster of current site
#     site.stats$cluster <- network[match(site.stats$site, network$Name), level]
#     
#     # Site contingency table
#     tot.contin <- as.data.frame.matrix(table(db[, species.field],
#                                              db[, site.field]))
#     tot.contin[tot.contin > 1] <- 1
#     
#     # Site richness
#     rich <- colSums(tot.contin)
#     site.stats$richness <- rich[match(site.stats$site, names(rich))]
#     
#     if(stats.per.cluster)
#     {
#       site.stats.per.cluster <- data.frame(site = unique(db[, site.field]))
#       rownames(site.stats.per.cluster) <- site.stats.per.cluster$site
#       
#       # Cluster of current site
#       site.stats.per.cluster$cluster <- network[match(site.stats.per.cluster$site, network$Name), level]
#     }
#     
#     # unknown.sites <- NULL
#     
#     cur_col <- 1
#     for (site in site.stats$site)
#     {
#       sp.in.site <- unique(db[which(db[, site.field] == site), species.field])
#       
#       # Species endemic to the current cluster
#       endemic.sp <- sp.in.site[which(sp.in.site %in% names(res$nb.regions.per.sp)[res$nb.regions.per.sp == 1])]
#       site.stats[site, "end.richness"]<- length(endemic.sp)
#       
#       if(site %in% network$Name) # In case new sites are investigated but were not classified in the initial clusters
#       {
#         
#         # Characteristic species
#         characteristic.sp <- sp.in.site[which(sp.in.site %in%
#                                                 network$Name[which(network[, level] == site.stats[site, "cluster"])])]
#         
#         site.stats[site, "char.richness"]<- length(characteristic.sp)
#         
#         
#         # Non characteristic species
#         noncharacteristic.sp <- sp.in.site[which(!(sp.in.site %in%
#                                                      network$Name[which(network[, level] == site.stats[site, "cluster"])]))]
#         
#         
#         # Occurrence-based robustness
#         if(species.metrics) {
#           site.stats[site, "Occ.Rg"] <- sum(sp.stats$Occ.IndVal[which(sp.stats$species %in% characteristic.sp)]) -
#             sum(sp.stats$Occ.DilVal[which(sp.stats$species %in% noncharacteristic.sp)])
#           # Occurrence-based relative robustness (= Rg/ species richness)
#           site.stats[site, "Occ.RRg"] <- site.stats[site, "Occ.Rg"] / length(sp.in.site)
#         } else {
#           warning("Site robustness cannot be computed without species metrics")
#         }
#         
#         
#         if(!is.null(site.area))
#         {
#           # area-based robustness
#           if(species.metrics) {
#             site.stats[site, "Rg"] <- sum(sp.stats$IndVal[which(sp.stats$species %in% characteristic.sp)]) -
#               sum(sp.stats$DilVal[which(sp.stats$species %in% noncharacteristic.sp)])
#             # area-based relative robustness (= Rg/ species richness)
#             site.stats[site, "RRg"] <- site.stats[site, "Rg"] / length(sp.in.site)
#           }
#         }
#       } # else
#       # {
#       #   unknown.sites <- c(unknown.sites,
#       #                      site)
#       # }
#       if(stats.per.cluster)
#       {
#         for(focal.cluster in region.stats$cluster)
#         {
#           site.cluster <- site.stats[site, "cluster"]
#           #                focal cluster
#           #                  yes     no
#           # site      yes    A       B
#           # cluster   no     C       D
#           
#           # Characteristic species
#           A <- sp.in.site[which(sp.in.site %in%
#                                   network$Name[which(network[, level] == site.cluster)] &
#                                   sp.in.site %in%
#                                   network$Name[which(network[, level] == focal.cluster)])]
#           
#           B <- sp.in.site[which(sp.in.site %in%
#                                   network$Name[which(network[, level] == site.cluster)] &
#                                   !(sp.in.site %in%
#                                       network$Name[which(network[, level] == focal.cluster)]))]
#           
#           C <- sp.in.site[which(sp.in.site %in%
#                                   network$Name[which(network[, level] == focal.cluster)] &
#                                   !(sp.in.site %in%
#                                       network$Name[which(network[, level] == site.cluster)]))]
#           
#           D <- sp.in.site[which(!(sp.in.site %in%
#                                     network$Name[which(network[, level] == focal.cluster)]) &
#                                   !(sp.in.site %in%
#                                       network$Name[which(network[, level] == site.cluster)]))]
#           
#           
#           # Occurrence-based robustness
#           if(site.cluster == focal.cluster)
#           {
#             site.stats.per.cluster[site, paste0("Occ.Rg.", focal.cluster)] <- sum(sp.stats$Occ.IndVal[which(sp.stats$species %in% A)]) -
#               sum(sp.stats$Occ.DilVal[which(sp.stats$species %in% D)])
#             # Occurrence-based relative robustness (= Rg/ species richness)
#             site.stats.per.cluster[site, paste0("Occ.RRg.", focal.cluster)] <- site.stats.per.cluster[site, paste0("Occ.Rg.", focal.cluster)] / length(sp.in.site)
#             
#             if(!is.null(site.area))
#             {
#               # area-based robustness
#               site.stats.per.cluster[site, paste0("Rg.", focal.cluster)] <- sum(sp.stats$IndVal[which(sp.stats$species %in% A)]) -
#                 sum(sp.stats$DilVal[which(sp.stats$species %in% D)])
#               # area-based relative robustness (= Rg/ species richness)
#               site.stats.per.cluster[site, paste0("RRg.", focal.cluster)] <- site.stats.per.cluster[site, paste0("Rg.", focal.cluster)] / length(sp.in.site)
#             }
#           } else
#           {
#             site.stats.per.cluster[site, paste0("Occ.Rg.", focal.cluster)] <- sum(sp.stats$Occ.IndVal[which(sp.stats$species %in% B)]) +
#               sum(sp.stats$Occ.DilVal[which(sp.stats$species %in% D)]) -
#               sum(sp.stats$Occ.DilVal[which(sp.stats$species %in% C)])
#             # Occurrence-based relative robustness (= Rg/ species richness)
#             site.stats.per.cluster[site, paste0("Occ.RRg.", focal.cluster)] <- site.stats.per.cluster[site, paste0("Occ.Rg.", focal.cluster)] / length(sp.in.site)
#             
#             if(!is.null(site.area))
#             {
#               # area-based robustness
#               site.stats.per.cluster[site, paste0("Rg.", focal.cluster)] <- 1 - sum(sp.stats$IndVal[which(sp.stats$species %in% B)]) +
#                 sum(sp.stats$DilVal[which(sp.stats$species %in% D)]) -
#                 sum(sp.stats$DilVal[which(sp.stats$species %in% C)])
#               # area-based relative robustness (= Rg/ species richness)
#               site.stats.per.cluster[site, paste0("RRg.", focal.cluster)] <- site.stats.per.cluster[site, paste0("Rg.", focal.cluster)] / length(sp.in.site)
#             }
#           }
#         }
#         cat(".")
#         if((cur_col) == length(site.stats$site))
#         {
#           cat("Complete\n")
#         } else if((cur_col) %% 50 == 0)
#         {
#           cat(paste0("Site N°", cur_col, " ~ ", round(100 * (cur_col) / length(site.stats$site), 2), "% complete\n"))
#         }
#         
#         cur_col <- cur_col + 1
#       }
#     }
#     res$site.stats <- site.stats 
#   }
#   # if(length(unknown.sites))
#   # {
#   #   warning(paste0(length(unknown.sites), " sites do not exist in the cluster table.
#   #   Site indices will not be calculated for these sites: \n",
#   #                  paste(unknown.sites, collapse = "\n")))
#   # }
#   if(stats.per.cluster)
#   {
#     res$site.stats.per.cluster <- site.stats.per.cluster
#   }
#   return(res)
# }