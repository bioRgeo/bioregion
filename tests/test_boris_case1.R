library(bioRgeo)

contin <- t(as.matrix(readRDS("./data-raw/benthic_invertebrates.RDS")))
site_loc <- readRDS("./data-raw/benthic_invertebrates_sites.RDS")

#### Distance-based hierarchical clustering
# Calculate similarity
similarity_inv <- spproject(contin)
# Convert to distance
dist_inv <- similarity_to_distance(similarity_inv)

# Compute the hierarchical tree
inv_hclust <- clustering_hierarchical(dist_inv)

# Compute the tree and cut it with 12 clusters
inv_hclust <- clustering_hierarchical(dist_inv, n_clust = 12)

plot(inv_hclust)

# Cut the tree with dynamic cutting
inv_hclust_cutdyn <- cut_tree(inv_hclust,
                            dynamic_tree_cut = TRUE)
plot(inv_hclust_cut1)

# Search for optimal number of clusters
optim_n_1 <- find_nclust_tree(inv_hclust)
optim_n_2 <- find_nclust_tree(inv_hclust, eval_metric = "anosim")
optim_n_3 <- find_nclust_tree(inv_hclust, eval_metric = "avg_endemism",
                            sp_site_table = contin,
                            criterion = "elbow")
optim_n_4 <- find_nclust_tree(inv_hclust, eval_metric = "tot_endemism",
                            sp_site_table = contin,
                            criterion = "elbow")
optim_n_5 <- find_nclust_tree(inv_hclust, eval_metric = "pc_distance",
                              criterion = "elbow")


inv_hclust_cut1 <- cut_tree(inv_hclust,
                            n_clust = optim_n_1)
plot(inv_hclust_cut1)
inv_hclust_cut5 <- cut_tree(inv_hclust,
                            n_clust = optim_n_5)
plot(inv_hclust_cut5)

inv_hclust_cut1 <- cut_tree(inv_hclust,
                            n_clust = c(14, 18, 22))
plot(inv_hclust_cut1)
