# library(virtualspecies)
#
# a <- matrix(rep(dnorm(1:100, 50, sd = 25)),
#             nrow = 100, ncol = 100, byrow = TRUE)
# env <- stack(raster(a * dnorm(1:100, 50, sd = 25)),
#              raster(a * 1:100))
# names(env) <- c("variable1", "variable2")
#
#
# # Generation of the virtual species
#
#
# centers <- data.frame(variable1 = c(2.5e-4,
#                                     5e-6,
#                                     1.3e-4),
#                       variable2 = c(0.75,
#                                     0,
#                                     0.75))
#
# parameters <-
#   lapply(1:3,
#          function(x)
#            formatFunctions(variable1 = c(fun = 'dnorm',
#                                          mean = centers[x, 1],
#                                          sd = 1e-04),
#                            variable2 = c(fun = 'dnorm',
#                                          mean = centers[x, 2],
#                                          sd = 0.5)))
#
#
# sp <- lapply(1:3,
#              function(x, env, parameters)
#                generateSpFromFun(env, parameters[[x]]),
#              env = env, parameters = parameters)
#
# sp1.lowprev <- convertToPA(sp1, species.prevalence = 0.25, alpha = -.05)
# plotSuitabilityToProba(sp1.lowprev)
#
#
# sp.stack <- stack()
# for (species in 1:200)
# {
#   cent <- sample(1:3, size = 1)
#
#   sp.stack <- stack(sp.stack,
#                     convertToPA(sp[[cent]],
#                                 species.prevalence = sample(seq(0.05, 0.25, length = 100),
#                                                             size = 1),
#                                 alpha = -.05)$pa.raster)
# }
#
# sp.df <- getValues(sp.stack)
# sp.df <- apply(sp.df, 2, as.numeric)
# colnames(sp.df) <- paste0("sp", 1:200)
# rownames(sp.df) <- paste0("cell", 1:nrow(sp.df))
# saveRDS(sp.df, "./data-raw/virtual_community_3clusters.RDS")

sp.df <- readRDS("./data-raw/virtual_community_3clusters.RDS")
similarity_vs <- spproject(sp.df, metric = "Simpson")

dist_vs <- similarity_to_distance(similarity_vs)

dist_vs <- dist_vs[-which(is.na(dist_vs$Simpson)), ]
# Compute the hierarchical tree
vs_hclust <- clustering_hierarchical(dist_vs,
                                     randomize = FALSE)



# Compute the tree and cut it with 12 clusters
vs_hclust <- clustering_hierarchical(dist_vs, n_clust = 3)

plot(inv_hclust)
