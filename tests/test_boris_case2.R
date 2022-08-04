download.file("https://www.biorxiv.org/content/biorxiv/early/2019/06/12/319566/DC1/embed/media-1.csv",
              destfile = "data-raw/fish.csv", mode = "wb")
download.file("https://borisleroy.com/permanent/basinshapefile.RDS",
              destfile = "data-raw/basinshapefile.RDS", mode = "wb")


library(rgdal)
library(raster)
library(bioRgeo)
library(sf)
library(ggplot2)
install_binaries()

fish <- read.csv("data-raw/fish.csv")
basins <- readRDS("data-raw/basinshapefile.RDS")
basins <- st_as_sf(basins)

fish_comat <- net_to_mat(fish)

fish_dist <- distance(fish_comat)


##### INFOMAP #####
cl_infomap <- netclu_infomap(fish,
                             weight = FALSE,
                             numtrials = 10,
                             twolevel = FALSE,
                             direct = FALSE,
                             bipartite = TRUE,
                             primary_col = "X1.Basin.Name",
                             feature_col = "X6.Fishbase.Valid.Species.Name",
                             remove_feature = FALSE)

basins$cl_infomap <- cl_infomap[match(basins$BasinName,
                                      cl_infomap$ID), 3]
basins$cl_infomap <- as.factor(basins$cl_infomap)

ggplot() +
  geom_sf(data = basins,
          aes(fill = cl_infomap)) +
  scale_fill_discrete()

##### Beckett #####
cl_beckett <- netclu_beckett(fish,
                             weight = FALSE,
                             primary_col = "X1.Basin.Name",
                             feature_col = "X6.Fishbase.Valid.Species.Name",
                             remove_feature = FALSE)

basins$cl_beckett <- cl_beckett[match(basins$BasinName,
                                      cl_beckett$ID), 2]
basins$cl_beckett <- as.factor(basins$cl_beckett)

ggplot() +
  geom_sf(data = basins,
          aes(fill = cl_infomap)) +
  scale_fill_discrete()


##### greedy #####
cl_greedy <- netclu_greedy(fish,
                           weight = FALSE,
                           bipartite = TRUE,
                           primary_col = "X1.Basin.Name",
                           feature_col = "X6.Fishbase.Valid.Species.Name",
                           remove_feature = FALSE)

basins$cl_greedy <- cl_greedy[match(basins$BasinName,
                                      cl_greedy$ID), 2]
basins$cl_greedy <- as.factor(basins$cl_greedy)

ggplot() +
  geom_sf(data = basins,
          aes(fill = cl_infomap)) +
  scale_fill_discrete()

##### labelprop #####
cl_labelprop <- netclu_labelprop(fish,
                             weight = FALSE,
                             bipartite = TRUE,
                             primary_col = "X1.Basin.Name",
                             feature_col = "X6.Fishbase.Valid.Species.Name",
                             remove_feature = FALSE)

basins$cl_labelprop <- cl_labelprop[match(basins$BasinName,
                                      cl_labelprop$ID), 2]
basins$cl_labelprop <- as.factor(basins$cl_labelprop)

ggplot() +
  geom_sf(data = basins,
          aes(fill = cl_infomap)) +
  scale_fill_discrete()

##### leadingeigen #####
cl_leadingeigen <- netclu_leadingeigen(fish,
                             weight = FALSE,
                             bipartite = TRUE,
                             primary_col = "X1.Basin.Name",
                             feature_col = "X6.Fishbase.Valid.Species.Name",
                             remove_feature = FALSE)

basins$cl_leadingeigen <- cl_leadingeigen[match(basins$BasinName,
                                      cl_leadingeigen$ID), 2]
basins$cl_leadingeigen <- as.factor(basins$cl_leadingeigen)

ggplot() +
  geom_sf(data = basins,
          aes(fill = cl_infomap)) +
  scale_fill_discrete()

##### louvain #####
cl_louvain <- netclu_louvain(fish,
                             weight = FALSE,
                             bipartite = TRUE,
                             primary_col = "X1.Basin.Name",
                             feature_col = "X6.Fishbase.Valid.Species.Name",
                             remove_feature = FALSE)

basins$cl_louvain <- cl_louvain[match(basins$BasinName,
                                      cl_louvain$ID), 2]
basins$cl_louvain <- as.factor(basins$cl_louvain)

ggplot() +
  geom_sf(data = basins,
          aes(fill = cl_infomap)) +
  scale_fill_discrete()

##### oslom #####
cl_oslom <- netclu_oslom(fish,
                             weight = FALSE,
                         bipartite = TRUE,
                             primary_col = "X1.Basin.Name",
                             feature_col = "X6.Fishbase.Valid.Species.Name",
                             remove_feature = FALSE)

basins$cl_oslom <- cl_oslom[match(basins$BasinName,
                                      cl_oslom$ID), 2]
basins$cl_oslom <- as.factor(basins$cl_oslom)

ggplot() +
  geom_sf(data = basins,
          aes(fill = cl_infomap)) +
  scale_fill_discrete()

##### walktrap #####
cl_walktrap <- netclu_walktrap(fish,
                             weight = FALSE,
                             bipartite = TRUE,
                             primary_col = "X1.Basin.Name",
                             feature_col = "X6.Fishbase.Valid.Species.Name",
                             remove_feature = FALSE)

basins$cl_walktrap <- cl_walktrap[match(basins$BasinName,
                                      cl_walktrap$ID), 2]
basins$cl_walktrap <- as.factor(basins$cl_walktrap)

ggplot() +
  geom_sf(data = basins,
          aes(fill = cl_infomap)) +
  scale_fill_discrete()

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
