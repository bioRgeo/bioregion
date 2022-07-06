# devtools::source_url("https://raw.githubusercontent.com/Farewe/Cours_Bioregionalisation/master/scripts/modified_recluster_tree_function.R")
#
# library(bioRgeo)
# comat <- matrix(sample(0:1000, size = 50, replace = TRUE, prob = 1/1:1001), 5, 10)
# rownames(comat) <- paste0("Site",1:5)
# colnames(comat) <- paste0("Species",1:10)
#
# data("vegemat")
#
# simil <- spproject(vegemat, metric = "all")
#
# distances <- similarity_to_distance(simil)
#
# require(cluster)
#
# simHclust
#
#
# clusterNonHierarch <- function()
# {
#   if(!(all(method %in% c("kmeans", "meanshift", "dbscan",
#                          "gmm", "diana", "pam")))){
#     stop("Hierarchical clustering method chosen is not available.
#      Please chose among the followings:
#          'kmeans', 'meanshift',
#          'ward.D', 'ward.D2', 'single', 'complete', 'average',
#          'mcquitty', 'median', 'centroid', 'dbscan', 'gmm', 'diana' or
#          'pam'.")
#   }
#
#
#
#   ## 3. Clustering ----
#   # Initiate empty list
#   cluster_list <- list()
#   # for-loop over the selected algorithms
#   for(i in 1:length(method)){
#
#     if(method[i] == "meanshift"){
#       require(meanShiftR)
#       tmp <- meanShift(queryData = dat, algorithm = "LINEAR", alpha = 0,
#                        iterations = 100)
#
#       # Data.frame of results
#       res <- data.frame(site = rownames(dat),
#                         cluster = tmp$assignment)
#
#     } else
#   }
# }
#
# # dist.df <- distances
# method = "ward.D2"
# optim_method = "firstSEmax"
# n_clust = NULL
# nstart = 25
# B = 50
# K.max = 20
#
#
# dat <- .dfToDist(dist.df, metric = "Jaccard")
# dist1 <- as.matrix(a)
#
#
# a <- clustering_hierarchical(distances, n_clust = 7)
# a <- clustering_hierarchical(distances, cut_height = c(.3, .4))
#
# plot(a$hclust$final.tree)
#
#
#
# a <- clustering_hierarchical(distances, n_clust = 7)
# a <- clustering_hierarchical(distances, cut_height = 0.3)
#
#
# #   ## 1. Controls ----
# #   if(!is.matrix(dat) & !(class(dat) == "dist")){
# #     stop("dat should be either a matrix with sites as rows and species as
# #          columns or a 'dist' object.")
# #   }
# #
# #   if(!(all(method %in% c("ward.D", "ward.D2", "single", "complete",
# #                          "average", "mcquitty", "median", "centroid", "dbscan",
# #                          "gmm", "diana", "pam")))){
# #     stop("Hierarchical clustering method chosen is not available.
# #      Please chose among the followings:
# #          'kmeans', 'meanshift',
# #          'ward.D', 'ward.D2', 'single', 'complete', 'average',
# #          'mcquitty', 'median', 'centroid', 'dbscan', 'gmm', 'diana' or
# #          'pam'.")
# #   }
# #
# #   if(!(optim_method %in% c("globalmax", "firstmax", "Tibs2001SEmax",
# #                            "firstSEmax", "globalSEmax."))){
# #     stop("Chosen gap statistic to determine the optimal number of cluster is
# #     not available.
# #      Please chose among the followings:
# #          globalmax, firstmax, Tibs2001SEmax, firstSEmax or globalSEmax.")
# #   }
# #
# #   if(!is.null(n_clust)){
# #     if(!(abs(n_clust - round(n_clust)) < .Machine$double.eps^0.5)){
# #       stop("n_clust must be an integer determining the number of clusters.")
# #     }
# #   }
# #
# #   if(is.null(n_clust) & any(
# #     method %in%
# #     c("kmeans", "ward.D", "ward.D2", "single", "complete",
# #       "average", "mcquitty", "median", "centroid", "diana", "pam"))){
# #     warning("One of the chosen method is a supervised algorithm that needs a
# #     number of clusters. Since 'n_clust = NULL', an optimization algorithm is
# #     first executed to determine the optimal numbers of clusters.
# #     This step may take a while.")
# #   }
# #
# #   if(!(abs(nstart - round(nstart)) < .Machine$double.eps^0.5)){
# #     stop("nstart must be an integer determining the number of random centroids
# #          to start k-means analysis.")
# #   }
# #
# #   if(!(abs(B - round(B)) < .Machine$double.eps^0.5)){
# #     stop("B must be an integer determining the number of Monte Carlo bootstrap
# #          samples.")
# #   }
# #
# #   if(!is.numeric(K.max)){
# #     stop("K.max must be a numeric determining the maximum number of clusters
# #          to consider.")
# #   }
# #
# #   if(is.matrix(dat)){
# #     if(K.max > nrow(dat)){
# #       stop("K.max should not be superior to the number of sites.")
# #     }
# #   }else if(class(dat) == "dist"){
# #     if(K.max > length(labels(dat))){
# #       warning("K.max should not be superior to the number of sites, reducing to number of sites.")
# #       K.max <- length(labels(dat))
# #     }
# #   }
# #
# #   if(length(method) > 1 &
# #      any(c("meanshift", "dbscan", "diana") %in% method) &
# #      !(all(method %in% c("meanshift", "dbscan", "diana")))){
# #     stop("'meanshift', 'dbscan' and 'diana' need a different input than the
# #     other clustering techniques. You should run cluster() separately with
# #          these methods.")
# #   }
# #
# #   ## 2. Input conversion ----
# #   if(!(all(method %in% c("meanshift", "dbscan", "diana")))){
# #     if(is.null(n_clust)){
# #       # Storing the matrix input, necessary to find the optimal nb of clusters
# #       dat1 <- dat
# #     }
# #     if(class(dat) != "dist"){
# #       # Project dat using simil() function with simpson metric
# #       dat <- simil(dat, metric = "simpson", input = "matrix", output = "dist")
# #     }
# #   }
# #
# #   ## 3. Clustering ----
# #   # Initiate empty list
# #   cluster_list <- list()
# #   # for-loop over the selected algorithms
# #   for(i in 1:length(method)){
# #
# #      if(method[i] == "dbscan"){
# #       require(dbscan)
# #       # Size of the epsilon neighborhood: normally determined by observing the
# #       # knee-plot
# #       eps <- quantile(kNNdist(dat, k = 5)[
# #         kNNdist(dat, k = 5) > 0], 0.25)
# #       # eps <- mean(kNNdist(dat, k = 5))
# #       db_clust <- dbscan(x = dat, eps = eps, minPts = 5)
# #
# #       # Data.frame of results
# #       res <- data.frame(site = rownames(dat),
# #                         cluster = as.character(db_clust$cluster))
# #
# #     } else if(method[i] == "gmm"){
# #       # Conversion to dist object
# #       dat <- as.dist(dat)
# #       require(mclust)
# #       gmm_mclust <- Mclust(dat)
# #
# #       # Data.frame of results
# #       res <- data.frame(site = names(gmm_mclust$classification),
# #                         cluster = as.character(gmm_mclust$classification))
# #
# #     } else{
# #       # Conversion to dist object
# #       # dat <- as.dist(dat)
# #
# #       # Number of clusters for supervised algorithms
# #       if(!is.null(n_clust)){
# #         optim_k <- n_clust
# #       } else{
# #         # Determine optimal numbers of clusters
# #         if(method[i] == "diana"){
# #           gap_stat <- cluster::clusGap(dat, FUN = kmeans, nstart = nstart, K.max = K.max,
# #                               B = B)
# #         } else{
# #           gap_stat <- cluster::clusGap(dat1, FUN = stats::kmeans, nstart = nstart, K.max = K.max,
# #                               B = B)
# #         }
# #
# #         optim_k <- maxSE(f = gap_stat$Tab[, "gap"],
# #                          SE.f = gap_stat$Tab[, "SE.sim"],
# #                          method = optim_method)
# #       }
# #
# #       if(!(method[i] %in% c("kmeans", "pam"))){
# #         if(method[i] == "diana"){
# #           h <- cluster::diana(dat)
# #         } else{
# #           # h <- hclust(dat, method = method[i])
# #           require(fastcluster)
# #           h <- fastcluster::hclust(dat, method = method[i])
# #           # plot(h)
# #         }
# #
# #         # Cut the tree with optim_k numbers
# #         dend <- stats::as.dendrogram(h)
# #         # Data.frame of results
# #         res <- data.frame(site = names(dendextend::cutree(dend,
# #                                                           k = optim_k)),
# #                           cluster = as.character(dendextend::cutree(dend,
# #                                                                     k = optim_k)))
# #       } else if(method[i] == "kmeans"){
# #         h <- kmeans(dat, centers = optim_k, iter.max = 10, nstart = nstart)
# #         res <- data.frame(site = names(h$cluster),
# #                           cluster = as.numeric(h$cluster))
# #       } else if(method[i] == "pam"){
# #         h <- cluster::pam(dat, k = optim_k, metric = "euclidean")
# #         res <- data.frame(site = names(h$clustering),
# #                           cluster = as.numeric(h$clustering))
# #       }
# #     }
# #
# #     # Convert cluster column into character
# #     res$cluster <- as.character(res$cluster)
# #
# #     # Changing column names of res: paste(cluster, method[i])
# #     colnames(res)[colnames(res) == "cluster"] <- paste0("cluster_", method[i])
# #
# #     # return(res)
# #     cluster_list[[i]] <- res
# #   }
# #   # Combine all list elements into one data.frame
# #   cluster_list <- Reduce(merge, cluster_list)
# #   return(cluster_list)
# #
# # }
#
