# Check controls ---------------------------------------------------------------
test_that("type", {
  
  expect_error(
    controls(args=NULL, 
             data=NULL, 
             type = "zz"),
    "Control type not defined!",
    fixed = TRUE)
 
})

test_that("input_nhandhclu", {
  
  data <- 1
  expect_error(
    controls(args=NULL, 
             data=data, 
             type = "input_nhandhclu"),
    "^data is not a bioregion.pairwise object")
  
})

test_that("input_similarity", {
  
  comat <- matrix(sample(0:1000, 
                         size = 50, 
                         replace = TRUE, 
                         prob = 1 / 1:1001),
                  5, 10)
  rownames(comat) <- paste0("Site", 1:5)
  colnames(comat) <- paste0("Species", 1:10)
  simil <- similarity(comat, metric = "all", formula = "a + b")
  attr(simil, "type") <- NULL
  expect_message(
    controls(args=NULL, 
             data=simil, 
             type = "input_similarity"),
    "^simil is a bioregion.pairwise object but")
  attr(simil, "type") <- "dissimilarity"
  expect_error(
    controls(args=NULL, 
             data=simil, 
             type = "input_similarity"),
    "^simil seems to be a dissimilarity object")
  
})

test_that("input_dissimilarity", {
  
  comat <- matrix(sample(0:1000, 
                         size = 50, 
                         replace = TRUE, 
                         prob = 1 / 1:1001),
                  5, 10)
  rownames(comat) <- paste0("Site", 1:5)
  colnames(comat) <- paste0("Species", 1:10)
  dissimil <- dissimilarity(comat, metric = "all", formula = "a + b")
  attr(dissimil, "type") <- NULL
  expect_message(
    controls(args=NULL, 
             data=dissimil, 
             type = "input_dissimilarity"),
    "^dissimil is a bioregion.pairwise object but")
  attr(dissimil, "type") <- "similarity"
  expect_error(
    controls(args=NULL, 
             data=dissimil, 
             type = "input_dissimilarity"),
    "^dissimil seems to be a similarity object")
  
})

test_that("input_conversion_similarity", {
  
  comat <- matrix(sample(0:1000, 
                         size = 50, 
                         replace = TRUE, 
                         prob = 1 / 1:1001),
                  5, 10)
  rownames(comat) <- paste0("Site", 1:5)
  colnames(comat) <- paste0("Species", 1:10)
  simil <- similarity(comat, metric = "all", formula = "a + b")
  attr(simil, "type") <- NULL
  test <- 1
  expect_error(
    controls(args=NULL, 
             data=test, 
             type = "input_conversion_similarity"),
    "^test should be a bioregion.pairwise object created by")
  expect_error(
    controls(args=NULL, 
             data=simil, 
             type = "input_conversion_similarity"),
    "^simil is a bioregion.pairwise object but")
  attr(simil, "type") <- "dissimilarity"
  expect_error(
    controls(args=NULL, 
             data=simil, 
             type = "input_conversion_similarity"),
    "^simil is already composed of dissimilarity metrics.")
  
})

test_that("input_conversion_dissimilarity", {
  
  comat <- matrix(sample(0:1000, 
                         size = 50, 
                         replace = TRUE, 
                         prob = 1 / 1:1001),
                  5, 10)
  rownames(comat) <- paste0("Site", 1:5)
  colnames(comat) <- paste0("Species", 1:10)
  dissimil <- dissimilarity(comat, metric = "all", formula = "a + b")
  attr(dissimil, "type") <- NULL
  test <- 1
  expect_error(
    controls(args=NULL, 
             data=test, 
             type = "input_conversion_dissimilarity"),
    "^test should be a bioregion.pairwise object created by")
  expect_error(
    controls(args=NULL, 
             data=dissimil, 
             type = "input_conversion_dissimilarity"),
    "^dissimil is a bioregion.pairwise object but")
  attr(dissimil, "type") <- "similarity"
  expect_error(
    controls(args=NULL, 
             data=dissimil, 
             type = "input_conversion_dissimilarity"),
    "^dissimil is already composed of similarity metrics.")
  
})

test_that("input_net", {
  
  net <- 1
  expect_error(
    controls(args=NULL, 
             data=net, 
             type = "input_net"),
    "net must be a data.frame.")

  net <- data.frame(ID=c(1,2))
  expect_error(
    controls(args=NULL, 
             data=net, 
             type = "input_net"),
    "net must be a data.frame with at least two columns.")
  
  net <- data.frame(ID=c(1,1),ID2=c(1,1))
  expect_error(
    controls(args=NULL, 
             data=net, 
             type = "input_net"),
    "^The first two columns of net")
  
  net <- data.frame(ID=c(1,2),ID2=c(NA,1))
  expect_error(
    controls(args=NULL, 
             data=net, 
             type = "input_net"),
    "^NA")
  
})

test_that("input_net_directed", {
  
  dir <- c(1,1)
  net <- 1
  expect_error(
    controls(args=dir, 
             data=net, 
             type = "input_net_directed"),
    "dir must be of length 1.")
  
  dir <- FALSE
  net <- data.frame(N1=c(1,2),N2=c(2,1))
  expect_error(
    controls(args=dir, 
             data=net, 
             type = "input_net_directed"),
    "^net should not be directed")
  
  net <- data.frame(N1=c(1,2),N2=c(2,1))
  expect_error(
    controls(args=NULL, 
             data=net, 
             type = "input_net_isdirected"),
    "^The network is directed")
  
  net <- data.frame(N1=c(1,2),N2=c(1,1))
  expect_error(
    controls(args=NULL, 
             data=net, 
             type = "input_net_isloop"),
    "^The network contains self")
  
})

test_that("input_net_weight", {
  
  w <- c(1,1)
  net <- 1
  expect_error(
    controls(args=w, 
             data=net, 
             type = "input_net_weight"),
    "w must be of length 1.")
  
  w <- 1
  net <- 1
  expect_error(
    controls(args=w, 
             data=net, 
             type = "input_net_weight"),
    "w must be a boolean.")
  
  w <- TRUE
  net <- data.frame(N1=c(1,2),N2=c(2,1))
  expect_error(
    controls(args=w, 
             data=net, 
             type = "input_net_weight"),
    "^net must be a data.frame with at least three columns")
  
})

test_that("input_net_index", {
  
  ind <- c(1,1)
  net <- 1
  expect_error(
    controls(args=ind, 
             data=net, 
             type = "input_net_index"),
    "ind must be of length 1.")
  
  ind <- "N1"
  net <- data.frame(N1=c(1,2),N2=c(3,4),W=c(5,5))
  expect_error(
    controls(args=ind, 
             data=net, 
             type = "input_net_index"),
    "^If ind is a character")
  
  ind <- 0.1
  net <- data.frame(N1=c(1,2),N2=c(3,4),W=c(5,5))
  expect_error(
    controls(args=ind, 
             data=net, 
             type = "input_net_index"),
    "If ind is numeric, it should be an integer.")
  
  ind <- 2
  net <- data.frame(N1=c(1,2),N2=c(3,4),W=c(5,5))
  expect_error(
    controls(args=ind, 
             data=net, 
             type = "input_net_index"),
    "ind should be strictly higher than 2.")
  
  ind <- 4
  net <- data.frame(N1=c(1,2),N2=c(3,4),W=c(5,5))
  expect_error(
    controls(args=ind, 
             data=net, 
             type = "input_net_index"),
    "ind should be lower or equal to 3.")
  
  ind <- factor("N2")
  net <- data.frame(N1=c(1,2),N2=c(3,4),W=c(5,5))
  expect_error(
    controls(args=ind, 
             data=net, 
             type = "input_net_index"),
    "^ind should be numeric or character.")
  
})

test_that("input_net_index_value", {
  
  ind <- 3
  net <- data.frame(N1=c(1,2),N2=c(3,4),W=c(NA,5))
  expect_error(
    controls(args=ind, 
             data=net, 
             type = "input_net_index_value"),
    "^NA")
  
  ind <- 3
  net <- data.frame(N1=c(1,2),N2=c(3,4),W=c("a","a"))
  expect_error(
    controls(args=ind, 
             data=net, 
             type = "input_net_index_value"),
    "The weight column must be numeric.")
  
})

test_that("input_net_index_positive_value", {
  
  ind <- 3
  net <- data.frame(N1=c(1,2),N2=c(3,4),W=c(NA,5))
  expect_error(
    controls(args=ind, 
             data=net, 
             type = "input_net_index_positive_value"),
    "^NA")
  
  ind <- 3
  net <- data.frame(N1=c(1,2),N2=c(3,4),W=c("a","a"))
  expect_error(
    controls(args=ind, 
             data=net, 
             type = "input_net_index_positive_value"),
    "The weight column must be numeric.")
  
  ind <- 3
  net <- data.frame(N1=c(1,2),N2=c(3,4),W=c(-1,1))
  expect_error(
    controls(args=ind, 
             data=net, 
             type = "input_net_index_positive_value"),
    "^The weight column should contain.")
  
})

test_that("input_net_bip", {
  
  net <- data.frame(N1=c(1,2),N2=c(1,2))
  expect_error(
    controls(args=NULL, 
             data=net, 
             type = "input_net_bip"),
    "The network is not bipartite!")
  
})

test_that("input_net_bip_col", {
  
  col <- "test"
  net <- data.frame(N1=c(1,2),N2=c(1,2))
  expect_error(
    controls(args=col, 
             data=net, 
             type = "input_net_bip_col"),
    "^If col is a character, it should be")
  
  col <- 3
  net <- data.frame(N1=c(1,2),N2=c(1,2))
  expect_error(
    controls(args=col, 
             data=net, 
             type = "input_net_bip_col"),
    "^If col is numeric, it should be")
  
  col <- TRUE
  net <- data.frame(N1=c(1,2),N2=c(1,2))
  expect_error(
    controls(args=col, 
             data=net, 
             type = "input_net_bip_col"),
    "col should be numeric or character.")
  
})

test_that("input_matrix", {
  
  mat <- 1
  expect_error(
    controls(args=NULL, 
             data=mat, 
             type = "input_matrix"),
    "mat must be a matrix.")
  
  mat <- matrix(1, 10, 10)
  rownames(mat) <- rep("1", 10)
  expect_error(
    controls(args=NULL, 
             data=mat, 
             type = "input_matrix"),
    "Duplicated rownames detected!")
  
  mat <- matrix(1, 10, 10)
  colnames(mat) <- rep("1", 10)
  expect_error(
    controls(args=NULL, 
             data=mat, 
             type = "input_matrix"),
    "Duplicated colnames detected!")
  
  mat <- matrix(NA, 10, 10)
  expect_error(
    controls(args=NULL, 
             data=mat, 
             type = "input_matrix"),
    "^NA")
  
})

test_that("input_dist", {
  
  mat <- 1
  expect_error(
    controls(args=NULL, 
             data=mat, 
             type = "input_dist"),
    "mat must be a dist object.")
  
  mat <- as.dist(matrix(1, 10, 10))
  mat[1] <- "1"
  expect_error(
    controls(args=NULL, 
             data=mat, 
             type = "input_dist"),
    "mat must be numeric.")
  
  mat <- as.dist(matrix(NA, 10, 10))
  expect_error(
    controls(args=NULL, 
             data=mat, 
             type = "input_dist"),
    "^NA")
  
})

test_that("input_data_frame_nhandhclu", {
  
  mat <- 1
  expect_error(
    controls(args=NULL, 
             data=mat, 
             type = "input_data_frame_nhandhclu"),
    "mat must be a data.frame.")
  
  mat <- data.frame(N1=c(0,1),N2=c(2,3))
  expect_error(
    controls(args=NULL, 
             data=mat, 
             type = "input_data_frame_nhandhclu"),
    "mat must be a data.frame with at least three columns.")
  
  mat <- data.frame(N1=c(NA,1),N2=c(2,3),W=c(5,6))
  expect_error(
    controls(args=NULL, 
             data=mat, 
             type = "input_data_frame_nhandhclu"),
    "^NA")
  
  mat <- data.frame(N1=c(2,1),N2=c(2,3),W=c(5,6))
  expect_error(
    controls(args=NULL, 
             data=mat, 
             type = "input_data_frame_nhandhclu"),
    "mat contains rows with the same site on both columns!")
  
  mat <- data.frame(N1=c(2,2),N2=c(1,1),W=c(5,6))
  expect_error(
    controls(args=NULL, 
             data=mat, 
             type = "input_data_frame_nhandhclu"),
    "^The first two columns of mat contain")
  
  mat <- data.frame(N1=c(2,1),N2=c(1,2),W=c(5,6))
  expect_error(
    controls(args=NULL, 
             data=mat, 
             type = "input_data_frame_nhandhclu"),
    "^The first two columns of mat contain ")
  
})

test_that("input_data_frame", {
  
  mat <- 1
  expect_error(
    controls(args=NULL, 
             data=mat, 
             type = "input_data_frame"),
    "mat must be a data.frame.")
  
  mat <- data.frame(N1=c(NA,1),N2=c(2,3),W=c(5,6))
  expect_error(
    controls(args=NULL, 
             data=mat, 
             type = "input_data_frame"),
    "^NA")
  
})

test_that("character", {
  
  test <- c(1,1)
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "character"),
    "test must be of length 1."
  )
  
  test <- 1
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "character"),
    "test must be a character."
  )
  
  test <- c(1,1)
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "character_vector"),
    "test must be a character."
  )
  
})

test_that("boolean", {
  
  test <- c(1,1)
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "boolean"),
    "test must be of length 1."
  )
  
  test <- 1
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "boolean"),
    "test must be a boolean."
  )
  
  test <- c(1,1)
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "boolean_vector"),
    "test must be a boolean."
  )
  
})

test_that("numeric", {
  
  test <- c("1","1")
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "numeric"),
    "test must be of length 1."
  )
  
  test <- "1"
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "numeric"),
    "test must be numeric."
  )
  
  test <- c("1","1")
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "numeric_vector"),
    "test must be numeric."
  )
  
})

test_that("positive_numeric", {
  
  test <- c("1","1")
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "positive_numeric"),
    "test must be of length 1."
  )
  
  test <- "1"
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "positive_numeric"),
    "test must be numeric."
  )
  
  test <- -1
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "positive_numeric"),
    "test must be higher than 0."
  )
  
  test <- c("1","1")
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "positive_numeric_vector"),
    "test must be numeric."
  )
  
  test <- c(-1,1)
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "positive_numeric_vector"),
    "test must be composed of values higher than 0."
  )
  
})

test_that("strict_positive_numeric", {
  
  test <- c("1","1")
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "strict_positive_numeric"),
    "test must be of length 1."
  )
  
  test <- "1"
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "strict_positive_numeric"),
    "test must be numeric."
  )
  
  test <- -1
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "strict_positive_numeric"),
    "test must be strictly higher than 0."
  )
  
  test <- c("1","1")
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "strict_positive_numeric_vector"),
    "test must be numeric."
  )
  
  test <- c(-1,1)
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "strict_positive_numeric_vector"),
    "test must be composed of values strictly higher than 0."
  )
  
})

test_that("integer", {
  
  test <- c("1","1")
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "integer"),
    "test must be of length 1."
  )
  
  test <- "1"
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "integer"),
    "test must be numeric."
  )
  
  test <- 0.1
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "integer"),
    "test must be an integer."
  )
  
  test <- c(1,0.1)
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "integer_vector"),
    "test must be composed of integers."
  )
  
})

test_that("positive_integer", {
  
  test <- c("1","1")
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "positive_integer"),
    "test must be of length 1."
  )
  
  test <- "1"
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "positive_integer"),
    "test must be numeric."
  )
  
  test <- 0.1
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "positive_integer"),
    "test must be an integer."
  )
  
  test <- -1
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "positive_integer"),
    "test must be higher than 0."
  )
  
  test <- c(1,0.1)
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "positive_integer_vector"),
    "test must be composed of integers."
  )
  
  test <- c(-1,1)
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "positive_integer_vector"),
    "test must be composed of values higher than 0."
  )
  
})

test_that("strict_positive_integer", {
  
  test <- c("1","1")
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "strict_positive_integer"),
    "test must be of length 1."
  )
  
  test <- "1"
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "strict_positive_integer"),
    "test must be numeric."
  )
  
  test <- 0.1
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "strict_positive_integer"),
    "test must be an integer."
  )
  
  test <- -1
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "strict_positive_integer"),
    "test must be strictly higher than 0."
  )
  
  test <- 0
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "strict_positive_integer"),
    "test must be strictly higher than 0."
  )
  
  test <- c(1,0.1)
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "strict_positive_integer_vector"),
    "test must be composed of integers."
  )
  
  test <- c(-1,1)
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "strict_positive_integer_vector"),
    "test must be composed of values strictly higher than 0."
  )
  
  test <- c(0,1)
  expect_error(
    controls(args = test, 
             data = NULL, 
             type = "strict_positive_integer_vector"),
    "test must be composed of values strictly higher than 0."
  )
  
})

# Tests for species/sites metrics ----------------------------------------------
test_that("invalid inputs", {
  
  # Inputs
  n_sites  <- sample(100:1000, 1)     
  n_species <- sample(100:1000, 1) 
  n_clusters <- sample(10, 1)
  clusters_g <- sample(n_clusters, n_sites, replace = TRUE) 
  clusters_s <- sample(n_clusters, n_species, replace = TRUE) 
  comat <- matrix(runif(n_sites*n_species), n_sites, n_species)
  rownames(comat) <- 1:n_sites
  colnames(comat) <- 1:n_species
  comat[comat < runif(1)] <- 0
  sim <- matrix(runif(n_sites*n_sites), n_sites, n_sites)
  diag(sim) <- 1
  temp <- t(sim)
  sim[lower.tri(sim)] <- temp[lower.tri(sim)]
  rownames(sim) <- 1:n_sites
  colnames(sim) <- 1:n_sites

  # sb
  sb <- sbgc(clusters_g, 
             bioregion_metrics = c("Specificity", "NSpecificity", 
                                   "Fidelity", 
                                   "IndVal", "NIndVal", 
                                   "Rho", 
                                   "CoreTerms"),
             bioregionalization_metrics = c("P"),
             comat,
             type = "sb", 
             data = "both")
  s <- sample(n_species, 1)
  b <- sample(n_clusters, 1)

  comat_bin <- comat>0
  tab <- cbind(clusters_g, as.numeric(comat_bin[,s]))
  agtab11 <- aggregate(tab[,2], list(tab[,1]), sum)
  agtab12 <- aggregate(tab[,2], list(tab[,1]), length)
  
  n_sb <- agtab11[agtab11[,1]==b, 2]
  n_s <- sum(tab[,2])
  n_b <- sum(clusters_g==b)
  
  Specificity_occ <- n_sb / n_s
  NSpecificity_occ <- (n_sb / n_b) / sum(agtab11[,2]/agtab12[,2])
  Fidelity_occ <- n_sb / n_b
  IndVal_occ <- Specificity_occ*Fidelity_occ
  NIndVal_occ <- NSpecificity_occ*Fidelity_occ
  num <- n_sb - n_s*(n_b / n_sites)
  den <- sqrt((n_b*(n_sites - n_b)/(n_sites - 1))*
              (n_s / n_sites)*
              (1 - (n_s / n_sites)))
  Rho_occ <- num / den
  
  tab <- cbind(clusters_g, as.numeric(comat[,s]))
  agtab21 <- aggregate(tab[,2], list(tab[,1]), sum)
  agtab22 <- aggregate(tab[,2], list(tab[,1]), length)
  
  w_sb <- agtab21[agtab21[,1]==b, 2]
  w_s <- sum(tab[,2])
  w_b <- sum(comat[clusters_g==b,])
  
  Specificity_abund <- w_sb / w_s
  NSpecificity_abund <- (w_sb / n_b) / sum(agtab21[,2]/agtab12[,2])
  Fidelity_abund <- w_sb / w_b
  IndVal_abund <- Specificity_abund*Fidelity_occ
  NIndVal_abund <- NSpecificity_abund*Fidelity_occ

  num <- (w_sb / n_b) - mean(tab[,2])
  den <- sqrt(((n_sites - n_b)/(n_sites - 1))*
                (var(tab[,2]) / n_b))
  Rho_abund <- num / den
  
  test <- c(s, b,
            n_sb, n_s, n_b,
            Specificity_occ, 
            NSpecificity_occ,
            Fidelity_occ,
            IndVal_occ,
            NIndVal_occ,
            Rho_occ,
            w_sb, w_s, w_b,
            Specificity_abund, 
            NSpecificity_abund,
            Fidelity_abund,
            IndVal_abund,
            NIndVal_abund,
            Rho_abund)
  
  check <- sb$bioregion
  check <- check[check[,1]==s & check[,2]==b,]
  
  expect_equal(as.numeric(test), as.numeric(check))
  
  P_occ <- 1 - sum((agtab11[,2]/sum(agtab11[,2]))*(agtab11[,2]/sum(agtab11[,2])))
  P_abund <- 1 - sum((agtab21[,2]/sum(agtab21[,2]))*(agtab21[,2]/sum(agtab21[,2])))
  
  test <- c(s, P_occ, P_abund)
  
  check <- sb$bioregionalization
  check <- check[check[,1]==s,]
  
  expect_equal(as.numeric(test), as.numeric(check))
  
  # gc
  gc <- sbgc(clusters_s, 
             bioregion_metrics = c("Specificity", "NSpecificity", 
                                   "Fidelity", 
                                   "IndVal", "NIndVal", 
                                   "Rho", 
                                   "CoreTerms"),
             bioregionalization_metrics = c("P"),
             comat,
             type = "gc", 
             data = "both")
  g <- sample(n_sites, 1)
  c <- sample(n_clusters, 1)
  
  comat_bin <- t(comat>0)
  tab <- cbind(clusters_s, as.numeric(comat_bin[,g]))
  agtab11 <- aggregate(tab[,2], list(tab[,1]), sum)
  agtab12 <- aggregate(tab[,2], list(tab[,1]), length)
  
  n_gc <- agtab11[agtab11[,1]==c, 2]
  n_g <- sum(tab[,2])
  n_c <- sum(clusters_s==c)
  
  Specificity_occ <- n_gc / n_g
  NSpecificity_occ <- (n_gc / n_c) / sum(agtab11[,2]/agtab12[,2])
  Fidelity_occ <- n_gc / n_c
  IndVal_occ <- Specificity_occ*Fidelity_occ
  NIndVal_occ <- NSpecificity_occ*Fidelity_occ
  num <- n_gc - n_g*(n_c / n_species)
  den <- sqrt((n_c*(n_species - n_c)/(n_species - 1))*
                (n_g / n_species)*
                (1 - (n_g / n_species)))
  Rho_occ <- num / den
  
  comat <- t(comat)
  tab <- cbind(clusters_s, as.numeric(comat[,g]))
  agtab21 <- aggregate(tab[,2], list(tab[,1]), sum)
  agtab22 <- aggregate(tab[,2], list(tab[,1]), length)
  
  w_gc <- agtab21[agtab21[,1]==c, 2]
  w_g <- sum(tab[,2])
  w_c <- sum(comat[clusters_s==c,])
  
  Specificity_abund <- w_gc / w_g
  NSpecificity_abund <- (w_gc / n_c) / sum(agtab21[,2]/agtab12[,2])
  Fidelity_abund <- w_gc / w_c
  IndVal_abund <- Specificity_abund*Fidelity_occ
  NIndVal_abund <- NSpecificity_abund*Fidelity_occ
  
  num <- (w_gc / n_c) - mean(tab[,2])
  den <- sqrt(((n_species - n_c)/(n_species - 1))*
                (var(tab[,2]) / n_c))
  Rho_abund <- num / den
  
  test <- c(g, c,
            n_gc, n_g, n_c,
            Specificity_occ, 
            NSpecificity_occ,
            Fidelity_occ,
            IndVal_occ,
            NIndVal_occ,
            Rho_occ,
            w_gc, w_g, w_c,
            Specificity_abund, 
            NSpecificity_abund,
            Fidelity_abund,
            IndVal_abund,
            NIndVal_abund,
            Rho_abund)
  
  check <- gc$bioregion
  check <- check[check[,1]==g & check[,2]==c,]
  
  expect_equal(as.numeric(test), as.numeric(check))
  
  P_occ <- 1 - sum((agtab11[,2]/sum(agtab11[,2]))*(agtab11[,2]/sum(agtab11[,2])))
  P_abund <- 1 - sum((agtab21[,2]/sum(agtab21[,2]))*(agtab21[,2]/sum(agtab21[,2])))
  
  test <- c(g, P_occ, P_abund)
  
  check <- gc$bioregionalization
  check <- check[check[,1]==g,]
  
  expect_equal(as.numeric(test), as.numeric(check))
  
  # gb
  gb <- gb(as.character(clusters_g), 
           bioregion_metrics = c("MeanSim", "SdSim"),
           bioregionalization_metrics = c("Silhouette"),
           comat = NULL,
           sim,
           include_cluster = FALSE)
  
  tab <- cbind(clusters_g[-g], as.numeric(sim[-g,g]))
  agtab11 <- aggregate(tab[,2], list(tab[,1]), mean)
  agtab12 <- aggregate(tab[,2], list(tab[,1]), sd)
  
  test <- c(g, b,
            agtab11[agtab11[,1]==b,2],
            agtab12[agtab12[,1]==b,2])
  
  check <- gb$bioregion
  check <- check[check[,1]==g & check[,2]==b,]
  
  expect_equal(as.numeric(test), as.numeric(check))
  
  if(n_clusters == 1){
    sil <- NA
  }else{
    sil <- (agtab11[agtab11[,1]==clusters_g[g],2] - max(agtab11[agtab11[,1]!=clusters_g[g],2])) /
      max(agtab11[agtab11[,1]==clusters_g[g],2], max(agtab11[agtab11[,1]!=clusters_g[g],2]))
    sil[is.infinite(sil)] <- NA  
  }
  
  test <- c(g,sil)
  
  check <- gb$bioregionalization
  check <- check[check[,1]==g,]
  
  expect_equal(round(as.numeric(test),digits=2), 
               round(as.numeric(check),digits=2))
  

})

################################################################################

# Test detect_data_type_from_metric --------------------------------------------
test_that("detect_data_type_from_metric works with occurrence metrics", {
  
  # Standard occurrence metrics
  expect_equal(detect_data_type_from_metric("Jaccard"), "occurrence")
  expect_equal(detect_data_type_from_metric("Simpson"), "occurrence")
  expect_equal(detect_data_type_from_metric("Sorensen"), "occurrence")
  expect_equal(detect_data_type_from_metric("Jaccardturn"), "occurrence")
  expect_equal(detect_data_type_from_metric("abc"), "occurrence")
  
})

test_that("detect_data_type_from_metric works with abundance metrics", {
  
  # Standard abundance metrics
  expect_equal(detect_data_type_from_metric("Bray"), "abundance")
  expect_equal(detect_data_type_from_metric("ABC"), "abundance")
  expect_equal(detect_data_type_from_metric("Brayturn"), "abundance")
  
})

test_that("detect_data_type_from_metric works with betapart occurrence metrics", {
  
  # Betapart occurrence metrics (case-insensitive)
  expect_equal(detect_data_type_from_metric("beta.sim"), "occurrence")
  expect_equal(detect_data_type_from_metric("BETA.SIM"), "occurrence")
  expect_equal(detect_data_type_from_metric("Beta.Sim"), "occurrence")
  expect_equal(detect_data_type_from_metric("beta.sne"), "occurrence")
  expect_equal(detect_data_type_from_metric("beta.sor"), "occurrence")
  expect_equal(detect_data_type_from_metric("beta.jtu"), "occurrence")
  expect_equal(detect_data_type_from_metric("beta.jne"), "occurrence")
  expect_equal(detect_data_type_from_metric("beta.jac"), "occurrence")
  
})

test_that("detect_data_type_from_metric works with betapart abundance metrics", {
  
  # Betapart abundance metrics (case-insensitive)
  expect_equal(detect_data_type_from_metric("beta.bray.bal"), "abundance")
  expect_equal(detect_data_type_from_metric("BETA.BRAY.BAL"), "abundance")
  expect_equal(detect_data_type_from_metric("Beta.Bray.Bal"), "abundance")
  expect_equal(detect_data_type_from_metric("beta.bray.gra"), "abundance")
  expect_equal(detect_data_type_from_metric("beta.bray"), "abundance")
  expect_equal(detect_data_type_from_metric("beta.ruz.bal"), "abundance")
  expect_equal(detect_data_type_from_metric("beta.ruz.gra"), "abundance")
  expect_equal(detect_data_type_from_metric("beta.ruz"), "abundance")
  
})

test_that("detect_data_type_from_metric works with unknown metrics", {
  
  # Euclidean is explicitly unknown
  expect_equal(detect_data_type_from_metric("Euclidean"), NA)
  
  # NA and NULL return unknown
  expect_equal(detect_data_type_from_metric(NA), NA)
  expect_equal(detect_data_type_from_metric(NULL), NA)
  
  # Custom or unknown metrics
  expect_equal(detect_data_type_from_metric("custom_metric"), NA)
  expect_equal(detect_data_type_from_metric("a/(a+b+c)"), NA)
  expect_equal(detect_data_type_from_metric(NA), NA)
  
})


# Test betapart integration ----------------------------------------------------
test_that("detect_data_type_from_metric works with betapart occurrence metrics (beta.pair)", {
  
  skip_if_not_installed("betapart")
  
  # Create a small binary matrix
  comat_bin <- matrix(sample(0:1, 50, replace = TRUE), 5, 10)
  rownames(comat_bin) <- paste0("Site", 1:5)
  colnames(comat_bin) <- paste0("Species", 1:10)
  
  # Test with betapart::beta.pair (occurrence-based)
  beta_result <- betapart::beta.pair(comat_bin, index.family = "jaccard")
  
  # Convert to bioregion format
  dissim_beta <- as_bioregion_pairwise(beta_result, pkg = "betapart")
  
  # Check that betapart metric names are correctly detected as occurrence
  expect_equal(detect_data_type_from_metric("beta.jac"), "occurrence")
  expect_equal(detect_data_type_from_metric("beta.jtu"), "occurrence")
  expect_equal(detect_data_type_from_metric("beta.jne"), "occurrence")
  
  # Test with Sorensen family
  beta_result_sor <- betapart::beta.pair(comat_bin, index.family = "sorensen")
  dissim_beta_sor <- as_bioregion_pairwise(beta_result_sor, pkg = "betapart")
  
  expect_equal(detect_data_type_from_metric("beta.sor"), "occurrence")
  expect_equal(detect_data_type_from_metric("beta.sim"), "occurrence")
  expect_equal(detect_data_type_from_metric("beta.sne"), "occurrence")
  
})

test_that("detect_data_type_from_metric works with betapart abundance metrics (beta.pair.abund)", {
  
  skip_if_not_installed("betapart")
  
  # Create a small abundance matrix
  comat <- matrix(sample(0:100, 50, replace = TRUE), 5, 10)
  rownames(comat) <- paste0("Site", 1:5)
  colnames(comat) <- paste0("Species", 1:10)
  
  # Test with betapart::beta.pair.abund (abundance-based)
  beta_result <- betapart::beta.pair.abund(comat, index.family = "bray")
  
  # Convert to bioregion format
  dissim_beta <- as_bioregion_pairwise(beta_result, pkg = "betapart")
  
  # Check that betapart abundance metric names are correctly detected
  expect_equal(detect_data_type_from_metric("beta.bray"), "abundance")
  expect_equal(detect_data_type_from_metric("beta.bray.bal"), "abundance")
  expect_equal(detect_data_type_from_metric("beta.bray.gra"), "abundance")
  
  # Test with Ruzicka family
  beta_result_ruz <- betapart::beta.pair.abund(comat, index.family = "ruzicka")
  dissim_beta_ruz <- as_bioregion_pairwise(beta_result_ruz, pkg = "betapart")
  
  expect_equal(detect_data_type_from_metric("beta.ruz"), "abundance")
  expect_equal(detect_data_type_from_metric("beta.ruz.bal"), "abundance")
  expect_equal(detect_data_type_from_metric("beta.ruz.gra"), "abundance")
  
})

test_that("betapart metrics are case-insensitive", {
  
  # Test case insensitivity for occurrence metrics
  expect_equal(detect_data_type_from_metric("BETA.JAC"), "occurrence")
  expect_equal(detect_data_type_from_metric("Beta.Jac"), "occurrence")
  expect_equal(detect_data_type_from_metric("beta.JAC"), "occurrence")
  expect_equal(detect_data_type_from_metric("BETA.SOR"), "occurrence")
  expect_equal(detect_data_type_from_metric("Beta.Sim"), "occurrence")
  
  # Test case insensitivity for abundance metrics
  expect_equal(detect_data_type_from_metric("BETA.BRAY"), "abundance")
  expect_equal(detect_data_type_from_metric("Beta.Bray"), "abundance")
  expect_equal(detect_data_type_from_metric("beta.BRAY"), "abundance")
  expect_equal(detect_data_type_from_metric("BETA.RUZ.BAL"), "abundance")
  expect_equal(detect_data_type_from_metric("Beta.Ruz.Gra"), "abundance")
  
})

test_that("betapart integration with clustering functions preserves data_type", {
  
  skip_if_not_installed("betapart")
  
  # Create test matrices
  comat_bin <- matrix(sample(0:1, 100, replace = TRUE), 10, 10)
  rownames(comat_bin) <- paste0("Site", 1:10)
  colnames(comat_bin) <- paste0("Species", 1:10)
  
  comat_abund <- matrix(sample(0:50, 100, replace = TRUE), 10, 10)
  rownames(comat_abund) <- paste0("Site", 1:10)
  colnames(comat_abund) <- paste0("Species", 1:10)
  
  # Test occurrence-based betapart metrics
  beta_occ <- betapart::beta.pair(comat_bin, index.family = "jaccard")
  dissim_occ <- as_bioregion_pairwise(beta_occ, pkg = "betapart")
  
  # Run clustering with occurrence-based betapart metrics
  clust_occ <- nhclu_pam(dissim_occ, index = "beta.jac", n_clust = 3)
  expect_equal(clust_occ$inputs$data_type, "occurrence")
  expect_equal(clust_occ$inputs$pairwise_metric, "beta.jac")
  
  # Test abundance-based betapart metrics
  beta_abund <- betapart::beta.pair.abund(comat_abund, index.family = "bray")
  dissim_abund <- as_bioregion_pairwise(beta_abund, pkg = "betapart")
  
  # Run clustering with abundance-based betapart metrics
  clust_abund <- nhclu_pam(dissim_abund, index = "beta.bray", n_clust = 3)
  expect_equal(clust_abund$inputs$data_type, "abundance")
  expect_equal(clust_abund$inputs$pairwise_metric, "beta.bray")
  
})

test_that("betapart.core and betapart.core.abund work correctly", {
  
  skip_if_not_installed("betapart")
  
  # Create test matrices
  comat_bin <- matrix(sample(0:1, 100, replace = TRUE), 10, 10)
  rownames(comat_bin) <- paste0("Site", 1:10)
  colnames(comat_bin) <- paste0("Species", 1:10)
  
  comat_abund <- matrix(sample(0:50, 100, replace = TRUE), 10, 10)
  rownames(comat_abund) <- paste0("Site", 1:10)
  colnames(comat_abund) <- paste0("Species", 1:10)
  
  # Test betapart.core (occurrence) - converts to a, b, c format
  beta_core_occ <- betapart::betapart.core(comat_bin)
  dissim_core_occ <- as_bioregion_pairwise(beta_core_occ, pkg = "betapart")
  
  # Verify that a, b, c columns exist (converted from betapart occurrence format)
  expect_true("a" %in% colnames(dissim_core_occ))
  expect_true("b" %in% colnames(dissim_core_occ))
  expect_true("c" %in% colnames(dissim_core_occ))
  expect_true("min(b,c)" %in% colnames(dissim_core_occ))
  
  # Test betapart.core.abund (abundance) - converts to A and derived columns
  beta_core_abund <- betapart::betapart.core.abund(comat_abund)
  dissim_core_abund <- as_bioregion_pairwise(beta_core_abund, pkg = "betapart")
  
  # Verify that A and derived columns exist (converted from betapart abundance format)
  expect_true("A" %in% colnames(dissim_core_abund))
  expect_true("min(B,C)" %in% colnames(dissim_core_abund))
  expect_true("max(B,C)" %in% colnames(dissim_core_abund))
  expect_true("sum(B,C)" %in% colnames(dissim_core_abund))
  
})














