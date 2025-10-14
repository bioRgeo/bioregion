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
  expect_equal(detect_data_type_from_metric("Euclidean"), "unknown")
  
  # NA and NULL return unknown
  expect_equal(detect_data_type_from_metric(NA), "unknown")
  expect_equal(detect_data_type_from_metric(NULL), "unknown")
  
  # Custom or unknown metrics
  expect_equal(detect_data_type_from_metric("custom_metric"), "unknown")
  expect_equal(detect_data_type_from_metric("a/(a+b+c)"), "unknown")
  expect_equal(detect_data_type_from_metric("unknown"), "unknown")
  
})

# Test determine_weight_usage --------------------------------------------------
test_that("determine_weight_usage respects user_weight priority (Priority 1)", {
  
  # Create a mock bioregionalization object
  bio <- list(
    inputs = list(
      data_type = "occurrence",
      bipartite = FALSE
    ),
    args = list(index = "Jaccard")
  )
  
  # User explicitly sets weight = TRUE (abundance)
  result <- determine_weight_usage(
    bioregionalization = bio,
    user_weight = TRUE,
    user_index = NULL,
    net = NULL,
    verbose = FALSE
  )
  expect_equal(result$use_weight, TRUE)
  expect_equal(result$source, "user_argument")
  
  # User explicitly sets weight = FALSE (occurrence)
  result <- determine_weight_usage(
    bioregionalization = bio,
    user_weight = FALSE,
    user_index = NULL,
    net = NULL,
    verbose = FALSE
  )
  expect_equal(result$use_weight, FALSE)
  expect_equal(result$source, "user_argument")
  
})

test_that("determine_weight_usage uses user_weight with user_index", {
  
  bio <- list(
    inputs = list(data_type = "occurrence", bipartite = FALSE),
    args = list(index = "Jaccard")
  )
  
  # User provides both weight and index
  result <- determine_weight_usage(
    bioregionalization = bio,
    user_weight = TRUE,
    user_index = "Weight",
    net = NULL,
    verbose = FALSE
  )
  expect_equal(result$use_weight, TRUE)
  expect_equal(result$weight_col, "Weight")
  expect_equal(result$source, "user_argument")
  
})

test_that("determine_weight_usage uses user_weight with net", {
  
  bio <- list(
    inputs = list(data_type = "occurrence", bipartite = FALSE),
    args = list(index = "Jaccard")
  )
  
  net <- data.frame(
    Site = c("A", "A", "B"),
    Species = c("sp1", "sp2", "sp1"),
    Weight = c(10, 20, 15)
  )
  
  # User provides weight with net
  result <- determine_weight_usage(
    bioregionalization = bio,
    user_weight = TRUE,
    user_index = NULL,
    net = net,
    verbose = FALSE
  )
  expect_equal(result$use_weight, TRUE)
  expect_equal(result$weight_col, 3)
  expect_equal(result$source, "user_argument")
  
})

test_that("determine_weight_usage respects user_index priority (Priority 2)", {
  
  bio <- list(
    inputs = list(data_type = "occurrence", bipartite = FALSE),
    args = list(index = "Jaccard")
  )
  
  # User provides index (implies weight = TRUE)
  result <- determine_weight_usage(
    bioregionalization = bio,
    user_weight = NULL,
    user_index = "Abundance",
    net = NULL,
    verbose = FALSE
  )
  expect_equal(result$use_weight, TRUE)
  expect_equal(result$weight_col, "Abundance")
  expect_equal(result$source, "user_index")
  
  # User provides numeric index
  result <- determine_weight_usage(
    bioregionalization = bio,
    user_weight = NULL,
    user_index = 3,
    net = NULL,
    verbose = FALSE
  )
  expect_equal(result$use_weight, TRUE)
  expect_equal(result$weight_col, 3)
  expect_equal(result$source, "user_index")
  
})

test_that("determine_weight_usage uses data_type occurrence (Priority 3)", {
  
  bio <- list(
    inputs = list(data_type = "occurrence", bipartite = FALSE),
    args = list(index = "Jaccard")
  )
  
  # No user arguments, should use data_type
  result <- determine_weight_usage(
    bioregionalization = bio,
    user_weight = NULL,
    user_index = NULL,
    net = NULL,
    verbose = FALSE
  )
  expect_equal(result$use_weight, FALSE)
  expect_equal(result$weight_col, NULL)
  expect_equal(result$source, "bioregionalization_data_type")
  
})

test_that("determine_weight_usage uses data_type abundance (Priority 3)", {
  
  bio <- list(
    inputs = list(data_type = "abundance", bipartite = FALSE),
    args = list(index = "Bray")
  )
  
  # No user arguments, should use data_type
  result <- determine_weight_usage(
    bioregionalization = bio,
    user_weight = NULL,
    user_index = NULL,
    net = NULL,
    verbose = FALSE
  )
  expect_equal(result$use_weight, TRUE)
  expect_equal(result$source, "bioregionalization_data_type")
  
})

test_that("determine_weight_usage uses data_type abundance with bipartite", {
  
  bio <- list(
    inputs = list(data_type = "abundance", bipartite = TRUE),
    args = list(index = "Weight")
  )
  
  # Bipartite with abundance should extract weight column from args$index
  result <- determine_weight_usage(
    bioregionalization = bio,
    user_weight = NULL,
    user_index = NULL,
    net = NULL,
    verbose = FALSE
  )
  expect_equal(result$use_weight, TRUE)
  expect_equal(result$weight_col, "Weight")
  expect_equal(result$source, "bioregionalization_data_type")
  
})

test_that("determine_weight_usage uses data_type abundance with net", {
  
  bio <- list(
    inputs = list(data_type = "abundance", bipartite = FALSE),
    args = list(index = "Bray")
  )
  
  net <- data.frame(
    Site = c("A", "A", "B"),
    Species = c("sp1", "sp2", "sp1"),
    Abund = c(10, 20, 15)
  )
  
  # Should default to column 3 when net is provided
  result <- determine_weight_usage(
    bioregionalization = bio,
    user_weight = NULL,
    user_index = NULL,
    net = net,
    verbose = FALSE
  )
  expect_equal(result$use_weight, TRUE)
  expect_equal(result$weight_col, 3)
  expect_equal(result$source, "bioregionalization_data_type")
  
})

test_that("determine_weight_usage handles unknown data_type (Priority 4)", {
  
  bio <- list(
    inputs = list(data_type = "unknown", bipartite = FALSE),
    args = list(index = "Euclidean")
  )
  
  net <- data.frame(
    Site = c("A", "A", "B"),
    Species = c("sp1", "sp2", "sp1"),
    Weight = c(10, 20, 15)
  )
  
  # Unknown data_type should fall through to auto-detection
  result <- determine_weight_usage(
    bioregionalization = bio,
    user_weight = NULL,
    user_index = NULL,
    net = net,
    verbose = FALSE
  )
  expect_equal(result$use_weight, TRUE)
  expect_equal(result$weight_col, 3)
  expect_equal(result$source, "auto_detect_net")
  
})

test_that("determine_weight_usage auto-detects from net (Priority 4)", {
  
  bio <- list(
    inputs = list(bipartite = FALSE),
    args = list(index = "custom")
  )
  
  net_numeric <- data.frame(
    Site = c("A", "A", "B"),
    Species = c("sp1", "sp2", "sp1"),
    Count = c(10, 20, 15)
  )
  
  # Should detect numeric 3rd column
  result <- determine_weight_usage(
    bioregionalization = bio,
    user_weight = NULL,
    user_index = NULL,
    net = net_numeric,
    verbose = FALSE
  )
  expect_equal(result$use_weight, TRUE)
  expect_equal(result$weight_col, 3)
  expect_equal(result$source, "auto_detect_net")
  
  net_non_numeric <- data.frame(
    Site = c("A", "A", "B"),
    Species = c("sp1", "sp2", "sp1"),
    Notes = c("x", "y", "z")
  )
  
  # Should NOT detect non-numeric 3rd column
  expect_message(
    result <- determine_weight_usage(
      bioregionalization = bio,
      user_weight = NULL,
      user_index = NULL,
      net = net_non_numeric,
      verbose = TRUE
    ),
    "Could not determine if original data was occurrence or abundance-based"
  )
  expect_equal(result$use_weight, FALSE)
  expect_equal(result$source, "default")
  
})

test_that("determine_weight_usage defaults appropriately (Priority 5)", {
  
  bio <- list(
    inputs = list(bipartite = FALSE),
    args = list(index = "custom")
  )
  
  # No information available, should default to FALSE
  expect_message(
    result <- determine_weight_usage(
      bioregionalization = bio,
      user_weight = NULL,
      user_index = NULL,
      net = NULL,
      verbose = TRUE
    ),
    "Could not determine if original data was occurrence or abundance-based"
  )
  expect_equal(result$use_weight, FALSE)
  expect_equal(result$weight_col, NULL)
  expect_equal(result$source, "default")
  
})

test_that("determine_weight_usage handles NULL bioregionalization inputs", {
  
  bio <- list(
    inputs = list(),
    args = list()
  )
  
  net <- data.frame(
    Site = c("A", "A", "B"),
    Species = c("sp1", "sp2", "sp1"),
    Weight = c(10, 20, 15)
  )
  
  # Should fall back to auto-detection
  result <- determine_weight_usage(
    bioregionalization = bio,
    user_weight = NULL,
    user_index = NULL,
    net = net,
    verbose = FALSE
  )
  expect_equal(result$use_weight, TRUE)
  expect_equal(result$weight_col, 3)
  expect_equal(result$source, "auto_detect_net")
  
})

test_that("determine_weight_usage verbose messages work correctly", {
  
  bio <- list(
    inputs = list(data_type = "occurrence", bipartite = FALSE),
    args = list(index = "Jaccard")
  )
  
  # Test verbose = TRUE
  expect_message(
    determine_weight_usage(
      bioregionalization = bio,
      user_weight = TRUE,
      user_index = NULL,
      net = NULL,
      verbose = TRUE
    ),
    "Using weight specification from user argument"
  )
  
  expect_message(
    determine_weight_usage(
      bioregionalization = bio,
      user_weight = NULL,
      user_index = "Weight",
      net = NULL,
      verbose = TRUE
    ),
    "Weight column specified via 'index' argument"
  )
  
  expect_message(
    determine_weight_usage(
      bioregionalization = bio,
      user_weight = NULL,
      user_index = NULL,
      net = NULL,
      verbose = TRUE
    ),
    "Original data was occurrence-based"
  )
  
})

test_that("determine_weight_usage detects conflict between user and data_type", {
  
  bio <- list(
    inputs = list(data_type = "occurrence", bipartite = FALSE),
    args = list(index = "Jaccard")
  )
  
  # User sets weight = TRUE but data_type is occurrence
  expect_message(
    determine_weight_usage(
      bioregionalization = bio,
      user_weight = TRUE,
      user_index = NULL,
      net = NULL,
      verbose = TRUE
    ),
    "User-specified weight.*is different from.*original data type"
  )
  
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














