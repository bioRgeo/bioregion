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
    "^data is not a bioregion.pairwise.metric object")
  
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
    "^simil is a bioregion.pairwise.metric object but")
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
    "^dissimil is a bioregion.pairwise.metric object but")
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
    "^test should be a bioregion.pairwise.metric object created by")
  expect_error(
    controls(args=NULL, 
             data=simil, 
             type = "input_conversion_similarity"),
    "^simil is a bioregion.pairwise.metric object but")
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
    "^test should be a bioregion.pairwise.metric object created by")
  expect_error(
    controls(args=NULL, 
             data=dissimil, 
             type = "input_conversion_dissimilarity"),
    "^dissimil is a bioregion.pairwise.metric object but")
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












