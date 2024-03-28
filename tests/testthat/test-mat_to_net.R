# Inputs -----------------------------------------------------------------------
mat <- matrix(sample(1000, 50), 5, 10)
rownames(mat) <- paste0("Site", 1:5)
colnames(mat) <- paste0("Species", 1:10)

mat2 <- matrix(sample(100), 10, 10)
mat2[1,1] <- 0
mat2[1,2] <- 0

mat3 <- matrix(sample(1000, 50), 5, 10)
rownames(mat3) <- paste0("Site", 1:5)
colnames(mat3) <- paste0("Species", 1:10)
mat3[1,1] <- NA 

mat4 <- mat3
colnames(mat4)[2] <- colnames(mat4)[1]

mat5 <- mat3
rownames(mat5)[2] <- rownames(mat5)[1]


# Tests for valid outputs ------------------------------------------------------
test_that("valid outputs", {
  
  net <- mat_to_net(mat, weight = TRUE)
  expect_equal(class(net), "data.frame")
  expect_equal(dim(net)[1], 50)
  expect_equal(dim(net)[2], 3)
  
  net <- mat_to_net(mat, weight = FALSE)
  expect_equal(class(net), "data.frame")
  expect_equal(dim(net)[1], 50)
  expect_equal(dim(net)[2], 2)
  
  net2 <- mat_to_net(mat2, 
                     weight = TRUE,
                     remove_zeroes = TRUE,
                     include_diag = TRUE,
                     include_lower = TRUE)
  expect_equal(dim(net2)[1], 98)
  expect_equal(dim(net2)[2], 3)
  
  net2 <- mat_to_net(mat2, 
                     weight = FALSE,
                     remove_zeroes = FALSE,
                     include_diag = FALSE,
                     include_lower = TRUE)
  expect_equal(dim(net2)[1], 90)
  expect_equal(dim(net2)[2], 2)
  
  net2 <- mat_to_net(mat2, 
                     weight = FALSE,
                     remove_zeroes = FALSE,
                     include_diag = FALSE,
                     include_lower = FALSE)
  expect_equal(dim(net2)[1], 45)
  expect_equal(dim(net2)[2], 2)
  
  net2 <- mat_to_net(mat2, 
                     weight = TRUE,
                     remove_zeroes = TRUE,
                     include_diag = FALSE,
                     include_lower = FALSE)
  expect_equal(dim(net2)[1], 44)
  expect_equal(dim(net2)[2], 3)
  
  net2 <- mat_to_net(mat2, 
                     weight = FALSE,
                     remove_zeroes = FALSE,
                     include_diag = TRUE,
                     include_lower = FALSE)
  expect_equal(dim(net2)[1], 55)
  expect_equal(dim(net2)[2], 2)
  

})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    mat_to_net(mat, weight = "zz"),
    "weight must be a boolean", 
    fixed = TRUE)
  
  expect_error(
    mat_to_net(mat, remove_zeroes = "zz"),
    "remove_zeroes must be a boolean", 
    fixed = TRUE)
  
  expect_error(
    mat_to_net(mat, include_diag = "zz"),
    "include_diag must be a boolean", 
    fixed = TRUE)
  
  expect_error(
    mat_to_net(mat, include_lower = "zz"),
    "include_lower must be a boolean", 
    fixed = TRUE)
  
  expect_message(
    mat_to_net(mat, include_diag = FALSE),
    "include_diag is only used with squared matrix.", 
    fixed = TRUE)
  
  expect_message(
    mat_to_net(mat, include_lower = FALSE),
    "include_lower is only used with squared matrix.", 
    fixed = TRUE)
  
  expect_error(
    mat_to_net("1"),
    "mat must be a matrix", 
    fixed = TRUE)
  
  expect_error(
    mat_to_net(mat3),
    "NA(s) detected in the matrix!", 
    fixed = TRUE)
  
  expect_error(
    mat_to_net(mat4),
    "Duplicated colnames detected!", 
    fixed = TRUE)
  
  expect_error(
    mat_to_net(mat5),
    "Duplicated rownames detected!", 
    fixed = TRUE)

})
