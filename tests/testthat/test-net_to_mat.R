# Inputs -----------------------------------------------------------------------
net <- data.frame(
  Site = c(rep("A", 2), rep("B", 3), rep("C", 2)),
  Species = c("a", "b", "a", "c", "d", "b", "d"),
  Weight = c(10, 100, 1, 20, 50, 10, 0)
)

net2 <- data.frame(
  Site = c(rep("A", 2), rep("B", 3), rep("C", 2)),
  Species = c("a", "b", "a", "c", "d", "b", "d"),
  Weight = c("a", "b", "a", "c", "d", "b", "d")
)

net3 <- rbind(net, net)

net4 <- net
net4[1,1] <- NA

net5 <- net
net5[1,3] <- NA

net6 <- data.frame(
  Obj1 = c(1,2,3,3),
  Obj2 = c(1,2,3,4),
  Weight = c(1,-1,0,2)
)

# Tests for valid outputs ------------------------------------------------------
test_that("valid outputs", {
  
  mat <- net_to_mat(net)
  expect_equal(class(mat), c("matrix",  "array"))
  expect_equal(dim(mat)[1], 3)
  expect_equal(dim(mat)[2], 4)
  
  mat <- net_to_mat(net,                        
                    weight = FALSE, 
                    squared = FALSE,
                    symmetrical = FALSE, 
                    missing_value = 0)
  expect_equal(dim(mat)[1], 3)
  expect_equal(dim(mat)[2], 4)
  expect_equal(min(mat), 0)
  expect_equal(max(mat), 1)
  
  mat <- net_to_mat(net,                        
                    weight = TRUE, 
                    squared = FALSE,
                    symmetrical = FALSE, 
                    missing_value = -1)
  expect_equal(dim(mat)[1], 3)
  expect_equal(dim(mat)[2], 4)
  expect_equal(min(mat), -1)
  expect_equal(max(mat), max(net[,3]))
  
  mat <- net_to_mat(net,                        
                    weight = TRUE, 
                    squared = TRUE,
                    symmetrical = FALSE, 
                    missing_value = -1)
  expect_equal(dim(mat)[1], 7)
  expect_equal(dim(mat)[2], 7)
  expect_equal(min(mat), -1)
  expect_equal(max(mat), max(net[,3]))
  
  mat <- net_to_mat(net,                        
                    weight = TRUE, 
                    squared = TRUE,
                    symmetrical = TRUE, 
                    missing_value = -1)
  expect_equal(dim(mat)[1], 7)
  expect_equal(dim(mat)[2], 7)
  expect_equal(min(mat), -1)
  expect_equal(max(mat), max(net[,3]))
  expect_equal(sum(sort(mat[upper.tri(mat)]) != sort(mat[upper.tri(mat)])), 0)
  
  mat <- net_to_mat(net6,                        
                    weight = TRUE, 
                    squared = TRUE,
                    symmetrical = FALSE, 
                    missing_value = -1)
  expect_equal(dim(mat)[1], 4)
  expect_equal(dim(mat)[2], 4)
  expect_equal(min(mat), -1)
  expect_equal(max(mat), max(net6[,3]))
  
  mat <- net_to_mat(net6,                        
                    weight = TRUE, 
                    squared = TRUE,
                    symmetrical = TRUE, 
                    missing_value = -1)
  expect_equal(dim(mat)[1], 4)
  expect_equal(dim(mat)[2], 4)
  expect_equal(min(mat), -1)
  expect_equal(max(mat), max(net6[,3]))
  expect_equal(sum(sort(mat[upper.tri(mat)]) != sort(mat[upper.tri(mat)])), 0)
  
})

# Tests for invalid inputs -----------------------------------------------------
test_that("indalid inputs", {
  
  expect_error(
    net_to_mat(net, weight = "zz"),
    "weight must be a boolean.", 
    fixed = TRUE)
  
  expect_error(
    net_to_mat(net, weight = c("zz",1)),
    "weight must be of length 1.", 
    fixed = TRUE)
  
  expect_error(
    net_to_mat(net, squared = "zz"),
    "squared must be a boolean.", 
    fixed = TRUE)
  
  expect_error(
    net_to_mat(net, squared = c("zz",1)),
    "squared must be of length 1.", 
    fixed = TRUE)
  
  expect_error(
    net_to_mat(net, symmetrical = "zz"),
    "symmetrical must be a boolean.", 
    fixed = TRUE)
  
  expect_error(
    net_to_mat(net, symmetrical = c("zz",1)),
    "symmetrical must be of length 1.", 
    fixed = TRUE)
  
  expect_error(
    net_to_mat(net, missing_value = "zz"),
    "missing_value must be numeric.", 
    fixed = TRUE)
  
  expect_error(
    net_to_mat(net, missing_value = c("zz",1)),
    "missing_value must be of length 1.", 
    fixed = TRUE)
  
  expect_error(
    net_to_mat(net, squared = FALSE, symmetrical = TRUE),
    "symmetrical only for squared matrix!", 
    fixed = TRUE)
  
  expect_error(
    net_to_mat("1"),
    "net must be a data.frame.", 
    fixed = TRUE)
  
  expect_error(
    net_to_mat(data.frame(net$Site)),
    "net must be a data.frame with at least two columns.", 
    fixed = TRUE)
  
  expect_error(
    net_to_mat(net3),
    "The first two columns of net contain duplicated pairs of nodes!", 
    fixed = TRUE)
  
  expect_error(
    net_to_mat(net4),
    "NA(s) detected in net.", 
    fixed = TRUE)
  
  expect_error(
    net_to_mat(net[,-3], weight = TRUE),
    "net must be a data.frame with at least three columns if weight equal 
        TRUE.", 
    fixed = TRUE)
  
  expect_error(
    net_to_mat(net5, weight = TRUE),
    "NA(s) detected in the weight column.", 
    fixed = TRUE)
  
  expect_error(
    net_to_mat(net2, weight = TRUE),
    "The weight column must be numeric.", 
    fixed = TRUE)

})
