# Inputs -----------------------------------------------------------------------
comat <- matrix(sample(0:1000, size = 50, replace = TRUE,
                       prob = 1 / 1:1001), 5, 10)
rownames(comat) <- paste0("Site", 1:5)
colnames(comat) <- paste0("Species", 1:10)

comat0 <- matrix(sample(0:1000, size = 50, replace = TRUE,
                       prob = 1 / 1:1001), 5, 10)

comat1 <- comat
rownames(comat1)[2] <- rownames(comat1)[3]

comat2 <- comat
colnames(comat2)[2] <- colnames(comat2)[3]

comat3 <- comat
comat3[1,1] <- NA

comat4 <- comat
comat4[1,1] <- -1


# Tests for valid outputs ------------------------------------------------------
test_that("valid output", {
  simil <- similarity(comat0, metric = c("abc", "ABC", "Euclidean"))
  expect_equal(dim(simil)[1], 10)
  expect_equal(dim(simil)[2], 9)
  expect_equal(simil$Site1, c(rep("1",4), rep("2",3), rep("3",2),"4"))

  simil <- similarity(comat, metric = c("abc", "ABC", "Euclidean"))
  expect_equal(dim(simil)[1], 10)
  expect_equal(dim(simil)[2], 9)
  expect_equal(simil$Site1, paste0("Site", 
                                   c(rep("1",4), rep("2",3), rep("3",2),"4")))
  
  
  
})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  expect_error(
    similarity(comat, metric = 1),
    "metric must be a character.",
    fixed = TRUE)
  
  expect_error(
    similarity(comat, method = 1),
    "method must be a character.",
    fixed = TRUE)
  
  expect_error(
    similarity(comat, formula = 1),
    "formula must be a character.",
    fixed = TRUE)
  
  expect_error(
    similarity(comat, metric = NULL, formula = NULL),
    "metric or formula should be used.",
    fixed = TRUE)
  
  expect_error(
    similarity(comat, method = "zzz"),
    "The method is not available.
     Please chose among the followings:
         prodmat or loops.",
    fixed = TRUE)
  
  expect_error(
    similarity(comat, metric = "zzz"),
    "One or several metric(s) chosen is not available.
     Please chose among the followings:
         abc, Jaccard, Jaccardturn, Sorensen, Simpson, ABC, Bray, Brayturn or
         Euclidean.",
    fixed = TRUE)
  
  expect_error(
    similarity(comat=1),
    "comat must be a matrix", 
    fixed = TRUE)
  
  expect_error(
    similarity(comat3),
    "NA(s) detected in the matrix!", 
    fixed = TRUE)
  
  expect_error(
    similarity(comat2),
    "Duplicated colnames detected!", 
    fixed = TRUE)
  
  expect_error(
    similarity(comat1),
    "Duplicated rownames detected!", 
    fixed = TRUE)
  
  expect_error(
    similarity(comat4),
    "comat should contains only positive real: negative
         value detected!", 
    fixed = TRUE)
  
  expect_message(
    similarity(comat4, metric = "Euclidean"),
    "Negative value(s) detected in comat!", 
    fixed = TRUE)
  
  expect_message(
    similarity(comat0),
    "No rownames detected, they have been assigned automatically.", 
    fixed = TRUE)
  
})