# Inputs -----------------------------------------------------------------------
comat <- matrix(sample(0:1000, size = 50, replace = TRUE,
                       prob = 1 / 1:1001), 5, 10)
rownames(comat) <- paste0("Site", 1:5)
colnames(comat) <- paste0("Species", 1:10)

vec1=as.numeric(comat[1,])
vec2=as.numeric(comat[2,])
a=sum((vec1*vec2)>0)
b=sum(vec1>0)-a
c=sum(vec2>0)-a
A=sum(pmin(vec1,vec2))
B=sum(vec1)-A
C=sum(vec2)-A

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
  expect_equal(inherits(simil, "bioregion.pairwise"), TRUE)
  expect_equal(dim(simil)[1], 10)
  expect_equal(dim(simil)[2], 9)
  expect_equal(simil$Site1, c(rep("1",4), rep("2",3), rep("3",2),"4"))

  simil <- similarity(comat, metric = c("abc", "ABC", "Euclidean"))
  expect_equal(dim(simil)[1], 10)
  expect_equal(dim(simil)[2], 9)
  expect_equal(simil$Site1, paste0("Site", 
                                   c(rep("1",4), rep("2",3), rep("3",2),"4")))
  
  # abc
  simil1 <- similarity(comat, metric = "abc", method =  "prodmat")
  simil2 <- similarity(comat, metric = "abc", method =  "loops")
  expect_equal(sum(simil1==simil2), 50)
  expect_equal(simil1$a[1], a)
  expect_equal(simil1$b[1], b)
  expect_equal(simil1$c[1], c)
  
  # ABC
  simil <- similarity(comat, metric = "ABC")
  expect_equal(simil$A[1], A)
  expect_equal(simil$B[1], B)
  expect_equal(simil$C[1], C)
  
  # Jaccard
  simil <- similarity(comat, metric = "Jaccard", formula = "1-(b+c)/(a+b+c)")
  expect_equal(simil[1,3], 1-(b+c)/(a+b+c))
  expect_equal(simil[1,4], 1-(b+c)/(a+b+c))
  expect_equal(sum(simil[,3]==simil[,4]), 10)
  
  # Jaccardturn
  simil <- similarity(comat, metric = "Jaccardturn", formula = "1-2*pmin(b,c)/(a+2*pmin(b,c))")
  expect_equal(simil[1,3], 1-2*min(b,c)/(a+2*min(b,c)))
  expect_equal(simil[1,4], 1-2*min(b,c)/(a+2*min(b,c)))
  expect_equal(sum(simil[,3]==simil[,4]), 10)
  
  # Sorensen
  simil <- similarity(comat, metric = "Sorensen", formula = "1-(b+c)/(2*a+b+c)")
  expect_equal(simil[1,3], 1-(b+c)/(2*a+b+c))
  expect_equal(simil[1,4], 1-(b+c)/(2*a+b+c))
  expect_equal(sum(simil[,3]==simil[,4]), 10)
  
  # Simpson
  simil <- similarity(comat, metric = "Simpson", formula = "1-pmin(b,c)/(a+pmin(b,c))")
  expect_equal(simil[1,3], 1-min(b,c)/(a+min(b,c)))
  expect_equal(simil[1,4], 1-min(b,c)/(a+min(b,c)))
  expect_equal(sum(simil[,3]==simil[,4]), 10)
  
  # Bray
  simil <- similarity(comat, metric = "Bray", formula = "1-(B+C)/(2*A+B+C)")
  expect_equal(simil[1,3], 1-(B+C)/(2*A+B+C))
  expect_equal(simil[1,4], 1-(B+C)/(2*A+B+C))
  expect_equal(sum(simil[,3]==simil[,4]), 10)
  
  # Brayturn
  simil <- similarity(comat, metric = "Brayturn", formula = "1-pmin(B,C)/(A+pmin(B,C))")
  expect_equal(simil[1,3], 1-min(B,C)/(A+min(B,C)))
  expect_equal(simil[1,4], 1-min(B,C)/(A+min(B,C)))
  expect_equal(sum(simil[,3]==simil[,4]), 10)
  
  # Euclidean
  simil <- similarity(comat, metric = "Euclidean")
  expect_equal(simil$Euclidean[1], 1/(1+sqrt(sum((vec1-vec2)^2))))
  
  # formula
  simil <- similarity(comat, metric = NULL, formula = c("a","b","c","A","B","C"))
  expect_equal(simil$a[1], a)
  expect_equal(simil$b[1], b)
  expect_equal(simil$c[1], c)
  expect_equal(simil$A[1], A)
  expect_equal(simil$B[1], B)
  expect_equal(simil$C[1], C)
  
})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    similarity(comat, 
               metric = 1),
    "metric must be a character.",
    fixed = TRUE)
  
  expect_error(
    similarity(comat, 
               method = 1),
    "method must be a character.",
    fixed = TRUE)
  
  expect_error(
    similarity(comat, 
               formula = 1),
    "formula must be a character.",
    fixed = TRUE)
  
  expect_error(
    similarity(comat, 
               metric = NULL, 
               formula = NULL),
    "metric or formula should be used.",
    fixed = TRUE)
  
  expect_error(
    similarity(comat, 
               method = "zzz"),
    "^Please choose method from the following:")
  
  expect_error(
    similarity(comat,
               metric = "zzz"),
    "^One or several metric")
  
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
    "Negative value detected in comat.")
  
  expect_message(
    similarity(comat4, metric = "Euclidean"),
    "Negative value(s) detected in comat!", 
    fixed = TRUE)
  
  expect_message(
    similarity(comat0),
    "No rownames detected, they have been assigned automatically.", 
    fixed = TRUE)
  
})
