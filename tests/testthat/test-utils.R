# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    controls(args=NULL, data=NULL, type = "zz"),
    "Control type not defined!",
    fixed = TRUE)
 
})