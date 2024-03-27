# Inputs -----------------------------------------------------------------------
infomap_versiondispo <- c("2.1.0", "2.6.0", "2.7.1")

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  expect_error(
    install_binaries(binpath = 1),
    "binpath must be a character.", 
    fixed = TRUE)
  
  expect_error(
    install_binaries(infomap_version = 1),
    "infomap_version must be a character.", 
    fixed = TRUE)
  
  expect_error(
    install_binaries(infomap_version = rep("1",100)),
    paste0(
      "Please choose versions of Infomap in the list: ",
      paste(infomap_versiondispo, collapse = " ")
    ), 
    fixed = TRUE)
  
  expect_error(
    install_binaries(infomap_version = "1"),
    paste0(
      "Please choose versions of Infomap in the list: ",
      paste(infomap_versiondispo, collapse = " ")
    )
    , 
    fixed = TRUE)
})  
