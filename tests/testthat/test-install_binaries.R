# Inputs -----------------------------------------------------------------------
infomap_versiondispo <- c("2.1.0", "2.6.0", "2.7.1", "2.8.0")

# Tests for valid outputs ------------------------------------------------------
test_that("valid output", {
  
  quietly(
    install_binaries(binpath = "tempdir",
                     download_only = FALSE,
                     infomap_version = c("2.1.0", 
                                         "2.6.0", 
                                         "2.7.1",
                                         "2.8.0"),
                     verbose = TRUE)
  )
  
  #quietly(
    install_binaries(binpath = "pkgfolder",
                     download_only = FALSE,
                     infomap_version = c("2.1.0", 
                                         "2.6.0", 
                                         "2.7.1",
                                         "2.8.0"),
                     verbose = FALSE)
  #)
  
})

# Tests for invalid inputs -----------------------------------------------------
test_that("invalid inputs", {
  
  expect_error(
    install_binaries(binpath = 1),
    "binpath must be a character.", 
    fixed = TRUE)
  
  expect_error(
    install_binaries(binpath = c(1,1)),
    "binpath must be of length 1.", 
    fixed = TRUE)
  
  expect_error(
    install_binaries(binpath = "test_271224"),
    "Impossible to access test_271224.", 
    fixed = TRUE)
  
  expect_error(
    install_binaries(binpath = "tempdir", 
                     download_only = 1),
    "download_only must be a boolean.", 
    fixed = TRUE)
  
  expect_error(
    install_binaries(binpath = "tempdir", 
                     download_only = c(TRUE,FALSE)),
    "download_only must be of length 1.", 
    fixed = TRUE)
  
  expect_error(
    install_binaries(binpath = "tempdir", 
                     download_only = TRUE),
    "download_only cannot be set to TRUE if binpath is tempdir or pkgfolder!", 
    fixed = TRUE)
  
  expect_error(
    install_binaries(binpath = "pkgfolder", 
                     download_only = TRUE),
    "download_only cannot be set to TRUE if binpath is tempdir or pkgfolder!", 
    fixed = TRUE)
  
  expect_error(
    install_binaries(infomap_version = 1),
    "infomap_version must be a character.", 
    fixed = TRUE)
  
  expect_error(
    install_binaries(infomap_version = paste0(1:100)),
    "^Please select a version of Infomap from")
  
  expect_error(
    install_binaries(infomap_version = "1"),
    "^Please select a version of Infomap from")
  
  expect_error(
    install_binaries(verbose = 1),
    "verbose must be a boolean.",
    fixed = TRUE)
  
  expect_error(
    install_binaries(verbose = c(TRUE, FALSE)),
    "verbose must be of length 1.",
    fixed = TRUE)
  
})  
