#' Download, unzip, check permission and test the bioRgeo's executable binary
#' files
#'
#' This function downloads and unzips the 'bin' folder needed to run some
#' functions of bioRgeo. It also checks if the files have the permissions to
#' be executed as programs. It finally tests if the executable binary files
#' are running properly.
#'
#' @param binpath a `character` indicating the path to the folder that will
#' host the 'bin' folder containing the executable binary files (bioRgeo's
#' package by default, if you use a different folder please be sure to update
#' the `binpath` in [netclu_infomap], [netclu_louvain] and [netclu_oslom]).
#' 
#' @param infomap_version a `character` vector indicating the Infomap
#' version(s) to install.
#' 
#' @note
#' Only the Infomap version 2.1.0 and 2.6.0 are available for now.
#' 
#' @author
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}),
#' Boris Leroy (\email{leroy.boris@gmail.com}) and
#' Pierre Denelle (\email{pierre.denelle@gmail.com})
#' 
#' @export

install_binaries <- function(binpath = NULL,
                             infomap_version = c("2.1.0", "2.6.0")) {
  
  # Set binpath
  if (is.null(binpath)) {
    # Identify bioRgeo directory on your computer
    biodir <- list.dirs(.libPaths(), recursive = FALSE)
    binpath <- biodir[grep("bioRgeo", biodir)]
    
    if(length(binpath) > 1){
      message("Several bioRgeo directories have been detected in your default
      package/library folder(s). 
      The first one will be used by default.
      Please use the binpath argument to manually set the path to the bin
      folder.")
      binpath <- binpath[1]
    }
  } else {
    # Control
    controls(args = binpath, data = NULL, type = "character")
    if(!file.exists(binpath)) {
      stop(paste0("Impossible to access ", binpath), call. = FALSE)
    }
  }
  
  # Control infomap_version
  infomap_versiondispo <- c("2.1.0", "2.6.0")
  if(!is.character(infomap_version)) {
    stop("infomap_version must be a character", call. = FALSE)
  }
  
  infomap_version <- infomap_version[!duplicated(infomap_version)]
  
  if(length(infomap_version) > length(infomap_versiondispo)) {
    stop(paste0(
      "Please choose versions of Infomap in the list: ",
      paste(infomap_versiondispo, collapse = " ")
    ), call. = FALSE)
  }
  
  if (length(setdiff(infomap_version, infomap_versiondispo)) > 0) {
    stop(paste0(
      "Please choose versions of Infomap in the list: ",
      paste(infomap_versiondispo, collapse = " ")
    ), call. = FALSE)
  }
  nbversion <- length(infomap_version)
  
  # Check if bin.zip and bin already exists and remove them
  if (file.exists(paste0(binpath, "/bin.zip"))) {
    unlink(paste0(binpath, "/bin.zip"))
  }
  if (dir.exists(paste0(binpath, "/bin"))) {
    unlink(paste0(binpath, "/bin"), recursive = TRUE)
  }
  
  # Download bin.zip
  message(" ")
  message("1. Download bin.zip")
  message(" ")
  utils::download.file(
    "https://www.mmmycloud.com/index.php/s/DtZqrXAora6SzLo/download",
    paste0(binpath, "/bin.zip"),
    mode = "wb")
  
  # Unzip folder
  message(" ")
  message("2. Unzip folder")
  message(" ")
  utils::unzip(zipfile = paste0(binpath, "/bin.zip"), exdir = binpath)
  
  # Delete bin.zip
  unlink(paste0(binpath, "/bin.zip"))
  
  # Check presence files
  nboslom <- length(list.files(paste0(binpath, "/bin/OSLOM")))
  nbinfomap <- 0
  for(vinf in 1:nbversion) {
    nbinfomap <- nbinfomap +
      length(list.files(paste0(binpath, "/bin/INFOMAP/",
                               infomap_version[vinf])))
  }
  nblouvain <- length(list.files(paste0(binpath, "/bin/LOUVAIN")))
  
  if(nboslom == 8 & (nbinfomap == 8 * nbversion) & nblouvain == 8) {
    message(paste0(
      "The folder has been successfully downloaded and dezipped in ",
      binpath
    ))
  } else {
    unlink(paste0(binpath, "/bin"), recursive = TRUE)
    stop(paste0("An error occurred, download and/or dezip failed"),
         call. = FALSE)
  }
  
  # Identify OS
  os <- Sys.info()[["sysname"]]
  if (os == "Linux") {
    osid <- "lin"
  }
  if (os == "Windows") {
    osid <- "win"
  }
  if (os == "Darwin") {
    osid <- "mac"
  }
  
  # List files
  files <- NULL
  for (vinf in 1:nbversion) {
    files <- c(
      files,
      paste0(binpath, "/bin/INFOMAP/", infomap_version[vinf],
             "/infomap_omp_", osid),
      paste0(binpath, "/bin/INFOMAP/", infomap_version[vinf],
             "/infomap_noomp_", osid)
    )
  }
  files <- c(files,
             paste0(binpath, "/bin/LOUVAIN/convert_", osid),
             paste0(binpath, "/bin/LOUVAIN/louvain_", osid),
             paste0(binpath, "/bin/OSLOM/oslom_dir_", osid),
             paste0(binpath, "/bin/OSLOM/oslom_undir_", osid))
  
  if (osid == "win") {
    files <- paste0(files, ".exe")
  }
  nbfiles <- length(files)
  
  # Check permissions
  message(" ")
  message("3. Check permissions")
  message(" ")
  
  perm <- rep(-1, nbfiles)
  for(f in 1:nbfiles) {
    file <- files[f]
    perm[f] <- file.access(file, mode = 1)
    if (perm[f] == -1) {
      message(paste0("Permission to execute ", file, " as program: denied"))
    } else {
      perm[f] == 10
      message(paste0("Permission to execute ", file, " as program: granted"))
    }
  }
  
  if (sum(perm == -1) > 0) {
    message(" ")
    message("Try to change permissions automatically")
    message(" ")
    
    for(f in 1:nbfiles) {
      file <- files[f]
      if (perm[f] == -1) {
        if (osid == "lin") { # Linux
          system(paste0("chmod +x ", file))
        }
        if (osid == "win") { # Windows
          # system(paste0("chmod +x ", file))
        }
        if (osid == "mac") { # Mac
          system(paste0("chmod 755 ", file))
        }
        perm[f] <- file.access(file, mode = 1)
        if (perm[f] == -1) {
          message(paste0("Automatic change of permission of ", file,
                         " failed"))
        } else {
          perm[f] == 10
          message(paste0(
            "Automatic change of permission succeed, ", file,
            " can now be executed as progam"
          ))
        }
      }
    }
  }
  
  if (sum(perm == -1) > 0) {
    maxtry <- 10000
    nbtry <- 0
    while (nbtry < maxtry) {
      if (sum(perm == -1) > 0) {
        nopermfiles <- files[perm == -1]
      } else {
        message(" ")
        message("All permissions granted!")
        break
      }
      
      message(" ")
      message("Permission to execute the following files as program denied")
      message(" ")
      for (file in nopermfiles) {
        message(file)
      }
      
      message(" ")
      ask <- utils::menu(c(
        "I've just tried to change the permission manually and I want to
        check the permission again",
        "I would like to continue the execution of the function without
        checking the permission",
        "I would like to stop the function"
      ),
      title = "You can now try to change the permission of the above files
      manually and check the permission again")
      
      if(ask == 1) {
        for (f in 1:nbfiles) {
          file <- files[f]
          perm[f] <- file.access(file, mode = 1)
          if (perm[f] == -1) {
            #
          } else {
            perm[f] == 10
          }
        }
      } else if (ask == 2) {
        break
      } else {
        unlink(paste0(binpath, "/bin"), recursive = TRUE)
        message(" ")
        stop("Function install_binaries() stopped", call. = FALSE)
      }
    }
  }
  
  # Test INFOMAP
  message(" ")
  message("4. Test Infomap")
  message(" ")
  
  for(vinf in 1:nbversion) {
    version <- infomap_version[vinf]
    path <- paste0(binpath, "/bin/INFOMAP/", version, "/")
    files <- c(paste0("infomap_omp_", osid), paste0("infomap_noomp_", osid))
    if (osid == "lin") {
      cmd <- paste0(
        path, files[1], " -N 10 --two-level --tree --markov-time 0.5 ",
        path, "example.txt ", path
      )
      cmd <- paste0(cmd, " >/dev/null 2>&1")
      system(cmd)
    }
    if (osid == "mac") {
      cmd <- paste0(
        path, files[1], " -N 10 --two-level --tree --markov-time 0.5 ",
        path, "example.txt ", path
      )
      cmd <- paste0(cmd, " >/dev/null 2>&1")
      system(cmd)
    }
    if (osid == "win") {
      files <- paste0(files, ".exe")
      cmd <- paste0(
        path, files[1], " -N 10 --two-level --tree --markov-time 0.5 ",
        path, "example.txt ", path
      )
      system(cmd, show.output.on.console = FALSE)
    }
    testopm <- TRUE
    if (!("example.tree" %in% list.files(path))) {
      testopm <- FALSE
    }
    if (file.exists(paste0(path, "/example.tree"))) {
      unlink(paste0(path, "/example.tree"))
    }
    
    if (osid == "lin") {
      cmd <- paste0(
        path, files[2], " -N 10 --two-level --tree --markov-time 0.5 ",
        path, "example.txt ", path
      )
      cmd <- paste0(cmd, " >/dev/null 2>&1")
      system(cmd)
    }
    if (osid == "mac") {
      cmd <- paste0(
        path, files[2], " -N 10 --two-level --tree --markov-time 0.5 ",
        path, "example.txt ", path
      )
      cmd <- paste0(cmd, " >/dev/null 2>&1")
      system(cmd)
    }
    if (osid == "win") {
      cmd <- paste0(
        path, files[2], " -N 10 --two-level --tree --markov-time 0.5 ",
        path, "example.txt ", path
      )
      system(cmd, show.output.on.console = FALSE)
    }
    testnoopm <- TRUE
    if (!("example.tree" %in% list.files(path))) {
      testnoopm <- FALSE
    }
    if (file.exists(paste0(path, "/example.tree"))) {
      unlink(paste0(path, "/example.tree"))
    }
    
    if (!(testopm | testnoopm)) {
      message(" ")
      message("Infomap is not installed...")
      message("Please have a look at
      https//biorgeo.github.io/bioRgeo/articles/a3_1_install_executable_binary_files.html
              for more details")
    } else {
      if (testopm) {
        message(
          "Congratulation, you successfully install the ",
          version, " OpenMP version of Infomap!"
        )
        file.copy(
          paste0(path, files[1]),
          paste0(path, "infomap_", substr(files[1], 13, nchar(file[1])))
        )
      } else {
        message(" ")
        message("Congratulation, you successfully install the ", version,
                " no OpenMP version of Infomap!"
        )
        file.copy(
          paste0(path, files[2]),
          paste0(path, "infomap_", substr(files[1], 13, nchar(file[1])))
        )
        message(" ")
        message("A library is probably missing to install the OpenMP
                version...")
        message("Please have a look at 
https//biorgeo.github.io/bioRgeo/articles/a3_1_install_executable_binary_files.html
                for more details")
      }
      utils::write.table(1, paste0(path, "check.txt"))
    }
  }
  
  # Test LOUVAIN
  message(" ")
  message("5. Test Louvain")
  message(" ")
  
  path <- paste0(binpath, "/bin/LOUVAIN/")
  version <- list.files(path)[substr(list.files(path), 1, 7) == "version"]
  version <- substr(version, 9, nchar(version))
  files <- c(paste0("convert_", osid), paste0("louvain_", osid))
  if (osid == "lin") {
    cmd <- paste0(
      path, files[1], " -i ", path, "example.txt -o ", path,
      "example.bin"
    )
    cmd <- paste0(cmd, " >/dev/null 2>&1")
    system(cmd)
  }
  if (osid == "mac") {
    cmd <- paste0(
      path, files[1], " -i ", path, "example.txt -o ", path,
      "example.bin"
    )
    cmd <- paste0(cmd, " >/dev/null 2>&1")
    system(cmd)
  }
  if (osid == "win") {
    files <- paste0(files, ".exe")
    cmd <- paste0(
      path, files[1], " -i ", path, "example.txt -o ", path,
      "example.bin"
    )
    system(cmd, show.output.on.console = FALSE)
  }
  testconvert <- TRUE
  if (!("example.bin" %in% list.files(path))) {
    testconvert <- FALSE
  }
  
  if (testconvert) {
    cmd <- paste0(path, files[2], " ", path, "example.bin -l -1 -q id_qual")
    tree <- system(cmd, intern = TRUE)
    testlouvain <- TRUE
    if (tree[1] != "0 0") {
      testlouvain <- FALSE
    }
  }
  
  if (file.exists(paste0(path, "/example.bin"))) {
    unlink(paste0(path, "/example.bin"))
  }
  
  if (!testlouvain) {
    message(" ")
    message("Louvain is not installed...")
    message("Please have a look at 
https//biorgeo.github.io/bioRgeo/articles/a3_1_install_executable_binary_files.html
            for more details")
  } else {
    message(
      "Congratulation, you successfully install the version ",
      version, " of Louvain!"
    )
    utils::write.table(1, paste0(path, "check.txt"))
  }
  
  # Test OSLOM
  message(" ")
  message("6. Test OSLOM")
  message(" ")
  
  path <- paste0(binpath, "/bin/OSLOM/")
  version <- list.files(path)[substr(list.files(path), 1, 7) == "version"]
  version <- substr(version, 9, nchar(version))
  files <- c(paste0("oslom_undir_", osid), paste0("oslom_dir_", osid))
  if (osid == "lin") {
    cmd <- paste0(path, files[1], " -f ", path, "example.txt -uw")
    cmd <- paste0(cmd, " >/dev/null 2>&1")
    system(cmd)
  }
  if (osid == "mac") {
    cmd1 <- paste0("cd ", path, " >/dev/null 2>&1")
    cmd2 <- paste0("./oslom_undir_mac -f example.txt -uw > /dev/null 2>&1")
    cmd <- paste0(cmd1, " && ", cmd2)
    system(cmd)
  }
  if (osid == "win") {
    files <- paste0(files, ".exe")
    cmd <- paste0(path, files[1], " -f ", path, "example.txt -uw")
    dir.create(paste0(path, "example.txt_oslo_files"),
               showWarnings = FALSE,
               recursive = TRUE
    )
    system(cmd, show.output.on.console = FALSE)
  }
  
  testundir <- TRUE
  if (!("tp" %in% list.files(paste0(path, "example.txt_oslo_files")))) {
    testundir <- FALSE
  }
  if (dir.exists(paste0(path, "example.txt_oslo_files"))) {
    unlink(paste0(path, "example.txt_oslo_files"), recursive = TRUE)
  }
  if (file.exists(paste0(path, "tp"))) {
    unlink(paste0(path, "tp"))
  }
  if (file.exists("tp")) {
    unlink("tp")
  }
  if (file.exists(paste0(path, "time_seed.dat"))) {
    unlink(paste0(path, "time_seed.dat"))
  }
  if (file.exists("time_seed.dat")) {
    unlink("time_seed.dat")
  }
  
  if (osid == "lin") {
    cmd <- paste0(path, files[2], " -f ", path, "example.txt -uw")
    cmd <- paste0(cmd, " >/dev/null 2>&1")
    system(cmd)
  }
  if (osid == "mac") {
    cmd1 <- paste0("cd ", path, " >/dev/null 2>&1")
    cmd2 <- paste0("./oslom_undir_mac -f example.txt -uw > /dev/null 2>&1")
    cmd <- paste0(cmd1, " && ", cmd2)
    system(cmd)
  }
  if (osid == "win") {
    cmd <- paste0(path, files[2], " -f ", path, "example.txt -uw")
    dir.create(paste0(path, "example.txt_oslo_files"),
               showWarnings = FALSE,
               recursive = TRUE
    )
    system(cmd, show.output.on.console = FALSE)
  }
  testdir <- TRUE
  if (!("tp" %in% list.files(paste0(path, "example.txt_oslo_files")))) {
    testdir <- FALSE
  }
  if (dir.exists(paste0(path, "example.txt_oslo_files"))) {
    unlink(paste0(path, "example.txt_oslo_files"), recursive = TRUE)
  }
  if (file.exists(paste0(path, "tp"))) {
    unlink(paste0(path, "tp"))
  }
  if (file.exists("tp")) {
    unlink("tp")
  }
  if (file.exists(paste0(path, "time_seed.dat"))) {
    unlink(paste0(path, "time_seed.dat"))
  }
  if (file.exists("time_seed.dat")) {
    unlink("time_seed.dat")
  }
  
  if (!testundir) {
    message(" ")
    message("OSLOM is not installed...")
    message("Please have a look at 
https//biorgeo.github.io/bioRgeo/articles/a3_1_install_executable_binary_files.html
            for more details")
  } else {
    message("Congratulation, you successfully install the version ", version,
            " of OSLOM!")
    utils::write.table(1, paste0(path, "check.txt"))
    
    if (!testdir) {
      message("Warning: only the undirected version of OSLOM has been
              installed...")
      message("Please have a look at 
https//biorgeo.github.io/bioRgeo/articles/a3_1_install_executable_binary_files.html
              for more details")
    } else {
      utils::write.table(1, paste0(path, "checkdir.txt"))
    }
  }
  
  # Remove unnecessary files in INFOMAP
  for (vinf in 1:length(infomap_versiondispo)) {
    version <- infomap_versiondispo[vinf]
    if (file.exists(paste0(binpath, "/bin/INFOMAP/", version, "/check.txt"))){
      unlink(paste0(binpath, "/bin/INFOMAP/", version, "/infomap_noomp_mac"))
      unlink(paste0(binpath, "/bin/INFOMAP/", version, "/infomap_omp_mac"))
      unlink(paste0(binpath, "/bin/INFOMAP/", version,
                    "/infomap_noomp_win.exe"))
      unlink(paste0(binpath, "/bin/INFOMAP/", version, "/infomap_omp_win.exe"))
      unlink(paste0(binpath, "/bin/INFOMAP/", version, "/infomap_noomp_lin"))
      unlink(paste0(binpath, "/bin/INFOMAP/", version, "/infomap_omp_lin"))
      unlink(paste0(binpath, "/bin/INFOMAP/", version, "/example.txt"))
    } else {
      unlink(paste0(binpath, "/bin/INFOMAP/", version), recursive = TRUE)
    }
  }
  if (length(list.files(paste0(binpath, "/bin/INFOMAP/"))) == 0) {
    unlink(paste0(binpath, "/bin/INFOMAP"), recursive = TRUE)
  }
  
  # Remove unnecessary files in LOUVAIN
  if (file.exists(paste0(binpath, "/bin/LOUVAIN/check.txt"))) {
    if (osid == "lin") {
      unlink(paste0(binpath, "/bin/LOUVAIN/convert_mac"))
      unlink(paste0(binpath, "/bin/LOUVAIN/louvain_mac"))
      unlink(paste0(binpath, "/bin/LOUVAIN/convert_win.exe"))
      unlink(paste0(binpath, "/bin/LOUVAIN/louvain_win.exe"))
    }
    if (osid == "mac") {
      unlink(paste0(binpath, "/bin/LOUVAIN/convert_lin"))
      unlink(paste0(binpath, "/bin/LOUVAIN/louvain_lin"))
      unlink(paste0(binpath, "/bin/LOUVAIN/convert_win.exe"))
      unlink(paste0(binpath, "/bin/LOUVAIN/louvain_win.exe"))
    }
    if (osid == "win") {
      unlink(paste0(binpath, "/bin/LOUVAIN/convert_mac"))
      unlink(paste0(binpath, "/bin/LOUVAIN/louvain_mac"))
      unlink(paste0(binpath, "/bin/LOUVAIN/convert_lin"))
      unlink(paste0(binpath, "/bin/LOUVAIN/louvain_lin"))
    }
  } else {
    unlink(paste0(binpath, "/bin/LOUVAIN"), recursive = TRUE)
  }
  
  # Remove unnecessary files in OSLOM
  if (file.exists(paste0(binpath, "/bin/OSLOM/check.txt"))) {
    if (osid == "lin") {
      unlink(paste0(binpath, "/bin/OSLOM/oslom_dir_mac"))
      unlink(paste0(binpath, "/bin/OSLOM/oslom_undir_mac"))
      unlink(paste0(binpath, "/bin/OSLOM/oslom_dir_win.exe"))
      unlink(paste0(binpath, "/bin/OSLOM/oslom_undir_win.exe"))
    }
    if (osid == "mac") {
      unlink(paste0(binpath, "/bin/OSLOM/oslom_dir_lin"))
      unlink(paste0(binpath, "/bin/OSLOM/oslom_undir_lin"))
      unlink(paste0(binpath, "/bin/OSLOM/oslom_dir_win.exe"))
      unlink(paste0(binpath, "/bin/OSLOM/oslom_undir_win.exe"))
    }
    if (osid == "win") {
      unlink(paste0(binpath, "/bin/OSLOM/oslom_dir_mac"))
      unlink(paste0(binpath, "/bin/OSLOM/oslom_undir_mac"))
      unlink(paste0(binpath, "/bin/OSLOM/oslom_dir_lin"))
      unlink(paste0(binpath, "/bin/OSLOM/oslom_undir_lin"))
    }
  } else {
    unlink(paste0(binpath, "/bin/OSLOM"), recursive = TRUE)
  }
}
