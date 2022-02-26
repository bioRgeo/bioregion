#' Download, unzip and check permission of bioRgeo's executables
#'
#' This function download and unzip the bin folder needed to run some bioRgeo's functions. It also checks if the files
#' have the permissions to be executed as programs.
#'
#' @param binpath a string indicating the path to the folder that will host the bin folder (bioRgeo's package by default,
#' if you use a different folder please be sure to update the \code{binpath} in \link{infomap}, \link{louvain}
#' and \link{oslom})
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}),
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) and
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' @export
install_binaries <- function(binpath = NULL) {

  # Set binpath
  if (is.null(binpath)) {
    # Identify bioRgeo directory on your computer
    biodir <- list.dirs(.libPaths(), recursive = FALSE)
    binpath <- biodir[grep("bioRgeo", biodir)]
  } else {
    # Control
    if (!is.character(binpath)) {
      stop("path must be a string")
    }
    if (!file.exists(binpath)) {
      stop(paste0("Impossible to access ", binpath))
    }
  }

  # Check if bin.zip and bin already exists and remove them
  if (file.exists(paste0(binpath, "/bin.zip"))) {
    unlink(paste0(binpath, "/bin.zip"))
  }
  if (dir.exists(paste0(binpath, "/bin"))) {
    unlink(paste0(binpath, "/bin"), recursive = TRUE)
  }
  if (dir.exists(paste0(binpath, "/bin2"))) { # TO REMOVE AFTER TEST
    unlink(paste0(binpath, "/bin2"), recursive = TRUE)
  }

  # Download file
  message(" ")
  message("1. Download bin.zip")
  utils::download.file("https://www.mmmycloud.com/index.php/s/DtZqrXAora6SzLo/download", paste0(binpath, "/bin.zip"), mode = "wb")

  # Unzip file
  message(" ")
  message("2. Unzip file")
  utils::unzip(zipfile = paste0(binpath, "/bin.zip"), exdir = binpath)

  # Delete bin.zip
  unlink(paste0(binpath, "/bin.zip"))

  # Check presence files
  nboslom <- length(list.files(paste0(binpath, "/bin/OSLOM")))
  nbinfomap <- length(list.files(paste0(binpath, "/bin/INFOMAP")))
  nblouvain <- length(list.files(paste0(binpath, "/bin/LOUVAIN")))

  if (nboslom == 8 & nbinfomap == 8 & nblouvain == 8) {
    message(paste0("The folder has been successfully downloaded and dezipped in ", binpath))
  } else {
    unlink(paste0(binpath, "/bin"), recursive = TRUE)
    stop(paste0("An error occurred, download and/or dezip failed"), call.=FALSE)
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
  files <- c(
    paste0(binpath, "/bin/INFOMAP/infomap_opm_", osid), paste0(binpath, "/bin/INFOMAP/infomap_noopm_", osid),
    paste0(binpath, "/bin/LOUVAIN/convert_", osid), paste0(binpath, "/bin/LOUVAIN/louvain_", osid),
    paste0(binpath, "/bin/OSLOM/oslom_dir_", osid), paste0(binpath, "/bin/OSLOM/oslom_undir_", osid)
  )
  if (osid == "win") {
    files <- paste0(files, ".exe")
  }
  nbfiles=length(files)

  # Check permissions
  message(" ")
  message("3. Check permissions")
  message(" ")

  perm <- rep(-1,nbfiles)
  for (f in 1:nbfiles) {
    file=files[f]
    perm[f] <- file.access(file, mode = 1)
    if (perm[f] == -1) {
      message(paste0("Permission to execute ", file, " as program: denied"))
    } else {
      perm[f] == 10
      message(paste0("Permission to execute ", file, " as program: granted"))
    }
  }

  if(sum(perm == -1)>0){

    message(" ")
    message("Try to change permissions automatically")
    message(" ")

    for (f in 1:nbfiles) {
      file=files[f]
      if(perm[f] == -1){
        if (osid == "lin") { # Linux
          system(paste0("chmod +x ", file))
        }
        if (osid == "win") { # Windows
          #system(paste0("chmod +x ", file))
        }
        if (osid == "mac") { # Mac
          #system(paste0("chmod +x ", file))
        }
        perm[f] <- file.access(file, mode = 1)
        if (perm[f] == -1) {
          message(paste0("Automatic change of permission of ", file, " failed"))
        } else {
          perm[f] == 10
          message(paste0("Automatic change of permission succeed, ", file, " can now be executed as progam"))
        }
      }
    }
  }

  if(sum(perm == -1)>0){

    maxtry=10000
    nbtry <- 0
    while (nbtry < maxtry){

      if(sum(perm == -1)>0){
        nopermfiles=files[perm==-1]
      }else{
        message(" ")
        message("All permissions granted!")
        break
      }

      message(" ")
      message("Permission to execute the following files as program denied")
      message(" ")
      for(file in nopermfiles){
        message(file)
      }

      message(" ")
      ask <- utils::menu(c("I've just tried to change the permission manually and I want to check the permissions again",
                    "I would like to continue the execution of the function without checking the permissions",
                    "I would like to stop the function"),
                  title = "You can now try to change the permission of the above files manually and check the permission again")

      if(ask == 1){
        for (f in 1:nbfiles) {
          file=files[f]
          perm[f] <- file.access(file, mode = 1)
          if (perm[f] == -1) {
            #
          } else {
            perm[f] == 10
          }
        }

      }else if(ask == 2){
        break
      }else{
        unlink(paste0(binpath, "/bin"), recursive = TRUE)
        message(" ")
        stop("Function install_binaries() stopped", call.=FALSE)
      }
    }
  }

  # Test INFOMAP
  path=paste0(binpath, "/bin/INFOMAP/")
  version=list.files(path)[substr(list.files(path),1,7)=="version"]
  version=substr(version, 9, nchar(version))
  files <- c(paste0("infomap_opm_", osid), paste0("infomap_noopm_", osid))
  if (osid == "win") {
    files <- paste0(files, ".exe")
  }

  cmd=paste0(path, files[1], " -N 10 --two-level --tree --markov-time 0.5 ", path, "example.txt ", path)
  if(osid == "lin"){
    cmd=paste0(cmd, " >/dev/null")
  }
  system(cmd)
  testopm=TRUE
  if(!("example.tree" %in% list.files(path))){
    testopm=FALSE
  }
  if (file.exists(paste0(path,"/example.tree"))) {
    unlink(paste0(path,"/example.tree"))
  }

  cmd=paste0(path, files[2], " -N 10 --two-level --tree --markov-time 0.5 ", path, "example.txt ", path)
  if(osid == "lin"){
    cmd=paste0(cmd, " >/dev/null")
  }
  system(cmd)
  testnoopm=TRUE
  if(!("example.tree" %in% list.files(path))){
    testnoopm=FALSE
  }
  if (file.exists(paste0(path,"/example.tree"))) {
    unlink(paste0(path,"/example.tree"))
  }

  if(!(testopm | testnoopm)){
    message(" ")
    message("Infomap is not running!")
  }else{
    if(testopm){
      message(" ")
      message("Congratulation, you successfully install the ", version, " OpenMP version of Infomap!")
      file.copy(paste0(path, files[1]), paste0(path, "infomap_", substr(files[1],13,nchar(file[1]))))
    }else{
      message(" ")
      message("Congratulation, you successfully install the ", version, " no OpenMP version of Infomap!")
      file.copy(paste0(path, files[2]), paste0(path, "infomap_", substr(files[1],13,nchar(file[1]))))
    }
    utils::write.table(1, paste0(path, "check.txt"))
  }

  # Test LOUVAIN
  path=paste0(binpath, "/bin/LOUVAIN/")
  version=list.files(path)[substr(list.files(path),1,7)=="version"]
  version=substr(version, 9, nchar(version))
  files <- c(paste0("convert_", osid), paste0("louvain_", osid))
  if (osid == "win") {
    files <- paste0(files, ".exe")
  }

  cmd=paste0(path, files[1], " -i ", path, "example.txt -o ", path, "example.bin")
  if(osid == "lin"){
    cmd=paste0(cmd, " >/dev/null")
  }
  system(cmd)
  testconvert=TRUE
  if(!("example.bin" %in% list.files(path))){
    testconvert=FALSE
  }

  if(testconvert){
    cmd=paste0(path, files[2], " ", path, "example.bin -l -1 -q id_qual")
    tree=system(cmd, intern=TRUE)
    testlouvain=TRUE
    if(tree[1]!="0 0"){
      testlouvain=FALSE
    }
  }

  if (file.exists(paste0(path,"/example.bin"))) {
    unlink(paste0(path,"/example.bin"))
  }

  if(!testlouvain){
    message(" ")
    message("Louvain is not running!")
  }else{
    message(" ")
    message("Congratulation, you successfully install the version ", version, " of Louvain!")
    utils::write.table(1, paste0(path, "check.txt"))
  }

  # Test OSLOM
  path=paste0(binpath, "/bin/OSLOM/")
  version=list.files(path)[substr(list.files(path),1,7)=="version"]
  version=substr(version, 9, nchar(version))
  files <- c(paste0("oslom_undir_", osid), paste0("oslom_dir_", osid))
  if (osid == "win") {
    files <- paste0(files, ".exe")
  }

  cmd=paste0(path, files[1], " -f ", path, "example.txt -uw")
  if(osid == "lin"){
    cmd=paste0(cmd, " >/dev/null 2>&1")
  }
  system(cmd)
  testundir=TRUE
  if(!("tp" %in% list.files(paste0(path, "example.txt_oslo_files")))){
    testundir=FALSE
  }
  if (dir.exists(paste0(path,"example.txt_oslo_files"))) {
    unlink(paste0(path,"example.txt_oslo_files"), recursive=TRUE)
  }

  cmd=paste0(path, files[2], " -f ", path, "example.txt -uw")
  if(osid == "lin"){
    cmd=paste0(cmd, " >/dev/null 2>&1")
  }
  system(cmd)
  testdir=TRUE
  if(!("tp" %in% list.files(paste0(path, "example.txt_oslo_files")))){
    testdir=FALSE
  }
  if (dir.exists(paste0(path,"example.txt_oslo_files"))) {
    unlink(paste0(path,"example.txt_oslo_files"), recursive=TRUE)
  }

  if(!testundir){
    message(" ")
    message("OSLOM is not running!")
  }else{
    message(" ")
    message("Congratulation, you successfully install the version ", version, " of OSLOM!")
    utils::write.table(1, paste0(path, "check.txt"))
    if(!testdir){
      message("Warning: only the undirected version of OSLOM has been install!")
    }
  }

  # Remove unnecessary files


}
