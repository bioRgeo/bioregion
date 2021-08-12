#' Download, unzip and check permission of bioRgeo's executables
#'
#' This function download and unzip the bin folder needed to run some bioRgeo's functions. It also check if the files
#' have the permissions to be executed as program.
#'
#' @param binpath a string indicating the path to the folder that will host the bin folder (bioRgeo's package by default,
#' if you use a different folder please be sure to update the \code{binpath} in \link{infomap}, \link{louvain}
#' and \link{oslom})
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}),
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) and
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' @export
bin <- function(binpath = NULL){

  if(is.null(binpath)){
    # Identify bioRgeo directory on your computer
    biodir <- list.dirs(.libPaths(), recursive = FALSE)
    binpath <- biodir[grep("bioRgeo", biodir)]
  }else{
    # Control
    if(!is.character(binpath)){
      stop("path must be a string")
    }
    if(!file.exists(binpath)){
      stop(paste0("Impossible to access ", binpath))
    }
  }

  # Download file
  utils::download.file("https://www.mmmycloud.com/index.php/s/diAFHLajspZDjYx/download", paste0(binpath,"/bin.zip"))

  # Unzip file
  utils::unzip(zipfile=paste0(binpath,"/bin.zip"), exdir=binpath)

  # Delete bin.zip
  unlink(paste0(binpath,"/bin.zip"))

  # Check presence files
  nboslom=length(list.files(paste0(binpath,"/bin/OSLOM")))
  nbinfomap=length(list.files(paste0(binpath,"/bin/INFOMAP")))
  nblouvain=length(list.files(paste0(binpath,"/bin/LOUVAIN")))

  if(nboslom==4 & nbinfomap==2 & nblouvain==4){
    print(paste0("The folder has been successfully downloaded and dezipped in ", binpath))
  }else{
    print(paste0("An error occurred, download and/or dezipp failed"))
  }

  # Check OS
  os=Sys.info()[['sysname']]

  # Check permissions
  if(os == "Linux"){

    # INFOMAP/infomap_lin
    file="INFOMAP/infomap_lin"
    perm=file.access(paste0(binpath, "/bin/", file), mode = 1)
    if(perm == -1){
      system(paste0("chmod +x ", binpath, "/bin/", file))
    }
    perm=file.access(paste0(binpath, "/bin/", file), mode = 1)
    if(perm == -1){
      print(paste0("Automatic change of permission failed, please allow executing ", file, " as program manually"))
    }else{
      print(paste0("Automatic change of permission successed, ", file, " can now be executed as progam"))
    }

    # OSLOM/oslom_dir_lin
    file="OSLOM/oslom_dir_lin"
    perm=file.access(paste0(binpath, "/bin/", file), mode = 1)
    if(perm == -1){
      system(paste0("chmod +x ", binpath, "/bin/", file))
    }
    perm=file.access(paste0(binpath, "/bin/", file), mode = 1)
    if(perm == -1){
      print(paste0("Automatic change of permission failed, please allow executing ", file, " as program manually"))
    }else{
      print(paste0("Automatic change of permission successed, ", file, " can now be executed as progam"))
    }

    # OSLOM/oslom_undir_lin
    file="OSLOM/oslom_undir_lin"
    perm=file.access(paste0(binpath, "/bin/", file), mode = 1)
    if(perm == -1){
      system(paste0("chmod +x ", binpath, "/bin/", file))
    }
    perm=file.access(paste0(binpath, "/bin/", file), mode = 1)
    if(perm == -1){
      print(paste0("Automatic change of permission failed, please allow executing ", file, " as program manually"))
    }else{
      print(paste0("Automatic change of permission successed, ", file, " can now be executed as progam"))
    }

    # LOUVAIN/convert_lin
    file="LOUVAIN/convert_lin"
    perm=file.access(paste0(binpath, "/bin/", file), mode = 1)
    if(perm == -1){
      system(paste0("chmod +x ", binpath, "/bin/", file))
    }
    perm=file.access(paste0(binpath, "/bin/", file), mode = 1)
    if(perm == -1){
      print(paste0("Automatic change of permission failed, please allow executing ", file, " as program manually"))
    }else{
      print(paste0("Automatic change of permission successed, ", file, " can now be executed as progam"))
    }

    # LOUVAIN/louvain_lin
    file="LOUVAIN/louvain_lin"
    perm=file.access(paste0(binpath, "/bin/", file), mode = 1)
    if(perm == -1){
      system(paste0("chmod +x ", binpath, "/bin/", file))
    }
    perm=file.access(paste0(binpath, "/bin/", file), mode = 1)
    if(perm == -1){
      print(paste0("Automatic change of permission failed, please allow executing ", file, " as program manually"))
    }else{
      print(paste0("Automatic change of permission successed, ", file, " can now be executed as progam"))
    }

  }


}
