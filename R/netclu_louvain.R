#' Louvain community finding
#'
#' This function finds communities in a (un)weighted undirected network based on the Louvain algorithm.
#'
#' @param net a two- or three-column \code{data.frame} representing a network with the two first columns
#' as (un)directed links between pair of nodes and an optional third column indicating the weight of the link
#' @param weight a boolean indicating if the weights should be considered if there is a third column
#' @param lang a string indicating what version of Louvain should be used (igraph or Cpp, see Details)
#' @param q the quality function used to compute partition of the graph (modularity is chosen by default, see Details)
#' @param c the parameter for the Owsinski-Zadrozny quality function (between 0 and 1, 0.5 is chosen by default)
#' @param k the kappa_min value for the Shi-Malik quality function (it must be > 0, 1 is chosen by default)
#' @param delete_temp a boolean indicating if the temporary folder should be removed (see Details)
#' @param path_temp a string indicating the path to the temporary folder (see Details)
#' @param binpath a string indicating the path to the bin folder (see \link{install_binaries} and Details)
#' @export
#' @details
#' Louvain is a network community detection algorithm proposed in \insertCite{Blondel2008}{bioRgeo}. This function
#' proposed two implementations of the function (parameter \code{lang}): the \href{https://cran.r-project.org/web/packages/igraph/index.html}{igraph} implementation
#' (\link[igraph]{cluster_louvain}) and the C++ implementation (\url{https://sourceforge.net/projects/louvain/}, version 0.3).
#' The latest offers the possibility to choose among several quality functions, \code{q = 0} for the classical
#' Newman-Girvan criterion (also called "Modularity"), 1 for the Zahn-Condorcet criterion, 2 for
#' the Owsinski-Zadrozny criterion (you should specify the value of the parameter with option -c), 3	for the Goldberg
#' Density criterion, 4	for the A-weighted Condorcet criterion, 5 for the Deviation to Indetermination criterion, 6
#' for the Deviation to Uniformity criterion, 7 for the Profile Difference criterion, 8	for the Shi-Malik criterion
#' (you should specify the value of kappa_min with option -k) and 9	for the Balanced Modularity criterion.
#'
#' The C++ version of Louvain is based on the version 0.3 (\url{https://sourceforge.net/projects/louvain/}).
#' This function needs executable files to run. They can be installed with \link{install_binaries}. If you set the path to
#' the folder that will host the bin folder manually while running \link{install_binaries} please make sure to set \code{binpath}
#' accordingly.
#'
#' The C++ version of Louvain generates temporary folders and/or files that are stored in the \code{path_temp} folder
#' (folder "louvain_temp" in the working directory by default). This temporary folder is removed by default
#' (\code{delete_temp = TRUE}).
#'
#' @return A \code{data.frame} providing one community by node.
#'
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}),
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) and
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' @seealso \link{netclu_infomap}, \link{netclu_oslom}
#' @examples
#' comat=matrix(sample(1000,50),5,10)
#' rownames(comat)=paste0("Site",1:5)
#' colnames(comat)=paste0("Species",1:10)
#'
#' net=similarity(comat,metric="Simpson")
#' com=netclu_louvain(net, lang="igraph")
#' @references
#' \insertRef{Blondel2008}{bioRgeo}
#' @export
netclu_louvain <- function(net, weight = TRUE, lang="Cpp", q = 0, c = 0.5, k = 1, delete_temp = TRUE,
                    path_temp = "louvain_temp", binpath = NULL){

  # Control input net
  if (!is.data.frame(net)) {
    stop("net must be a two- or three-columns data.frame")
  }

  if (dim(net)[2] != 2 & dim(net)[2] != 3) {
    stop("net must be a two- or three-columns data.frame")
  }

  nbna <- sum(is.na(net))
  if (nbna > 0) {
    stop("NA(s) detected in the data.frame")
  }

  # Control parameters
  if (!is.logical(weight)) {
    stop("weight must be a boolean")
  }

  if (weight & dim(net)[2] == 2) {
    stop("net must be a three-columns data.frame if weight equal TRUE")
  }

  if (weight & dim(net)[2] == 3) {
    if(!is.numeric(net[,3])){
      stop("The third column of net must be numeric")
    }else{
      minet <- min(net[, 3])
      if (minet < 0) {
        stop("The third column of net should contains only positive real: negative value detected!")
      }
    }
  }

  if(!(lang %in% c("Cpp","igraph","all"))){
    stop("The Louvain version is not available.
     Please chose among the followings:
         Cpp, igraph, all")
  }

  # Prepare input for LOUVAIN
  idnode1=as.character(net[,1])
  idnode2=as.character(net[,2])
  idnode=c(idnode1,idnode2)
  idnode=idnode[!duplicated(idnode)]
  idnode=data.frame(ID=1:length(idnode),ID_NODE=idnode)

  netemp=data.frame(node1=idnode[match(net[,1],idnode[,2]),1],node2=idnode[match(net[,2],idnode[,2]),1])
  if(weight){
    netemp=cbind(netemp,net[,3])
    netemp=netemp[netemp[,3]>0,]
    colnames(netemp)[3]="weight"
  }

  if(lang=="igraph" | lang=="all"){

    net=igraph::graph_from_data_frame(netemp, directed=FALSE)
    comtemp=igraph::cluster_louvain(net)
    comtemp=cbind(as.numeric(comtemp$names),as.numeric(comtemp$membership))

    com=data.frame(ID=idnode[,2], Com=0)
    com[match(comtemp[,1],idnode[,1]),2]=comtemp[,2]

    comr=com

  }

  if(lang=="Cpp" | lang=="all"){

    # Set binpath
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

    # Check OS
    os=Sys.info()[['sysname']]

    # Check if LOUVAIN has successfully been installed
    if (!file.exists(paste0(binpath, "/bin/LOUVAIN/check.txt"))) {
      stop("Louvain is not installed... Please have a look at https//biorgeo.github.io/bioRgeo/articles/bin.html for more details.")
    }

    # Create temp folder
    dir.create(path_temp, showWarnings = FALSE, recursive = TRUE)
    if(!file.exists(path_temp)){
      stop(paste0("Impossible to create directory ", path_temp))
    }

    # Control parameters
    if(!is.numeric(q)){
      stop("q must be numeric")
    }else{
      if(q<0){
        stop("q must be positive")
      }
      if((q-floor(q))>0){
        stop("q must be an integer higher or equal to 0")
      }
    }

    if(!is.numeric(c)){
      stop("c must be numeric")
    }else{
      if(c<0 | c>1){
        stop("c must be in the interval (0,1)")
      }
    }

    if(!is.numeric(k)){
      stop("k must be numeric")
    }else{
      if(k<0){
        stop("k must be positive")
      }
    }

    if(!is.logical(delete_temp)){
      stop("delete_temp must be a boolean")
    }

    if(!is.character(path_temp)){
      stop("path_temp must be a string")
    }

    # Export input in LOUVAIN folder
    utils::write.table(netemp, paste0(path_temp, "/net.txt"), row.names=FALSE, col.names=FALSE, sep=" ")

    # Prepare command to run LOUVAIN
    # Convert net.txt with LOUVAIN
    if(weight){
      cmd=paste0("-i ", path_temp, "/net.txt -o ", path_temp, "/net.bin -w ", path_temp, "/net.weights")
    }else{
      cmd=paste0("-i ", path_temp, "/net.txt -o ", path_temp, "/net.bin")
    }

    if(os == "Linux"){
      cmd=paste0(binpath, "/bin/LOUVAIN/convert_lin ", cmd)
    }else if(os == "Windows"){
      cmd=paste0(binpath, "/bin/LOUVAIN/convert_win.exe ", cmd)
    }else if(os == "Darwin"){
      cmd=paste0(binpath, "/bin/LOUVAIN/convert_mac ", cmd)
    }else{
      stop("Linux, Windows or Mac distributions only.")
    }

    tree=system(command = cmd)

    # Run LOUVAIN
    if(weight){
      cmd=paste0(path_temp, "/net.bin -l -1 -q ", q, " -c ", c , " -k ", k," -w ", path_temp, "/net.weights")
    }else{
      cmd=paste0(path_temp, "/net.bin -l -1 -q ", q, " -c ", c , " -k ", k)
    }

    if(os == "Linux"){
      cmd=paste0(binpath, "/bin/LOUVAIN/louvain_lin ", cmd, " > ", path_temp, "/net.tree")
      system(command = cmd)
    }else if(os == "Windows"){
      cmd=paste0(binpath, "/bin/LOUVAIN/louvain_win.exe ", cmd)
      tree=system(command = cmd, intern=TRUE)
      cat(tree[1:(length(tree)-1)], file = paste0(path_temp, "/net.tree"), sep = "\n")
    }else if(os == "Darwin"){
      cmd=paste0(binpath, "/bin/LOUVAIN/louvain_mac ", cmd, " > ", path_temp, "/net.tree")
      system(command = cmd)
    }else{
      stop("Linux, Windows or Mac distributions only.")
    }

    # Control: if the command line did not work
    if(!("net.tree" %in% list.files(paste0(path_temp)))){
      stop("Command line was wrongly implemented. Louvain did not run.")
    }

    # Retrieve output from net.tree
    tree=utils::read.table(paste0(path_temp, "/net.tree"))

    id0=which(tree[,1]==0)
    tree=tree[(id0[1]+1):(id0[2]-1),]

    com=data.frame(ID=idnode[,2], Com=0)
    com[match(tree[,1],idnode[,1]),2]=tree[,2]

    comc=com

    # Remove temporary file
    if(delete_temp){
      unlink(paste0(path_temp), recursive = TRUE)
    }

  }

  if(lang=="all"){
    com=cbind(comr,comc[,-1])
  }

  # Rename and reorder columns
  com[, 1] <- as.character(com[, 1])
  com <- knbclu(com)

  # Return output
  return(com)

}