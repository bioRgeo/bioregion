#' Infomap community finding
#'
#' This function finds communities in a (un)weighted (un)directed network based on the Infomap algorithm
#' (\url{https://github.com/mapequation/infomap}, version 1.6.0).
#'
#' @param net a two- or three-column \code{data.frame} representing a network with the two first columns
#' as (un)directed links between pair of nodes and an optional third column indicating the weight of the link
#' @param weight a boolean indicating if the weights should be considered if there is a third column
#' @param bipartite a boolean indicating if the network is bipartite (see Details)
#' @param nbmod penalize solutions the more they differ from this number (0 by default for no preferred number of modules)
#' @param markovtime scales link flow to change the cost of moving between modules, higher values results
#' in fewer modules (default is 1)
#' @param seed for the random number generator
#' @param twolevel a boolean indicating if the algorithm should optimize a two-level partition of the network
#' (default is multi-level)
#' @param directed a boolean indicating if the network is directed (from column 1 to column 2)
#' @param delete_temp a boolean indicating if the temporary folder should be removed (see Details)
#' @param path_temp a string indicating the path to the temporary folder (see Details)
#' @param binpath a string indicating the path to the bin folder (see \link{bin} and Details)
#' @export
#' @details
#' Infomap is a network clustering algorithm based on the Map equation proposed in
#' \insertCite{Rosvall2008}{bioRgeo} that finds communities in (un)weighted and (un)directed networks.
#' Infomap has two ways to deal with bipartite networks. The first possibility is to consider the bipartite network
#' as unipartite network. The second possibility is to set the \code{bipartite} argument to TRUE in order to
#' approximate a two-step random walker (see \url{https://www.mapequation.org/infomap/} for more information).
#'
#' This function is based on the 1.6.0 C++ version of Infomap (\url{https://github.com/mapequation/infomap/releases}).
#' This function needs executables files to run. They can be installed with \link{bin}. If you set the the path to
#' the folder that will host the bin folder  manually while running \link{bin} please make sure to set \code{binpath}
#' accordingly.
#'
#' The C++ version of Infomap generates temporary folders and/or files that are stored in the \code{path_temp} folder
#' (folder "infomap_temp" in the workind directory by default). This temporary folder is removed by default
#' (\code{delete_temp = TRUE}).
#'
#' @return A \code{data.frame} providing one partition by hierarchical level.
#'
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}),
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) and
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' @seealso \link{louvain}, \link{oslom}
#' @examples
#' comat=matrix(sample(1000,50),5,10)
#' rownames(comat)=paste0("Site",1:5)
#' colnames(comat)=paste0("Species",1:10)
#'
#' net=spproject(comat,metric="Simpson")
#' #com=infomap(net) # run bin() to use this function
#' @references
#' \insertRef{Rosvall2008}{bioRgeo}
#' @export
infomap <- function(net, weight = TRUE, bipartite= FALSE, nbmod = 0, markovtime = 1, seed = 1, twolevel = FALSE,
                    directed = FALSE, delete_temp = TRUE, path_temp = "infomap_temp", binpath = NULL){

  # Remove warning for tidyr
  defaultW <- getOption("warn")
  options(warn = -1)

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

  # Check if INFOMAP is present
  if(!file.exists(paste0(binpath,"/bin/INFOMAP/infomap_lin"))){
    stop("Infomap is not in bioRgeo/bin directory...")
  }
  if(!file.exists(paste0(binpath,"/bin/INFOMAP/infomap_win.exe"))){
    stop("Infomap is not in bioRgeo/bin directory...")
  }
  #if(!file.exists(paste0(binpath,"/bin/INFOMAP/infomap_mac"))){
  #  stop("Infomap is not in bioRgeo/bin directory...")
  #}

  # Controls
  if(!is.data.frame(net)){
    stop("net must be a two- or three-columns data.frame")
  }

  if(dim(net)[2]!=2 & dim(net)[2]!=3){
    stop("net must be a two- or three-columns data.frame")
  }

  sco=sum(is.na(net))
  if(sco>0){
    stop("NA(s) detected in the data.frame")
  }

  if(weight & dim(net)[2]==2){
    stop("net must be a three-columns data.frame if weight equal TRUE")
  }

  if(weight & dim(net)[2]==3){
    if(class(net[,3])!="numeric" & class(net[,3])!="integer"){
      stop("The third column of df must be numeric")
    }
  }

  if(!is.logical(weight)){
    stop("weight must be a boolean")
  }

  if(!is.logical(bipartite)){
    stop("bipartite must be a boolean")
  }

  if(!is.numeric(nbmod)){
    stop("nbmod must be numeric")
  }else{
    if(nbmod<0){
      stop("nbmod must be positive")
    }
    if((nbmod-floor(nbmod))>0){
      stop("nbmod must be an integer higher than 0")
    }
  }

  if(!is.numeric(seed)){
    stop("seed must be numeric")
  }else{
    if(seed<=0){
      stop("seed must be strictly higher than 0")
    }
    if((seed-floor(seed))>0){
      stop("seed must be an integer higher or equal to 0")
    }
  }

  if(!is.numeric(markovtime)){
    stop("markovtime must be numeric")
  }else{
    if(markovtime<=0){
      stop("markovtime must be strictly higher than 0")
    }
  }

  if(!is.logical(twolevel)){
    stop("twolevel must be a boolean")
  }

  if(!is.logical(directed)){
    stop("directed must be a boolean")
  }

  if(!is.logical(delete_temp)){
    stop("delete_temp must be a boolean")
  }

  if(!is.character(path_temp)){
    stop("path_temp must be a string")
  }

  # Create temp folder
  dir.create(path_temp, showWarnings = FALSE, recursive = TRUE)
  if(!file.exists(path_temp)){
    stop(paste0("Impossible to create directory ", path_temp))
  }

  # Prepare input for INFOMAP
  idnode1=as.character(net[,1])
  idnode2=as.character(net[,2])
  if(bipartite){
    # Control bipartite
    if(length(intersect(idnode1,idnode2))>0){
      stop("If bipartite = TRUE primary and feature node should be different.")
    }
    idnode1=idnode1[!duplicated(idnode1)]
    idnode1=data.frame(ID=1:length(idnode1),ID_NODE=idnode1,Type=1) # Primary nodes
    idnode2=idnode2[!duplicated(idnode2)]
    idnode2=data.frame(ID=((dim(idnode1)[1]+1):(dim(idnode1)[1]+length(idnode2))),ID_NODE=idnode2,Type=2) # Feature nodes
    N=dim(idnode1)[1]+1 # First node id of the feature node type
    idnode=rbind(idnode1,idnode2)
  }else{
    idnode=c(idnode1,idnode2)
    idnode=idnode[!duplicated(idnode)]
    idnode=data.frame(ID=1:length(idnode),ID_NODE=idnode)
  }

  netemp=data.frame(node1=idnode[match(net[,1],idnode[,2]),1],node2=idnode[match(net[,2],idnode[,2]),1])
  if(weight){
    netemp=cbind(netemp,net[,3])
    netemp=netemp[netemp[,3]>0,]
  }

  # Export input in INFOMAP folder
  if(bipartite){ # Add tag if bipartite
    cat(paste0("*Bipartite ",N),"\n",file=paste0(path_temp,"/net.txt"))
    utils::write.table(netemp, paste0(path_temp,"/net.txt"), append=TRUE, row.names=FALSE, col.names=FALSE, sep=" ")
  }else{
    utils::write.table(netemp, paste0(path_temp,"/net.txt"), row.names=FALSE, col.names=FALSE, sep=" ")
  }

  # Prepare command to run INFOMAP
  cmd=paste0("--silent --seed ", seed," --preferred-number-of-modules ", nbmod, " --markov-time ", markovtime)
  if(twolevel){
    cmd=paste0(cmd, " --two-level")
  }
  if(directed){
    cmd=paste0(cmd, " --flow-model directed")
  }else{
    cmd=paste0(cmd, " --flow-model undirected")
  }
  if(bipartite){
    cmd=paste0("-i bipartite ", cmd, " ", path_temp, "/net.txt ", path_temp)
  }else{
    cmd=paste0("-i link-list ", cmd, " ", path_temp, "/net.txt ", path_temp)
  }

  if(os == "Linux"){
    cmd=paste0(binpath, "/bin/INFOMAP/infomap_lin ", cmd)
  }else if(os == "Windows"){
    cmd=paste0(binpath, "/bin/INFOMAP/infomap_win.exe ", cmd)
  }else if(os == "Darwin"){
    stop("TO IMPLEMENT")
  }else{
    stop("Linux, Windows or Mac distributions only.")
  }

  # Run INFOMAP
  system(command = cmd)

  # Control: if the command line did not work
  if(!("net.tree" %in% list.files(paste0(path_temp)))){
    stop("Command line was wrongly implemented. Infomap did not run.")
  }

  # Retrieve output from net.tree
  tree=utils::read.table(paste0(path_temp,"/net.tree"))

  # Reformat tree [TO COMMENT]
  idinf=as.numeric(tree[,4])  # INFOMAP node ids

  # require(tidyr) # Extract the modules from tree
  df <- data.frame(x = as.character(tree[,1]))
  cominf = tidyr::separate(df, "x", as.character(1:100), sep = ":") # Max 100 levels
  cominf[is.na(cominf)]=0
  for(k in 1:dim(cominf)[2]){ # Transform in numeric
    cominf[,k]=as.numeric(as.character(cominf[,k]))
  }
  cominf=cominf[,apply(cominf,2,sum)>0] # Data frame with information contains in the first column of tree (0 when no info)
  nblev=dim(cominf)[2] # The number of columns of cominf correspond to the number of levels + one dummy column

  for(k in 2:nblev){ # Set a value 0 for the dummy column
    cominf[cominf[,k]==0,(k-1)]=0
  }
  cominf[cominf[,nblev]>0,nblev]=0

  for(k in 2:nblev){ # Extract an real partition for lower levels
    cominf[,k]=as.numeric(factor(paste0(cominf[,k-1],"_",cominf[,k])))
  }

  cominf=cominf[,nblev:1] # Reverse column order

  com=data.frame(ID=idnode[,2], dum=0) # Dummy level
  com[match(idinf,idnode[,1]),2]=cominf[,1]

  if(nblev>=2){
    com$Com=0
    com[match(idinf,idnode[,1]),3]=cominf[,2]
  }
  if(nblev>=3){
    com$HCom=0
    com[match(idinf,idnode[,1]),4]=cominf[,3]
  }
  if(nblev>=4){
    com$HHCom=0
    com[match(idinf,idnode[,1]),5]=cominf[,4]
  }
  if(nblev>=5){
    com$HHHCom=0
    com[match(idinf,idnode[,1]),6]=cominf[,5]
  }
  if(nblev>=6){
    com$HHHHCom=0
    com[match(idinf,idnode[,1]),7]=cominf[,6]
  }

  com=com[,-2] # Remove dum

  # Remove temporary
  if(delete_temp){
    unlink(paste0(path_temp), recursive = TRUE)
  }

  # Rename and reorder columns
  com=com[,c(1,dim(com)[2]:2)]
  colnames(com)[2:dim(com)[2]]=paste0("Com",1:(dim(com)[2]-1))

  # Put the warning back
  options(warn = defaultW)

  # Return output
  com[,1]=as.character(com[,1])
  return(com)

}
