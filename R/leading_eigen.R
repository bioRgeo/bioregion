#' Finding communities based on leading eigenvector of the community matrix
#'
#' This function finds communities in a (un)weighted undirected network based on leading eigenvector
#' of the community matrix.
#'
#' @param net a two- or three-column \code{data.frame} representing a network with the two first columns
#' as (un)directed links between pair of nodes and an optional third column indicating the weight of the link
#' @param weight a boolean indicating if the weights should be considered if there is a third column
#' @export
#' @details
#' This function is based on leading eigenvector of the community matrix \insertCite{Newman2006}{bioRgeo}
#' as implemented in the \href{https://cran.r-project.org/web/packages/igraph/index.html}{igraph} package
#' (\link[igraph]{cluster_leading_eigen}).
#'
#' @return A \code{data.frame} providing one community by node.
#'
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}),
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) and
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' @seealso \link{infomap}, \link{oslom}
#' @examples
#' comat=matrix(sample(1000,50),5,10)
#' rownames(comat)=paste0("Site",1:5)
#' colnames(comat)=paste0("Species",1:10)
#'
#' net=similarity(comat,metric="Simpson")
#' com=leading_eigen(net)
#' @references
#' \insertRef{Newman2006}{bioRgeo}
#' @export
leading_eigen <- function(net, weight = TRUE){

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
      stop("The third column of net must be numeric")
    }
  }

  if(!is.logical(weight)){
    stop("weight must be a boolean")
  }

  # Prepare input
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

  # Run algo
  net=igraph::graph_from_data_frame(netemp, directed=FALSE)
  comtemp=igraph::cluster_leading_eigen(net)
  comtemp=cbind(as.numeric(comtemp$names),as.numeric(comtemp$membership))

  com=data.frame(ID=idnode[,2], Com=0)
  com[match(comtemp[,1],idnode[,1]),2]=comtemp[,2]

  # Return output
  com[,1]=as.character(com[,1])
  return(com)

}
