#' Community structure detection in weighted bipartite network via modularity optimisation
#'
#' This function takes a bipartite weighted graph and computes modules by applying Newman’s modularity
#' measure in a bipartite weighted version to it.
#'
#' @param net a two- or three-column \code{data.frame} representing a network with the two first columns
#' as undirected links between pair of nodes and an optional third column indicating the weight of the link
#' @param weight a boolean indicating if the weights should be considered if there is a third column
#' @param forceLPA a boolean indicating if the even faster pure LPA-algorithm of Beckett should be used? DIRT-
#' LPA, the default, is less likely to get trapped in a local minimum, but is slightly
#' slower. Defaults to FALSE.
#' @export
#' @details
#' This function is based on the modularity optimization algorithm provided by Stephen Beckett \insertCite{Beckett2016}{bioRgeo}
#' as implemented in the \href{https://cran.r-project.org/web/packages/bipartite/index.html}{bipartite} package
#' (\link[bipartite]{computeModules}).
#'
#' @return A \code{data.frame} providing one community by node.
#'
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}),
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) and
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' @seealso \link{infomap}, \link{oslom}
#' @examples
#' net <- data.frame(Site = c(rep("A", 2), rep("B", 3), rep("C", 2)),
#' Species = c("a", "b", "a", "c", "d", "b", "d"),
#' Weight = c(10, 100, 1, 20, 50, 10, 20))
#'
#' com=beckett(net)
#' @references
#' \insertRef{Beckett2016}{bioRgeo}
#' @export
beckett <- function(net, weight = TRUE, forceLPA = FALSE){

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

  if(!is.logical(forceLPA)){
    stop("forceLPA must be a boolean")
  }

  # Prepare input
  idnode1=as.character(net[,1])
  idnode2=as.character(net[,2])

  if(length(intersect(idnode1,idnode2))>0){ # Control bipartite
    stop("The network should be bipartite!")
  }

  idnode=c(idnode1,idnode2)
  idnode=idnode[!duplicated(idnode)]
  idnode=data.frame(ID=1:length(idnode),ID_NODE=idnode)

  netemp=data.frame(node1=idnode[match(net[,1],idnode[,2]),1],node2=idnode[match(net[,2],idnode[,2]),1])
  if(weight){
    netemp=cbind(netemp,net[,3])
    netemp=netemp[netemp[,3]>0,]
    colnames(netemp)[3]="weight"
  }

  # Transform netemp into a contingency table
  comat=df_to_contingency(netemp, weight=weight)

  # Run algo
    comtemp=bipartite::computeModules(comat, forceLPA = forceLPA)@modules[-1,-c(1,2)]
  comtemp=cbind(c(as.numeric(rownames(comat)),as.numeric(colnames(comat))),apply(comtemp,2,function(x) which(x>0)))

  com=data.frame(ID=idnode[,2], Com=0)
  com[match(comtemp[,1],idnode[,1]),2]=comtemp[,2]

  # Return output
  com[,1]=as.character(com[,1])
  return(com)

}
