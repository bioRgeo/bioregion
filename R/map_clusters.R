#' Create a map of bioregions
#'
#' This plot function can be used to visualise bioregions based on a 
#' bioRgeo.clusters object combined with a geometry (sf objects). 
#'
#' @param clusters an object of class \code{bioRgeo.clusters} or a 
#' \code{data.frame}. If a \code{data.frame} is used, the first column should
#'  represent the sites' ID, and the next column(s) the clusters.
#' @param geometry a spatial object that can be handled by the \code{sf}
#' package. The first attribute should correspond to the sites' ID (see Details).
#' @param write_clusters a \code{boolean} indicating if the \code{clusters} 
#' should be added in \code{geometry}.
#' @param plot a \code{boolean} indicating if the plot should be drawn.
#' @export
#' @details
#' The \code{clusters} and \code{geometry} site IDs should correspond. They should
#' have the same type (i.e. \code{character} is cluster is a 
#' \code{bioRgeo.clusters} object) and the site of \code{clusters} should be 
#' included in the sites of \code{geometry}. 
#' @return  One or several maps of bioregions if \code{plot = TRUE} and the 
#' geometry with additional clusters' attributes if \code{write_clusters = TRUE}. 
#' 
#' @author
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}),
#' Boris Leroy (\email{leroy.boris@gmail.com}) and
#' Pierre Denelle (\email{pierre.denelle@gmail.com})
#'
#' @examples
#' data(vegemat)
#' data(vegesp)
#' 
#' net=similarity(vegemat,metric=c("Simpson"))
#' clu=netclu_greedy(net)
#' map <- map_clusters(clu, vegesp, write_clusters = TRUE, plot = FALSE)
#' 
#' @export
map_clusters <- function(clusters, geometry, write_clusters = FALSE, plot = TRUE) {

  # Control clusters 
  if (inherits(clusters, "bioRgeo.clusters")) {
    clu=TRUE
    df=clusters$clusters
  }else{  
    # data.frame
    if (!is.data.frame(clusters)) {
      stop("If not a bioRgeo.clusters's object, clusters must be a data.frame.",
           call. = FALSE)
    }
    # at least two columns
    if (dim(clusters)[2] < 2) {
      stop("clusters must be a data.frame with at least two columns.", call. = FALSE)
    }
    # no duplictaed ID
    if (sum(duplicated(clusters[,1])) > 0) {
      message("Duplicated site ID detected!")
    }
    # no NAs
    nbna <- sum(is.na(clusters))
    if (nbna > 0) {
      stop("NA(s) detected in the data.frame!", call. = FALSE)
    }
    clu=FALSE
    df=clusters
  }
  
  # Control geometry
  if(class(geometry)[1]!="sf"){
    stop("It seems that the geometry used is not an sf object.",
         call. = FALSE)
  }
  
  # Control clusters in geometry
  idc <- df[,1]
  idg <- geometry[,1,drop=TRUE]
  if(length(setdiff(idc,idg))>0){
    stop("The site of clusters should be included in the sites of geometry",
         call. = FALSE)
  }
  
  # Control parameters
  controls(args=write_clusters, data=NULL, type="boolean")
  controls(args=plot, data=NULL, type="boolean")
  
  # Prepare geometry
  sp=geometry[match(idc,idg),]
  nbsp=dim(sp)[2]
  nbdf=dim(df)[2]
  sp=cbind(sp,df[,-1])
  colnames(sp)[nbsp:(nbsp+nbdf-2)]=colnames(df)[-1]
  
  # Plot
  if(plot){
    
    geomsp=sf::st_geometry(sp)
    plotsp=sp[,-(1:(nbsp-1))]
    nbplotsp=dim(plotsp)[2]-1
    
    if(nbplotsp==1){ 
      plot(plotsp)
    }else{ 
      mod4q=floor(nbplotsp/4)
      mod4r=nbplotsp-mod4q*4

      if(mod4q==0){
        plot(plotsp)
      }else{
        for(k in 1:mod4q){
          grDevices::dev.new()
          plot(plotsp[((k-1)*4+1):((k-1)*4+4)])
        }
        if(mod4r>0){
          grDevices::dev.new()
          plot(plotsp[((mod4q)*4+1):((mod4q)*4+mod4r)], key.pos=NULL)
        }
        
      }
    }
    
  }
  
  # Write
  if(write_clusters){
    return(sp)
  }

}
