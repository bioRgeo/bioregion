#' Create a map of bioregions
#'
#' This plot function can be used to visualize bioregions based on a 
#' `bioregion.clusters` object combined with a geometry (`sf` objects). 
#'
#' @param clusters An object of class `bioregion.clusters` or a `data.frame`. 
#' If a `data.frame` is used, the first column should represent the sites' 
#' ID, and the subsequent column(s) should represent the clusters.
#'  
#' @param geometry A spatial object that can be handled by the `sf` package. 
#' The first attribute should correspond to the sites' ID (see Details).
#' 
#' @param write_clusters A `boolean` indicating if the `clusters` should be 
#' added to the `geometry`.
#' 
#' @param plot A `boolean` indicating if the plot should be drawn.
#' 
#' @param ... Further arguments to be passed to `sf::plot()`.
#' 
#' @return One or several maps of bioregions if `plot = TRUE` and the 
#' geometry with additional clusters' attributes if `write_clusters = TRUE`.
#' 
#' @details
#' The `clusters` and `geometry` site IDs should correspond. They should 
#' have the same type (i.e., `character` if `clusters` is a 
#' `bioregion.clusters` object) and the sites of `clusters` should be 
#' included in the sites of `geometry`.
#' 
#' @author
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) \cr
#' Boris Leroy (\email{leroy.boris@gmail.com}) \cr
#' Pierre Denelle (\email{pierre.denelle@gmail.com})
#'
#' @examples
#' data(fishmat)
#' data(fishsf)
#' 
#' net <- similarity(fishmat, metric = "Simpson")
#' clu <- netclu_greedy(net)
#' map <- map_bioregions(clu, fishsf, write_clusters = TRUE, plot = FALSE)
#' 
#' @importFrom sf st_geometry
#' 
#' @export

map_bioregions <- function(clusters, 
                           geometry, 
                           write_clusters = FALSE,
                           plot = TRUE, 
                           ...) {

  controls(args = write_clusters, data = NULL, type = "boolean")
  controls(args = plot, data = NULL, type = "boolean")
  
  # Control clusters 
  if (inherits(clusters, "bioregion.clusters")) {
    clu <- TRUE
    df <- clusters$clusters
  }else{  
    # data.frame
    if (!is.data.frame(clusters)) {
      stop(
        "If not a bioregion.clusters's object, clusters must be a data.frame.",
        call. = FALSE)
    }
    # at least two columns
    if (dim(clusters)[2] < 2) {
      stop("clusters must be a data.frame with at least two columns.",
           call. = FALSE)
    }
    # no duplictaed ID
    if (sum(duplicated(clusters[,1])) > 0) {
      message("Duplicated site ID detected!")
    }
    # no NAs
    nbna <- sum(is.na(clusters))
    if (nbna > 0) {
      stop("NA(s) detected in the data.frame!", 
           call. = FALSE)
    }
    clu <- FALSE
    df <- clusters
  }
  
  # Control geometry
  if(class(geometry)[1] != "sf"){
    stop("It seems that the geometry used is not an sf object.",
         call. = FALSE)
  }
  
  # Control clusters in geometry
  idc <- df[,1]
  idg <- geometry[, 1, drop = TRUE]
  if(length(setdiff(idc, idg)) > 0){
    stop("The site of clusters should be included in the sites of geometry.",
         call. = FALSE)
  }
  
  # Control parameters
  controls(args = write_clusters, data = NULL, type = "boolean")
  controls(args = plot, data = NULL, type = "boolean")
  
  # Prepare geometry
  sp <- geometry[match(idc, idg), ]
  nbsp <- dim(sp)[2]
  nbdf <- dim(df)[2]
  sp <- cbind(sp, df[, -1])
  colnames(sp)[nbsp:(nbsp+nbdf-2)] <- colnames(df)[-1]
  
  # Plot
  if(plot){
    
    geomsp <- sf::st_geometry(sp)
    plotsp <- sp[, -(1:(nbsp-1))]
    nbplotsp <- dim(plotsp)[2]-1
    
    if(nbplotsp == 1){ 
      plot(plotsp)
    }else{ 
      mod4q <- floor(nbplotsp/4)
      mod4r <- nbplotsp-mod4q*4

      if(mod4q == 0){
        plot(plotsp, ...)
      }else{
        for(k in 1:mod4q){
          grDevices::dev.new()
          plot(plotsp[((k-1)*4+1):((k-1)*4+4)], ...)
        }
        if(mod4r>0){
          grDevices::dev.new()
          plot(plotsp[((mod4q)*4+1):((mod4q)*4+mod4r)], key.pos=NULL, ...)
        }
        
      }
    }
  }
  
  # Write
  if(write_clusters){
    return(sp)
  }
}
