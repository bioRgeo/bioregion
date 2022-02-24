#' Create a data.frame from a contingency table
#'
#' This function creates a two- or three-columns \code{data.frame} where
#' each row represents the interaction between two objects (site and species for example)
#' and an optional third column indicating the weight of the interaction (if \code{weight = TRUE})
#' from a contingency table (sites as rows and species as columns for example).
#'
#' @param comat a contingency table (i.e. \code{matrix})
#' @param weight a boolean indicating if the value are weights
#' @param remove_absent_objects a boolean determining whether absent
#' objects from the contingency table have to be removed from the output
#' @export
#' @return A \code{data.frame} where each row represents the interaction between
#' two objects and an optional third column indicating the weight of the interaction.
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}),
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) and
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' @seealso \link{net_to_mat}
#'
#' @examples
#' comat=matrix(sample(1000,50),5,10)
#' rownames(comat)=paste0("Site",1:5)
#' colnames(comat)=paste0("Species",1:10)
#'
#' df=mat_to_net(comat,weight=TRUE)
#' @export
mat_to_net <- function(comat, weight = FALSE, remove_absent_objects = TRUE){

  # Controls
  if(!is.matrix(comat)){
    stop("Contingency table should be a matrix")
  }

  sco=sum(is.na(comat))
  if(sco>0){
    stop("NA(s) detected in the contingency table")
  }

  if(!is.logical(weight)){
    stop("weight must be a boolean")
  }

  if(!is.logical(remove_absent_objects)){
    stop("remove_absent_objects must be a boolean")
  }

  # Conversion as data.frame
  df <- reshape2::melt(comat)
  colnames(df) <- c("Object1", "Object2", "Weight")

  # Remove interactions with weight equal 0
  if(remove_absent_objects == TRUE){
    df <- df[df$Weight != 0,]
  }

  # Remove the weight column if weight is set to FALSE
  if(!weight){
    df=df[,-3]
  }

  # Reorder by Object 1 and 2
  df$Object1=factor(df$Object1,levels=rownames(comat))
  df$Object2=factor(df$Object2,levels=colnames(comat))
  df=df[order(df$Object1,df$Object2),]

  # Transform the two first columns in character
  df[,1]=as.character(df[,1])
  df[,2]=as.character(df[,2])

  # Return the data.frame
  return(df)
}
