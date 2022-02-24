#' Create a contingency table from a data.frame
#'
#' This function creates a contingency table from a two- or three-columns \code{data.frame} where
#' each row represents the interaction between two objects (site and species for example)
#' and an optional third column indicating the weight of the interaction (if \code{weight = TRUE}).
#'
#' @param df a two- or three-columns \code{data.frame} where
#' each row represents the interaction between two objects (site and species for example)
#' and an optional third column indicating the weight of the interaction
#' @param weight a boolean indicating if the weight should be considered
#' @param squared a boolean indicating if the output matrix should but squared (same objects in rows and columns)
#' @param symmetrical a boolean indicating if the resulting matrix should be symmetrical (only if \code{squared = TRUE}).
#' Note that different weights associated with two opposite pairs already present in df will be preserved.
#' @param value the value to assign to the pairs of objects not present in df (0 by default)
#' @export
#' @return A \code{matrix} with the first objects (first column of \code{df}) as rows and
#' the second objects (second column of \code{df}) as columns. Note that if \code{squared = TRUE} the rows and columns
#' have the same number of elements corresponding to the concatenation of unique objects in  \code{df}'s first and second
#' columns. If \code{squared = TRUE} the matrix can be forced to be symetrical based on the upper triangular part of
#' the matrix.
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}),
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) and
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' @seealso \link{mat_to_net}
#' @examples
#' df <- data.frame(Site = c(rep("A", 2), rep("B", 3), rep("C", 2)),
#' Species = c("a", "b", "a", "c", "d", "b", "d"),
#' Weight = c(10, 100, 1, 20, 50, 10, 20))
#'
#' comat=net_to_mat(df,weight=TRUE)
#' @export
net_to_mat <- function(df, weight = FALSE, squared = FALSE, symmetrical = FALSE, value = 0){

  # Controls
  if(!is.data.frame(df)){
    stop("df must be a two- or three-columns data.frame")
  }

  sco=sum(is.na(df))
  if(sco>0){
    stop("NA(s) detected in the data.frame")
  }

  if(!is.logical(weight)){
    stop("weight must be a boolean")
  }

  if(!is.logical(squared)){
    stop("squared must be a boolean")
  }

  if(!is.logical(symmetrical)){
    stop("symmetrical must be a boolean")
  }

  if(!squared & symmetrical){
    stop("symmetrical only for squared matrix!")
  }

  if(dim(df)[2]!=2 & dim(df)[2]!=3){
    stop("df must be a two- or three-columns data.frame")
  }

  if(weight & dim(df)[2]==2){
    stop("df must be a three-columns data.frame if weight equal TRUE")
  }

  if(weight & dim(df)[2]==3){
    if(class(df[,3])!="numeric" & class(df[,3])!="integer"){
      stop("The third column of df must be numeric")
    }
  }

  # Rename columns
  colnames(df)[1:2]=c("Object1","Object2")

  # Transform objects in character
  df$Object1=as.character(df$Object1)
  df$Object2=as.character(df$Object2)

  # Extract unique objects in both columns
  idobj1=df$Object1[!duplicated(df$Object1)]
  idobj2=df$Object2[!duplicated(df$Object2)]

  # Manage weight
  if(weight){
    colnames(df)[3] <- "Weight"
  } else{
    df=df[,1:2]
    df$Weight <- 1
  }

  # Modify df to obtain a squared contingency table
  if(squared){
    diff12=setdiff(idobj1,idobj2) # Objects contains in the first column but not in the second
    diff21=setdiff(idobj2,idobj1) # Objects contains in the second column but not in the first

    # Create a new dataframe containing all objects in both columns
    obj12=data.frame(Object1=df$Object1[1], Object2=diff12, Weight=NA)
    obj21=data.frame(Object1=diff21, Object2=df$Object2[1], Weight=NA)
    df=rbind(df,obj12,obj21)

    # Extract unique objects in both columns (should be the same at this stage)
    idobj1=df$Object1[!duplicated(df$Object1)]
    idobj2=idobj1
  }

  # Create contingency table from df
  df$Object1 <- factor(df$Object1,levels=idobj1)
  df$Object2 <- factor(df$Object2,levels=idobj2)

  comat <- with(df, {
    out <- matrix(nrow = nlevels(Object1), ncol = nlevels(Object2),
                  dimnames = list(levels(Object1), levels(Object2)))
    out[cbind(Object1, Object2)] <- Weight
    out
  })

  # Force the matrix to be symmetrical if squared = TRUE
  if(squared & symmetrical){
    low=base::lower.tri(comat)
    up=base::upper.tri(comat)
    nacomat=is.na(comat)
    temp=t(comat)
    comat[low & nacomat]=temp[low & nacomat]
    comat[up & nacomat]=temp[up & nacomat]
  }

  # Replace NAs with 0s
  comat[is.na(comat)] <- value

  # Check for empty rows and columns if squared = FALSE
  #if(!squared){
  #  comat <- comat[rowSums(comat) > 0,]
  #  comat <- comat[, colSums(comat) > 0]
  #}

  # Return the contingency matrix
  return(comat)
}
