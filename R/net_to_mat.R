#' Create a contingency table from a data.frame
#'
#' This function creates a contingency table from a two- or three-columns 
#' \code{data.frame} where each row represents the interaction between two
#' nodes (site and species for example) and an optional third column indicating
#' the weight of the interaction (if \code{weight = TRUE}).
#'
#' @param net a two- or three-columns \code{data.frame} where each row
#' represents the interaction between two nodes (site and species for example)
#' and an optional third column indicating the weight of the interaction.
#' 
#' @param weight a \code{boolean} indicating if the weight should be considered
#' 
#' @param squared a \code{boolean} indicating if the output matrix should but 
#' squared (same nodes in rows and columns).
#' 
#' @param symmetrical a \code{boolean} indicating if the resulting matrix
#' should be symmetrical (only if \code{squared = TRUE}).
#' Note that different weights associated with two opposite pairs already
#' present in net will be preserved.
#' 
#' @param missing_value the value to assign to the pairs of nodes not present
#' in net (0 by default).
#' 
#' @return A \code{matrix} with the first nodes (first column of \code{net}) as 
#' rows and the second nodes (second column of \code{net}) as columns. Note
#' that if \code{squared = TRUE} the rows and columns have the same number of
#' elements corresponding to the concatenation of unique objects in
#' \code{net}'s first and second columns. If \code{squared = TRUE} the matrix
#' can be forced to be symmetrical based on the upper triangular part of the
#' matrix.
#' 
#' @author
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}),
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) and
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' 
#' @seealso \link{mat_to_net}
#' 
#' @examples
#' \dontrun{
#' net <- data.frame(
#'   Site = c(rep("A", 2), rep("B", 3), rep("C", 2)),
#'   Species = c("a", "b", "a", "c", "d", "b", "d"),
#'   Weight = c(10, 100, 1, 20, 50, 10, 20)
#' )
#'
#' mat <- net_to_mat(net, weight = TRUE)
#' }
#' 
#' @export
net_to_mat <- function(net, weight = FALSE, squared = FALSE,
                       symmetrical = FALSE, missing_value = 0){
  
  # Control input net
  controls(args = NULL, data = net, type = "input_net")
  
  # Control parameters
  controls(args = weight, data = net, type = "input_net_weight")
  controls(args = squared, data = NULL, type = "boolean")
  controls(args = symmetrical, data = NULL, type = "boolean")
  controls(args = missing_value, data = NULL, type = "numeric")
  if (!squared & symmetrical) {
    stop("symmetrical only for squared matrix!", call. = FALSE)
  }
  if (weight & dim(net)[2] == 3) {
    if (!is.numeric(net[, 3])) {
      stop("The third column of net must be numeric", call. = FALSE)
    }
  }
  
  # Rename columns
  colnames(net)[1:2] <- c("Object1", "Object2")
  
  # Transform objects in character
  net$Object1 <- as.character(net$Object1)
  net$Object2 <- as.character(net$Object2)
  
  # Extract unique objects in both columns
  idobj1 <- net$Object1[!duplicated(net$Object1)]
  idobj2 <- net$Object2[!duplicated(net$Object2)]
  
  # Manage weight
  if (weight) {
    colnames(net)[3] <- "Weight"
  } else {
    net <- net[, 1:2]
    net$Weight <- 1
  }
  
  # Modify net to obtain a squared contingency table
  if (squared) {
    diff12 <- setdiff(idobj1, idobj2) # Objects contains in the first column
    # but not in the second
    diff21 <- setdiff(idobj2, idobj1) # Objects contains in the second column
    # but not in the first
    
    # Create a new dataframe containing all objects in both columns
    obj12 <- data.frame(Object1 = net$Object1[1], Object2 = diff12,
                        Weight = NA)
    obj21 <- data.frame(Object1 = diff21, Object2 = net$Object2[1],
                        Weight = NA)
    net <- rbind(net, obj12, obj21)
    
    # Extract unique objects in both columns (should be the same at this stage)
    idobj1 <- net$Object1[!duplicated(net$Object1)]
    idobj2 <- idobj1
  }
  
  # Create contingency table from net
  net$Object1 <- factor(net$Object1, levels = idobj1)
  net$Object2 <- factor(net$Object2, levels = idobj2)
  
  mat <- with(net, {
    out <- matrix(
      nrow = nlevels(Object1), ncol = nlevels(Object2),
      dimnames = list(levels(Object1), levels(Object2))
    )
    out[cbind(Object1, Object2)] <- Weight
    out
  })
  
  # Force the matrix to be symmetrical if squared = TRUE
  if (squared & symmetrical) {
    low <- base::lower.tri(mat)
    up <- base::upper.tri(mat)
    namat <- is.na(mat)
    temp <- t(mat)
    mat[low & namat] <- temp[low & namat]
    mat[up & namat] <- temp[up & namat]
  }
  
  # Replace NAs with 0s
  mat[is.na(mat)] <- missing_value
  
  # Check for empty rows and columns if squared = FALSE
  # if(!squared){
  #  mat <- mat[rowSums(mat) > 0,]
  #  mat <- mat[, colSums(mat) > 0]
  # }
  
  # Return the contingency matrix
  return(mat)
}
