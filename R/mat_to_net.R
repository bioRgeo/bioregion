#' Create a data.frame from a contingency table
#'
#' This function creates a two- or three-columns \code{data.frame} where
#' each row represents the interaction between two nodes (site and species for 
#' example) and an optional third column indicating the weight of the interaction 
#' (if \code{weight = TRUE}) from a contingency table (sites as rows and species 
#' as columns for example).
#'
#' @param mat a contingency table (i.e. \code{matrix}).
#' @param weight a \code{boolean} indicating if the value are weights.
#' @param remove_zeroes a \code{boolean} determining whether interactions with 
#' weight equal to 0 should be removed from the output.
#' @export
#' @return A \code{data.frame} where each row represents the interaction between
#' two nodes and an optional third column indicating the weight of the interaction.
#' @author
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}),
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) and
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' @seealso \link{net_to_mat}
#'
#' @examples
#' mat <- matrix(sample(1000, 50), 5, 10)
#' rownames(mat) <- paste0("Site", 1:5)
#' colnames(mat) <- paste0("Species", 1:10)
#'
#' net <- mat_to_net(mat, weight = TRUE)
#' @export
mat_to_net <- function(mat, weight = FALSE, remove_zeroes = TRUE) {

  # Control input mat
  controls(args=NULL, data=mat, type="input_matrix")
  
  # Control parameters
  controls(args=weight, data=NULL, type="boolean")
  controls(args=remove_zeroes, data=NULL, type="boolean")

  # Conversion as data.frame
  net <- reshape2::melt(mat)
  colnames(net) <- c("Node1", "Node2", "Weight")

  # Remove interactions with weight equal 0
  if (remove_zeroes == TRUE) {
    net <- net[net$Weight != 0, ]
  }

  # Remove the weight column if weight is set to FALSE
  if (!weight) {
    net <- net[, -3]
  }

  # Reorder by Nodes 1 and 2
  net$Node1 <- factor(net$Node1, levels = rownames(mat))
  net$Node2 <- factor(net$Node2, levels = colnames(mat))
  net <- net[order(net$Node1, net$Node2), ]

  # Transform the two first columns in character
  net[, 1] <- as.character(net[, 1])
  net[, 2] <- as.character(net[, 2])

  # Return the data.frame
  return(net)

}
