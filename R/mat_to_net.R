#' Create a data.frame from a contingency table
#'
#' This function creates a two- or three-columns \code{data.frame} where
#' each row represents the interaction between two objects (site and species for example)
#' and an optional third column indicating the weight of the interaction (if \code{weight = TRUE})
#' from a contingency table (sites as rows and species as columns for example).
#'
#' @param mat a contingency table (i.e. \code{matrix})
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
#' mat <- matrix(sample(1000, 50), 5, 10)
#' rownames(mat) <- paste0("Site", 1:5)
#' colnames(mat) <- paste0("Species", 1:10)
#'
#' net <- mat_to_net(mat, weight = TRUE)
#' @export
mat_to_net <- function(mat, weight = FALSE, remove_absent_objects = TRUE) {

  # Controls
  if (!is.matrix(mat)) {
    stop("Contingency table should be a matrix")
  }

  sco <- sum(is.na(mat))
  if (sco > 0) {
    stop("NA(s) detected in the contingency table")
  }

  if (!is.logical(weight)) {
    stop("weight must be a boolean")
  }

  if (!is.logical(remove_absent_objects)) {
    stop("remove_absent_objects must be a boolean")
  }

  # Conversion as data.frame
  net <- reshape2::melt(mat)
  colnames(net) <- c("Object1", "Object2", "Weight")

  # Remove interactions with weight equal 0
  if (remove_absent_objects == TRUE) {
    net <- net[net$Weight != 0, ]
  }

  # Remove the weight column if weight is set to FALSE
  if (!weight) {
    net <- net[, -3]
  }

  # Reorder by Object 1 and 2
  net$Object1 <- factor(net$Object1, levels = rownames(mat))
  net$Object2 <- factor(net$Object2, levels = colnames(mat))
  net <- net[order(net$Object1, net$Object2), ]

  # Transform the two first columns in character
  net[, 1] <- as.character(net[, 1])
  net[, 2] <- as.character(net[, 2])

  # Return the data.frame
  return(net)

}
