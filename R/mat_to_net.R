#' Create a data.frame from a contingency table
#'
#' This function creates a two- or three-columns `data.frame` where
#' each row represents the interaction between two nodes (site and species for
#' example) and an optional third column indicating the weight of the
#' interaction (if `weight = TRUE`) from a contingency table (sites as
#' rows and species as columns for example).
#'
#' @param mat a contingency table (i.e. `matrix`).
#'
#' @param weight a `boolean` indicating if the value are weights.
#'
#' @param remove_zeroes a `boolean` determining whether interactions with
#' weight equal to 0 should be removed from the output.
#'
#' @param include_diag a `boolean` indicating whether the diagonal
#' should be included in the output. Only for squared matrix.
#'
#' @param include_lower a `boolean` indicating whether the lower triangular
#' matrix should be included in the output. Only for squared matrix.
#'
#' @return A `data.frame` where each row represents the interaction
#' between two nodes and an optional third column indicating the weight of the
#' interaction.
#'
#' @author
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) \cr
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) \cr
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' 
#' @seealso [net_to_mat]
#'
#' @examples
#' mat <- matrix(sample(1000, 50), 5, 10)
#' rownames(mat) <- paste0("Site", 1:5)
#' colnames(mat) <- paste0("Species", 1:10)
#'
#' net <- mat_to_net(mat, weight = TRUE)
#'
#' @importFrom tidyr pivot_longer
#'
#' @export

mat_to_net <- function(mat,
                       weight = FALSE,
                       remove_zeroes = TRUE,
                       include_diag = TRUE,
                       include_lower = TRUE) {
  
  # Control inputs
  controls(args = weight, data = NULL, type = "boolean")
  controls(args = remove_zeroes, data = NULL, type = "boolean")
  controls(args = include_diag, data = NULL, type = "boolean")
  controls(args = include_lower, data = NULL, type = "boolean")
  
  controls(args = NULL, data = mat, type = "input_matrix")
  if (dim(mat)[1] == dim(mat)[2]) { # Squared matrix
    if (!include_diag) {
      diag(mat) <- NA
    }
    if (!include_lower) {
      mat[lower.tri(mat)] <- NA
    }
  } else {
    if (include_diag == FALSE) {
      message("include_diag is only used with squared matrix.")
    }
    if (include_lower == FALSE) {
      message("include_lower is only used with squared matrix.")
    }
  }

  # Visible binding for global variable
  Node1 <- NULL

  # Conversion as data.frame
  mat <- as.data.frame(mat)
  mat$Node1 <- rownames(mat)
  net <- as.data.frame(tidyr::pivot_longer(
    data = as.data.frame(mat),
    cols = -Node1, names_to = "Node2",
    values_to = "Weight"
  ))

  # Remove interactions with weight equal 0
  if(!include_diag | !include_lower){
    net <- net[!is.na(net$Weight),]
  }
  if (remove_zeroes) {
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
