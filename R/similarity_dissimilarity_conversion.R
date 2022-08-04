#' Convert similarity indices to dissimilarity metrics
#'
#' This function converts a data.frame of similarity metrics between sites to
#'  dissimilarity metrics (= beta diversity).
#'
#' @param similaritydata the output object from \code{\link{similarity}} or a
#' \code{data.frame} with the first columns called "Site1" and "Site2", and the
#' other columns being the similarity metrics.
#' @export
#' @note
#' \loadmathjax
#' The behaviour of this function changes depending on column names. Columns
#' "Site1" and "Site2" are copied identically. If there are columns called
#' "a", "b", "c", "A", "B", "C", they will also be copied identically.
#'
#' If a column is called "Euclidean", its distance will be calculated based
#' on the following formula:
#'
#' \mjeqn{Euclidean distance = (1 - Euclidean similarity) / Euclidean similarity}{Euclidean distance = (1 - Euclidean similarity) / Euclidean similarity}
#'
#' Otherwise, all other columns will be transformed into dissimilarity with the
#' following formula:
#'
#' \mjeqn{dissimilarity = 1 - similarity}{dissimilarity = 1 - similarity}
#'
#'
#' @return A \code{data.frame} with additional class \code{bioRgeo.pairwise.metric}
#'
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}),
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) and
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' @seealso \link{dissimilarity_to_similarity}
#' @examples
#' comat <- matrix(sample(0:1000, size = 50, replace = TRUE, prob = 1 / 1:1001), 5, 10)
#' rownames(comat) <- paste0("Site", 1:5)
#' colnames(comat) <- paste0("Species", 1:10)
#'
#' simil <- similarity(comat, metric = "all")
#' simil
#'
#' dissimilarity <- similarity_to_dissimilarity(simil)
#' dissimilarity
similarity_to_dissimilarity <- function(similaritydata) {

  if (!inherits(similaritydata, "bioRgeo.pairwise.metric")) {
    stop("similaritydata should be a bioRgeo object created by similarity() or dissimilarity_to_similarity()")
  }
  if (attr(similaritydata, "type") == "dissimilarity") {
    stop("similaritydata is already composed of dissimilarity indices. If you want to convert it to similarity, use dissimilarity_to_similarity()")
  }

  dissimilaritydata <- similaritydata

  # Overwrite attribute
  attr(dissimilaritydata, "type") <- "dissimilarity"


  metrics <- colnames(similaritydata)[-which(colnames(similaritydata) %in% c("Site1", "Site2", "a", "b", "c", "A", "B", "C"))]

  # Special case for Euclidean distances
  if ("Euclidean" %in% metrics) {
    dissimilaritydata[, "Euclidean"] <- (1 - similaritydata[, "Euclidean"]) / similaritydata[, "Euclidean"]
    metrics <- metrics[-which(metrics == "Euclidean")]
  }
  # If there are other metrics than Euclidean, we use the same formula for all of them
  if (length(metrics)) {
    dissimilaritydata[, metrics] <- 1 - similaritydata[, metrics]
  }

  return(dissimilaritydata)
}

#' Convert dissimilarity indices to similarity indices
#'
#' This function converts a data.frame of dissimilarity indices (beta diversity)
#' between sites to similarity indices.
#'
#' @param dissimilaritydata the output object from \code{\link{similarity_to_dissimilarity}} or a
#' \code{data.frame} with the first columns called "Site1" and "Site2", and the
#' other columns being the similarity metrics.
#' @export
#' @note
#' The behaviour of this function changes depending on column names. Columns
#' "Site1" and "Site2" are copied identically. If there are columns called
#' "a", "b", "c", "A", "B", "C", they will also be copied identically.
#'
#' If a column is called "Euclidean", its similarity will be calculated based
#' on the following formula:
#'
#' \eqn{Euclidean similarity = (1 - Euclidean distance) / Euclidean distance}
#'
#' Otherwise, all other columns will be transformed into similarity with the
#' following formula:
#'
#' \eqn{similarity = 1 - dissimilarity}
#'
#'
#' @return A \code{data.frame} with additional class \code{bioRgeo.similarity}
#'
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}),
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) and
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' @seealso \link{dissimilarity_to_similarity}
#' @examples
#' comat <- matrix(sample(0:1000, size = 50, replace = TRUE, prob = 1 / 1:1001), 5, 10)
#' rownames(comat) <- paste0("Site", 1:5)
#' colnames(comat) <- paste0("Species", 1:10)
#'
#' simil <- similarity(comat, metric = "all")
#' simil
#'
#' dissimilarity <- similarity_to_dissimilarity(simil)
#' dissimilarity
#'
#' simil <- dissimilarity_to_similarity(dissimilarity)
#' simil
dissimilarity_to_similarity <- function(dissimilaritydata) {
  if (!inherits(dissimilaritydata, "bioRgeo.pairwise.metric")) {
    stop("dissimilaritydata should be a bioRgeo object created by dissimilarity() or similarity_to_dissimilarity()")
  }
  if (attr(dissimilaritydata, "type") == "similarity") {
    stop("dissimilaritydata is already composed of similarity indices. If you want to convert it to dissimilarity, use similarity_to_dissimilarity()")

  }

  similaritydata <- dissimilaritydata

  # Overwrite attribute
  attr(similaritydata, "type") <- "dissimilarity"

  metrics <- colnames(similaritydata)[-which(colnames(similaritydata) %in% c("Site1", "Site2", "a", "b", "c", "A", "B", "C"))]
  # Special case for Euclidean distances
  if ("Euclidean" %in% metrics) {
    similaritydata[, "Euclidean"] <- 1 / (1 + dissimilaritydata[, "Euclidean"])
    metrics <- metrics[-which(metrics == "Euclidean")]
  }

  # If there are other metrics than Euclidean, we use the same formula for all fo them
  if (length(metrics)) {
    similaritydata[, metrics] <- 1 - similaritydata[, metrics]
  }

  return(similaritydata)
}
