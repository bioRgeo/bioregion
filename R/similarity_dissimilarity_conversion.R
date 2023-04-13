#' Convert similarity indices to dissimilarity metrics
#'
#' This function converts a data.frame of similarity metrics between sites to
#'  dissimilarity metrics (= beta diversity).
#'
#' @param similarity the output object from [similarity()] or a
#' `data.frame` with the first columns called "Site1" and "Site2", and the
#' other columns being the similarity metrics.
#' 
#' @note
#' \loadmathjax
#' The behavior of this function changes depending on column names. Columns
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
#' @return A `data.frame` with additional class 
#' `bioregion.pairwise.metric`, providing one or several similarity
#' metric(s) between each pair of sites. The two first columns represent each 
#' pair of sites, and the other column represent similarity metrics. Columns
#' with names "a", "b", "c", "A", "B" and "C"  are not altered.
#'
#' @author
#' Boris Leroy (\email{leroy.boris@gmail.com}),
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) and
#' Pierre Denelle (\email{pierre.denelle@gmail.com})
#' 
#' @seealso [dissimilarity_to_similarity()]
#' 
#' @examples
#' comat <- matrix(sample(0:1000, size = 50, replace = TRUE,
#' prob = 1 / 1:1001), 5, 10)
#' rownames(comat) <- paste0("Site", 1:5)
#' colnames(comat) <- paste0("Species", 1:10)
#'
#' simil <- similarity(comat, metric = "all")
#' simil
#'
#' dissimilarity <- similarity_to_dissimilarity(simil)
#' dissimilarity
#' 
#' @export

similarity_to_dissimilarity <- function(similarity){
  
  if (!inherits(similarity, "bioregion.pairwise.metric")) {
    stop("similarity should be a bioregion object created by similarity() or
         dissimilarity_to_similarity()")
  }
  if (attr(similarity, "type") == "dissimilarity") {
    stop("similarity is already composed of dissimilarity indices. If you want
         to convert it to similarity, use dissimilarity_to_similarity()")
  }
  
  dissimilaritydata <- similarity
  
  # Overwrite attribute
  attr(dissimilaritydata, "type") <- "dissimilarity"
  
  metrics <- colnames(similarity)[-which(colnames(similarity) %in% 
                                           c("Site1", "Site2", "a", "b", "c",
                                             "A", "B", "C"))]
  
  # Special case for Euclidean distances
  if ("Euclidean" %in% metrics) {
    dissimilaritydata[, "Euclidean"] <- (1 - similarity[, "Euclidean"]) /
      similarity[, "Euclidean"]
    metrics <- metrics[-which(metrics == "Euclidean")]
  }
  # If there are other metrics than Euclidean, we use the same formula for all
  # of them
  if (length(metrics)) {
    dissimilaritydata[, metrics] <- 1 - similarity[, metrics]
  }
  
  return(dissimilaritydata)
}

#' Convert dissimilarity indices to similarity indices
#'
#' This function converts a data.frame of dissimilarity indices 
#' (beta diversity) between sites to similarity indices.
#'
#' @param dissimilaritydata the output object from
#' [similarity_to_dissimilarity()] or a `data.frame` with the first columns
#' called "Site1" and "Site2", and the other columns being the similarity
#' metrics.
#' 
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
#' @return A `data.frame` with additional class `bioregion.pairwise.metric`,
#' providing one or several dissimilarity metric(s) between each pair of sites.
#' The two first columns represent each pair of sites, and the other column
#' represent dissimilarity metrics. Columns with names "a", "b", "c", "A", "B"
#' and "C"  are not altered.
#' 
#' @author
#' Boris Leroy (\email{leroy.boris@gmail.com}),
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) and
#' Pierre Denelle (\email{pierre.denelle@gmail.com})
#' 
#' @seealso [dissimilarity_to_similarity()]
#' 
#' @examples
#' comat <- matrix(sample(0:1000, size = 50, replace = TRUE,
#' prob = 1 / 1:1001), 5, 10)
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
#' 
#' @export

dissimilarity_to_similarity <- function(dissimilaritydata) {
  if (!inherits(dissimilaritydata, "bioregion.pairwise.metric")) {
    stop("dissimilaritydata should be a bioregion object created by
         dissimilarity() or similarity_to_dissimilarity()")
  }
  if (attr(dissimilaritydata, "type") == "similarity") {
    stop("dissimilaritydata is already composed of similarity indices. If you
         want to convert it to dissimilarity, use
         similarity_to_dissimilarity()")
  }
  
  similarity <- dissimilaritydata
  
  # Overwrite attribute
  attr(similarity, "type") <- "dissimilarity"
  
  metrics <- colnames(similarity)[-which(colnames(similarity) %in%
                                           c("Site1", "Site2", "a", "b", "c",
                                             "A", "B", "C"))]
  # Special case for Euclidean distances
  if ("Euclidean" %in% metrics) {
    similarity[, "Euclidean"] <- 1 / (1 + dissimilaritydata[, "Euclidean"])
    metrics <- metrics[-which(metrics == "Euclidean")]
  }
  
  # If there are other metrics than Euclidean, we use the same formula for all
  # of them
  if (length(metrics)) {
    similarity[, metrics] <- 1 - similarity[, metrics]
  }
  
  return(similarity)
}
