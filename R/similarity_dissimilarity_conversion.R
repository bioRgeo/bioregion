#' Convert similarity metrics to dissimilarity metrics
#'
#' This function converts a `data.frame` of similarity metrics between sites to
#'  dissimilarity metrics (beta diversity).
#'
#' @param similarity the output object from [similarity()] or
#' [dissimilarity_to_similarity()].
#' 
#' @param include_formula a `boolean` indicating if the metrics based on your 
#' own formula(s) should be converted (see Details). This argument is set to 
#' `TRUE` by default.
#' 
#' @note
#' The behavior of this function changes depending on column names. Columns
#' `Site1` and `Site2` are copied identically. If there are columns called
#' `a`, `b`, `c`, `A`, `B`, `C` they will also be copied identically. If there 
#' are columns based on your own formula (argument `formula` in [similarity()])
#' or not in the original list of similarity metrics (argument `metrics` in 
#' [similarity()]) and if the argument `include_formula` is set to `FALSE`, 
#' they will also be copied identically. Otherwise there are going to be
#' converted like they other columns (default behavior).
#'
#' If a column is called `Euclidean`, its distance will be calculated based
#' on the following formula:
#'
#' Euclidean distance = (1 - Euclidean similarity) / Euclidean similarity
#'
#' Otherwise, all other columns will be transformed into dissimilarity with the
#' following formula:
#'
#' dissimilarity = 1 - similarity
#'
#' @return A `data.frame` with additional class 
#' `bioregion.pairwise.metric`, providing dissimilarity
#' metric(s) between each pair of sites based on a similarity object.
#'
#' @author
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) \cr
#' Boris Leroy (\email{leroy.boris@gmail.com}) \cr
#' Pierre Denelle (\email{pierre.denelle@gmail.com})
#' 
#' @seealso [dissimilarity_to_similarity()] [similarity()] [dissimilarity()]
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
similarity_to_dissimilarity <- function(similarity, include_formula = TRUE){
  
  # List metrics to change
  all <- c("a", "b", "c","A", "B", "C",
           "Jaccard", "Jaccardturn", "Sorensen", "Simpson",
           "Bray", "Brayturn",
           "Euclidean")
  noteucl <- c("Jaccard", "Jaccardturn", "Sorensen", "Simpson",
               "Bray", "Brayturn")
  eucl <- "Euclidean"
  
  
  # Initialize output
  output <- similarity
  
  # Controls
  controls(args = include_formula, data = NULL, type = "boolean")
  controls(args = NULL, data = similarity, type = "input_conversion_similarity")
  
  # Overwrite attribute
  attr(output, "type") <- "dissimilarity"
  
  # Identify columns
  metrics <- colnames(output)[-c(1,2)]
  
  # Euclidean
  poseucl <- which(metrics %in% eucl)
  if(length(poseucl) > 0){
    output[,(poseucl + 2)] = (1 - output[,(poseucl + 2)]) / 
      output[,(poseucl + 2)]
  }
  
  # Not Euclidean
  posnoteucl <- which(metrics %in% noteucl)
  if(length(posnoteucl) > 0){
    output[,(posnoteucl + 2)] =  1- output[,(posnoteucl + 2)]  
  }
  
  # Include formula ?
  if(include_formula){
    posnotall <- which(!(metrics %in% all))
    if(length(posnotall) > 0){
      output[,(posnotall + 2)] = 1 - output[,(posnotall + 2)]  
    }
  }

  # Return output
  return(output)
}

#' Convert dissimilarity metrics to similarity metrics
#'
#' This function converts a `data.frame` of dissimilarity metrics (beta diversity)
#' between sites to similarity metrics.
#'
#' @param dissimilarity the output object from [dissimilarity()] or
#' [similarity_to_dissimilarity()].
#' 
#' @param include_formula a `boolean` indicating if the metrics based on your 
#' own formula(s) should be converted (see Details). This argument is set to 
#' `TRUE` by default.
#' 
#' @note
#' The behavior of this function changes depending on column names. Columns
#' `Site1` and `Site2` are copied identically. If there are columns called
#' `a`, `b`, `c`, `A`, `B`, `C` they will also be copied identically. If there 
#' are columns based on your own formula (argument `formula` in 
#' [dissimilarity()]) or not in the original list of dissimilarity metrics 
#' (argument `metrics` in [dissimilarity()]) and if the argument 
#' `include_formula` is set to `FALSE`, they will also be copied identically. 
#' Otherwise there are going to be converted like they other columns (default 
#' behavior).
#'
#' If a column is called `Euclidean`, the similarity will be calculated based
#' on the following formula:
#'
#' Euclidean similarity = 1 / (1 - Euclidean distance)
#'
#' Otherwise, all other columns will be transformed into dissimilarity with the
#' following formula:
#'
#'similarity = 1 - dissimilarity
#'
#' @return A `data.frame` with additional class 
#' `bioregion.pairwise.metric`, providing similarity
#' metric(s) between each pair of sites based on a dissimilarity object.
#'
#' @author
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) \cr
#' Boris Leroy (\email{leroy.boris@gmail.com}) \cr
#' Pierre Denelle (\email{pierre.denelle@gmail.com})
#' 
#' @seealso [similarity_to_dissimilarity()] [similarity()] [dissimilarity()]
#' 
#' @examples
#' comat <- matrix(sample(0:1000, size = 50, replace = TRUE,
#' prob = 1 / 1:1001), 5, 10)
#' rownames(comat) <- paste0("Site", 1:5)
#' colnames(comat) <- paste0("Species", 1:10)
#'
#' dissimil <- dissimilarity(comat, metric = "all")
#' dissimil
#'
#' similarity <- dissimilarity_to_similarity(dissimil)
#' similarity
#' 
#' @export
dissimilarity_to_similarity <- function(dissimilarity, include_formula = TRUE){
  
  # List metrics to change
  all <- c("a", "b", "c","A", "B", "C",
           "Jaccard", "Jaccardturn", "Sorensen", "Simpson",
           "Bray", "Brayturn",
           "Euclidean")
  noteucl <- c("Jaccard", "Jaccardturn", "Sorensen", "Simpson",
               "Bray", "Brayturn")
  eucl <- "Euclidean"
  
  
  # Initialize output
  output <- dissimilarity
  
  # Controls
  controls(args = include_formula, data = NULL, type = "boolean")
  controls(args = NULL, data = dissimilarity, type = "input_conversion_dissimilarity")
  
  # Overwrite attribute
  attr(output, "type") <- "similarity"
  
  # Identify columns
  metrics <- colnames(output)[-c(1,2)]
  
  # Euclidean
  poseucl <- which(metrics %in% eucl)
  if(length(poseucl) > 0){
    output[,(poseucl + 2)] = 1 / (1 + output[,(poseucl + 2)])
  }
  
  # Not Euclidean
  posnoteucl <- which(metrics %in% noteucl)
  if(length(posnoteucl) > 0){
    output[,(posnoteucl + 2)] = 1- output[,(posnoteucl + 2)]  
  }
  
  # Include formula ?
  if(include_formula){
    posnotall <- which(!(metrics %in% all))
    if(length(posnotall) > 0){
      output[,(posnotall + 2)] = 1- output[,(posnotall + 2)]  
    }
  }
  
  # Return output
  return(output)
}
