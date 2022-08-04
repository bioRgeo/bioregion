#' Compute dissimilarity metrics (beta-diversity) between sites based on species composition
#'
#' This function creates a \code{data.frame} where each row provides one or several dissimilarity
#' metric(s) between each pair of sites from a co-occurence \code{matrix} with sites as rows and species as columns.
#'
#' @param comat a co-occurence \code{matrix} with sites as rows and species as columns.
#' @param metric a vector of string(s) indicating which metrics
#' metric(s) to chose (see Details). If \code{"all"} is specified, then all
#' metrics will be calculated. Can be set to \code{NULL} if \code{formula} is used.
#' @param formula a vector of string(s) with your own formula based on the \code{a}, \code{b}, \code{c}, \code{A}, \code{B},
#' and \code{C} quantities (see Details). \code{formula} is set to \code{NULL} by default.
#' @param method a string indicating what method should be used to compute \code{abc} (see Details).
#' \code{method = "prodmat"} by default is more efficient but can be greedy in memory and \code{method="loops"} is less efficient
#' but less greedy in memory.
#' @export
#' @details
#' \loadmathjax
#' With \code{a} the number of species shared by a pair of sites, \code{b} species only present in the first site
#' and \code{c} species only present in the second site.
#'
#' \mjeqn{Jaccard = (b + c) / (a + b + c)}{Jaccard = 1 - (b + c) / (a + b + c)}
#'
#' \mjeqn{Jaccardturn = 2min(b, c) / (a + 2min(b, c))}{Jaccardturn = 1 - 2min(b, c) / (a + 2min(b, c))}\insertCite{Baselga2012}{bioRgeo}
#'
#' \mjeqn{Sorensen = (b + c) / (2a + b + c)}{Sorensen = 1 - (b + c) / (2a + b + c)}
#'
#' \mjeqn{Simpson = min(b, c) / (a + min(b, c))}{Simpson = 1 - min(b, c) / (a + min(b, c))}
#'
#' If abundances data are available, Bray-Curtis and its turnover component can also be computed with the
#' following equation:
#'
#' \mjeqn{Bray = (B + C) / (2A + B + C)}{Bray = 1 - (B + C) / (2A + B + C)}
#'
#' \mjeqn{Brayturn = min(B, C)/(A + min(B, C))}{Brayturn = 1 - min(B, C)/(A + min(B, C))} \insertCite{Baselga2013}{bioRgeo}
#'
#' with A the sum of the lesser values for common species shared by a pair of sites.
#' B and C are the total number of specimens counted at both sites minus A.
#'
#' \code{formula} can be used to compute customized metrics with the terms \code{a}, \code{b}, \code{c}, \code{A}, \code{B},
#' and \code{C}. For example \code{formula = c("(b + c) / (a + b + c)", "(B + C) / (2*A + B + C)")} will
#' compute the Jaccard and Bray-Curtis dissimilarity metrics, respectively.
#'
#' Euclidean computes the Euclidean distance between each pair of sites.
#'
#' @return A \code{data.frame} providing one or several dissimilarity
#' metric(s) between each pair of sites. The two first columns represents each pair of sites.
#' One column per dissimilarity metric provided in \code{metric} and \code{formula} except for the metric \emph{abc} and \emph{ABC} that
#' are stored in three columns (one for each letter).
#' @seealso \link{similarity} \link{dissimilarity_to_similarity} \link{similarity_to_dissimilarity}
#' @author
#' Pierre Denelle (\email{pierre.denelle@gmail.com}),
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) and
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' @examples
#' comat <- matrix(sample(0:1000, size = 50, replace = TRUE, prob = 1 / 1:1001), 5, 10)
#' rownames(comat) <- paste0("Site", 1:5)
#' colnames(comat) <- paste0("Species", 1:10)
#'
#' dist <- dissimilarity(comat, metric = c("abc", "ABC", "Simpson", "Brayturn"))
#' dist
#'
#' simil <- dissimilarity(comat, metric = "all", formula = "1 - (b + c) / (a + b + c)")
#' dist
#' @references
#' \insertRef{Baselga2012}{bioRgeo}
#' 
#' \insertRef{Baselga2013}{bioRgeo}
#' @export
dissimilarity <- function(comat, metric = "Simpson", formula = NULL, method = "prodmat") {

  # Compute similarities
  res <- similarity(comat, metric = metric, formula = formula, method = method)

  # Compute dissimilarity
  res <- similarity_to_dissimilarity(res)

  # Return the output
  return(res)
}