#' Compute dissimilarity metrics (beta-diversity) between sites based on
#' species composition
#'
#' This function creates a `data.frame` where each row provides one or
#' several dissimilarity metric(s) between each pair of sites from a
#' co-occurrence `matrix` with sites as rows and species as columns.
#'
#' @param comat a co-occurrence `matrix` with sites as rows and species
#' as columns.
#' 
#' @param metric a `character` vector indicating which metrics to chose
#' (see Details). Available options are *abc*, *ABC*, *Jaccard*,
#' *Jaccardturn*, *Sorensen*, *Simpson*,  *Bray*,
#' *Brayturn* or *Euclidean*.\cr
#' If `"all"` is specified, then all metrics will be
#' calculated. Can be set to `NULL` if `formula` is used.
#' 
#' @param formula a `character` vector with your own formula(s) based on the
#' `a`, `b`, `c`, `A`, `B`, and `C` quantities
#' (see Details). `formula` is set to `NULL` by default.
#' 
#' @param method a `character` indicating what method should be used to compute
#' `abc` (see Details).
#' `method = "prodmat"` by default is more efficient but can be greedy
#' in memory and `method="loops"` is less efficient but less greedy in
#' memory.
#' 
#' @details
#' \loadmathjax
#' With `a` the number of species shared by a pair of sites, `b` species only
#' present in the first site  and `c` species only present in the second site.
#'
#' \mjeqn{Jaccard = (b + c) / (a + b + c)}{Jaccard = (b + c) / (a + b + c)}
#'
#' \mjeqn{Jaccardturn = 2min(b, c) / (a + 2min(b, c))}{Jaccardturn = 2min(b, c) / (a + 2min(b, c))}\insertCite{Baselga2012}{bioregion}
#'
#' \mjeqn{Sorensen = (b + c) / (2a + b + c)}{Sorensen = (b + c) / (2a + b + c)}
#'
#' \mjeqn{Simpson = min(b, c) / (a + min(b, c))}{Simpson = min(b, c) / (a + min(b, c))}
#'
#' If abundances data are available, Bray-Curtis and its turnover component
#' can also be computed with the following equation:
#'
#' \mjeqn{Bray = (B + C) / (2A + B + C)}{Bray = (B + C) / (2A + B + C)}
#'
#' \mjeqn{Brayturn = min(B, C)/(A + min(B, C))}{Brayturn = min(B, C)/(A + min(B, C))} \insertCite{Baselga2013}{bioregion}
#'
#' with A the sum of the lesser values for common species shared by a pair of
#' sites. B and C are the total number of specimens counted at both sites minus
#' A.
#'
#' `formula` can be used to compute customized metrics with the terms
#' `a`, `b`, `c`, `A`, `B`, and `C`. For example
#' `formula = c("pmin(b,c) / (a + pmin(b,c))", "(B + C) / (2*A + B + C)")`
#' will compute the Simpson and Bray-Curtis dissimilarity metrics, respectively. 
#' **Note that pmin is used in the Simpson formula because a, b, c, A, B and C 
#' are `numeric` vectors.**
#'
#' Euclidean computes the Euclidean distance between each pair of sites.
#'
#' @return 
#' 
#' A `data.frame` with additional class `bioregion.pairwise.metric`,
#' providing one or several dissimilarity
#' metric(s) between each pair of sites. The two first columns represent each
#' pair of sites.
#' One column per dissimilarity metric provided in `metric` and 
#' `formula` except for the metric *abc* and *ABC* that
#' are stored in three columns (one for each letter).
#' 
#' @seealso [similarity()] [dissimilarity_to_similarity] 
#' [similarity_to_dissimilarity]
#' 
#' @author
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}), 
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) and
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' 
#' @examples
#' comat <- matrix(sample(0:1000, size = 50, replace = TRUE,
#' prob = 1 / 1:1001), 5, 10)
#' rownames(comat) <- paste0("Site", 1:5)
#' colnames(comat) <- paste0("Species", 1:10)
#'
#' dissim <- dissimilarity(comat,
#' metric = c("abc", "ABC", "Simpson", "Brayturn"))
#'
#' dissim <- dissimilarity(comat, metric = "all",
#' formula = "1 - (b + c) / (a + b + c)")
#' 
#' @references
#' \insertRef{Baselga2012}{bioregion}
#' 
#' \insertRef{Baselga2013}{bioregion}
#' 
#' @export

dissimilarity <- function(comat, 
                          metric = "Simpson", 
                          formula = NULL,
                          method = "prodmat"){
  
# Compute similarities
  res <- similarity(comat, metric = metric, formula = formula, method = method)
  
  # Compute dissimilarity
  res <- similarity_to_dissimilarity(res, include_formula = FALSE)
  
  # Return the output
  return(res)
}
