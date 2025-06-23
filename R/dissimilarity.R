#' Compute dissimilarity metrics (beta-diversity) between sites based on
#' species composition
#' 
#' This function generates a `data.frame` where each row provides one or
#' several dissimilarity metrics between pairs of sites, based on a 
#' co-occurrence `matrix` with sites as rows and species as columns.
#'
#' @param comat A co-occurrence `matrix` with sites as rows and species as
#' columns.
#' 
#' @param metric A `character` vector or a single `character` string specifying
#' the metrics to compute (see Details). Available options are `"abc"`, `"ABC"`,
#' `"Jaccard"`, `"Jaccardturn"`, `"Sorensen"`, `"Simpson"`, `"Bray"`,
#' `"Brayturn"`, and `"Euclidean"`. If `"all"` is specified, all metrics will
#' be calculated. Can be set to `NULL` if `formula` is used.
#' 
#' @param formula A `character` vector or a single `character` string specifying 
#' custom formula(s) based on the `a`, `b`, `c`, `A`, `B`, and `C` quantities 
#' (see Details). The default is `NULL`.
#' 
#' @param method A `character` string specifying the method to compute `abc` 
#' (see Details). The default is `"prodmat"`, which is more efficient but 
#' memory-intensive. Alternatively, `"loops"` is less memory-intensive but 
#' slower.
#' 
#' @return 
#' A `data.frame` with the additional class 
#' `bioregion.pairwise.metric`, containing one or several dissimilarity
#' metrics between pairs of sites. The first two columns represent the pairs of 
#' sites. There is one column per similarity metric provided in `metric` and
#' `formula`, except for the `abc` and `ABC` metrics, which are stored in three 
#' separate columns (one for each letter).
#' 
#' @details
#' With `a` the number of species shared by a pair of sites, `b` species only
#' present in the first site  and `c` species only present in the second site.
#'
#' Jaccard = (b + c) / (a + b + c)
#'
#' Jaccardturn = 2min(b, c) / (a + 2min(b, c)) (Baselga, 2012)
#'
#' Sorensen = (b + c) / (2a + b + c)
#'
#' Simpson = min(b, c) / (a + min(b, c))
#'
#' If abundances data are available, Bray-Curtis and its turnover component
#' can also be computed with the following equation:
#'
#' Bray = (B + C) / (2A + B + C)
#' 
#' Brayturn = min(B, C)/(A + min(B, C)) (Baselga, 2013)
#'
#' with `A` the sum of the lesser values for common species shared by a pair of
#' sites. `B` and `C` are the total number of specimens counted at both sites 
#' minus `A`.
#'
#' `formula` can be used to compute customized metrics with the terms
#' `a`, `b`, `c`, `A`, `B`, and `C`. For example
#' `formula = c("pmin(b,c) / (a + pmin(b,c))", "(B + C) / (2*A + B + C)")`
#' will compute the Simpson and Bray-Curtis dissimilarity metrics, respectively. 
#' Note that `pmin` is used in the Simpson formula because `a`, `b`, `c`, `A`, 
#' `B` and `C` are `numeric` vectors.
#'
#' Euclidean computes the Euclidean distance between each pair of sites.
#' 
#' @seealso 
#' For more details illustrated with a practical example, 
#' see the vignette: 
#' \url{https://biorgeo.github.io/bioregion/articles/a3_pairwise_metrics.html}.
#' 
#' Associated functions: 
#' [similarity] [dissimilarity_to_similarity] 
#' 
#' @author
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) \cr
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) \cr
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' 
#' @examples
#' comat <- matrix(sample(0:1000, size = 50, replace = TRUE,
#' prob = 1 / 1:1001), 5, 10)
#' rownames(comat) <- paste0("s", 1:5)
#' colnames(comat) <- paste0("sp", 1:10)
#'
#' dissim <- dissimilarity(comat,
#' metric = c("abc", "ABC", "Simpson", "Brayturn"))
#'
#' dissim <- dissimilarity(comat, metric = "all",
#' formula = "1 - (b + c) / (a + b + c)")
#' 
#' @references
#' Baselga, A. (2012) The Relationship between Species Replacement,
#' Dissimilarity Derived from Nestedness, and Nestedness.
#' \emph{Global Ecology and Biogeography}, 21(12), 1223--1232.
#' 
#' Baselga, A. (2013) Separating the two components of abundance-based
#' dissimilarity: balanced changes in abundance vs. abundance gradients.
#' \emph{Methods in Ecology and Evolution}, 4(6), 552--557.
#' 
#' @export
dissimilarity <- function(comat, 
                          metric = "Simpson", 
                          formula = NULL,
                          method = "prodmat"){
  
# Compute similarities
  res <- similarity(comat, 
                    metric = metric, 
                    formula = formula, 
                    method = method)
  
  # Compute dissimilarity
  res <- similarity_to_dissimilarity(res, include_formula = FALSE)
  
  # Return the output
  return(res)
}
