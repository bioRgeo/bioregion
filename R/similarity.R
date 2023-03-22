#' Compute similarity metrics between sites based on species composition
#'
#' This function creates a `data.frame` where each row provides one or
#' several similarity metric(s) between each pair of sites from a co-occurrence
#' `matrix` with sites as rows and species as columns.
#'
#' @param comat a co-occurrence `matrix` with sites as rows and species as
#' columns.
#' 
#' @param metric a vector of string(s) indicating which similarity metric(s) to
#' chose (see Details). If `"all"` is specified, then all metrics will be
#' calculated. Can be set to `NULL` if `formula` is used.
#' 
#' @param formula a vector of string(s) with your own formula based on the
#' `a`, `b`, `c`, `A`, `B`, and `C` quantities
#' (see Details). `formula` is set to `NULL` by default.
#' 
#' @param method a string indicating what method should be used to compute
#' `abc` (see Details).
#' `method = "prodmat"` by default is more efficient but can be greedy in
#' memory and `method="loops"` is less efficient but less greedy in
#' memory.
#' 
#' @details
#' \loadmathjax
#' With `a` the number of species shared by a pair of sites, `b`
#' species only present in the first site and `c` species only present in
#' the second site.
#'
#' \mjeqn{Jaccard = 1 - (b + c) / (a + b + c)}{Jaccard = 1 - (b + c) / (a + b + c)}
#'
#' \mjeqn{Jaccardturn = 1 - 2min(b, c) / (a + 2min(b, c))}{Jaccardturn = 1 - 2min(b, c) / (a + 2min(b, c))} \insertCite{Baselga2012}{bioregion}
#'
#' \mjeqn{Sorensen = 1 - (b + c) / (2a + b + c)}{Sorensen = 1 - (b + c) / (2a + b + c)}
#'
#' \mjeqn{Simpson = 1 - min(b, c) / (a + min(b, c))}{Simpson = 1 - min(b, c) / (a + min(b, c))}
#'
#' If abundances data are available, Bray-Curtis and its turnover component can
#' also be computed with the following equation:
#'
#' \mjeqn{Bray = 1 - (B + C) / (2A + B + C)}{Bray = 1 - (B + C) / (2A + B + C)}
#'
#' \mjeqn{Brayturn = 1 - min(B, C)/(A + min(B, C))}{Brayturn = 1 - min(B, C)/(A + min(B, C))} \insertCite{Baselga2013}{bioregion}
#'
#' with A the sum of the lesser values for common species shared by a pair of
#' sites.
#' B and C are the total number of specimens counted at both sites minus A.
#'
#' `formula` can be used to compute customized metrics with the terms
#' `a`, `b`, `c`, `A`, `B`, and `C`. For example
#' `formula = c("1 - (b + c) / (a + b + c)", "1 - (B + C) / (2*A + B + C)")`
#' will compute the Jaccard and Bray-Curtis similarity metrics, respectively.
#'
#' Euclidean computes the Euclidean similarity between each pair of site
#' following this equation:
#'
#' \mjeqn{Euclidean = 1 / (1 + d_{ij})}{Euclidean = 1 / (1 + d_{ij})}
#'
#' Where \mjeqn{d_{ij}}{d_{ij}} is the Euclidean distance between site i and 
#' site j in terms of species composition.
#'
#' @return A `data.frame` with additional class 
#' `bioregion.pairwise.metric`, providing one or several similarity
#' metric(s) between each pair of sites. The two first columns represent each 
#' pair of sites.
#' One column per similarity metric provided in `metric` and
#' `formula` except for the metric *abc* and *ABC* that are
#' stored in three columns (one for each letter).
#' 
#' @seealso [dissimilarity] [dissimilarity_to_similarity] [similarity_to_dissimilarity]
#' 
#' @author
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}), 
#' Pierre Denelle (\email{pierre.denelle@gmail.com}) and
#' Boris Leroy (\email{leroy.boris@gmail.com})
#' 
#' @examples
#' \donttest{
#' comat <- matrix(sample(0:1000, size = 50, replace = TRUE,
#' prob = 1 / 1:1001), 5, 10)
#' rownames(comat) <- paste0("Site", 1:5)
#' colnames(comat) <- paste0("Species", 1:10)
#'
#' simil <- similarity(comat, metric = c("abc", "ABC", "Simpson", "Brayturn"))
#'
#' simil <- similarity(comat, metric = "all",
#' formula = "1 - (b + c) / (a + b + c)")
#' }
#' @references
#' \insertRef{Baselga2012}{bioregion}
#' 
#' \insertRef{Baselga2013}{bioregion}
#' 
#' @importFrom Matrix tcrossprod
#' @importFrom stats dist
#' 
#' @export

similarity <- function(comat, metric = "Simpson", formula = NULL,
                       method = "prodmat"){
  
  # list of metrics based on abc
  lsmetricabc <- c("abc", "Jaccard", "Jaccardturn", "Sorensen", "Simpson")
  
  # list of metrics based on ABC
  lsmetricABC <- c("ABC", "Bray", "Brayturn")
  
  # list of metrics based on other features
  lsmetrico <- c("Euclidean")
  
  if ("all" %in% metric) {
    metric <- c(lsmetricabc, lsmetricABC, lsmetrico)
  }
  
  # Controls
  if (is.null(metric) & is.null(formula)) {
    stop("metric or formula should be used")
  }
  
  if (length(intersect(c(lsmetricabc, lsmetricABC, lsmetrico), metric)) !=
      length(metric) & is.null(formula)) {
    stop("One or several similarity metric(s) chosen is not available.
     Please chose among the followings:
         abc, Jaccard, Jaccardturn, Sorensen, Simpson, ABC, Bray, Brayturn or
         Euclidean.")
  }
  
  if (!is.null(formula) & !is.character(formula)) {
    stop("formula should be a vector of characters if not NULL")
  } else { # Check if abc and ABC in formula
    abcinformula <- (sum(c(grepl("a", formula), grepl("b", formula),
                           grepl("c", formula))) > 0)
    ABCinformula <- (sum(c(grepl("A", formula), grepl("B", formula),
                           grepl("C", formula))) > 0)
  }
  
  if (!is.matrix(comat)) {
    stop("Co-occurrence matrix should be a matrix")
  }
  
  sco <- sum(is.na(comat))
  minco <- min(comat)
  if (sco > 0) {
    stop("Co-occurrence matrix should contains only positive real: NA(s)
         detected!")
  }
  if (minco < 0) {
    stop("Co-occurrence matrix should contains only positive real: negative
         value detected!")
  }
  
  if (!(method %in% c("prodmat", "loops"))) {
    stop("The method is not available.
     Please chose among the followings: prodmat, loops")
  }
  
  # Extract site id
  siteid <- rownames(comat)
  
  # Initialize output
  res <- NULL
  
  # abcp: compute abc for presence data
  testabc <- ((length(intersect(lsmetricabc, metric)) > 0) |
                abcinformula) # check if abc should be computed
  if (testabc) {
    comatp <- comat
    comatp[comatp != 0] <- 1
    if (method == "prodmat") {
      # Compute the number of species in common "a" with matrix multiplication 
      # comatp%*%t(comatp)
      sumrow <- apply(comatp, 1, sum)
      # abcp <- prodmat(comatp,t(comatp))
      abcp <- Matrix::tcrossprod(comatp)
      rownames(abcp) <- siteid
      colnames(abcp) <- siteid
      
      # Create a data.frame from the matrix with mat_to_net (little trick to
      # deal with 0s)
      abcp[abcp == 0] <- -1
      abcp[lower.tri(abcp, diag = TRUE)] <- 0
      abcp <- mat_to_net(abcp, weight = TRUE, remove_zeroes = TRUE)
      colnames(abcp) <- c("Site1", "Site2", "a")
      abcp[abcp[, 3] == -1, 3] <- 0
      
      # Compute b and c
      abcp$b <- 0
      abcp$c <- 0
      abcp[, 4] <- sumrow[match(abcp[, 1], siteid)] - abcp[, 3]
      abcp[, 5] <- sumrow[match(abcp[, 2], siteid)] - abcp[, 3]
    }
    if (method == "loops") {
      # Use abc Rcpp function in src (three loops)
      abcp <- abc(comatp)
      abcp <- data.frame(Site1 = siteid[abcp[, 1]], Site2 = siteid[abcp[, 2]],
                         a = abcp[, 3], b = abcp[, 4], c = abcp[, 5])
    }
    
    # Update res if NULL
    if (is.null(res)) {
      res <- abcp[, 1:2]
    }
    
    # Compute metrics based on abc
    if ("Jaccard" %in% metric) {
      res$Jaccard <- 1 - (abcp$b + abcp$c) / (abcp$a + abcp$b + abcp$c)
    }
    if ("Jaccardturn" %in% metric) {
      res$Jaccardturn <- 1 - 2 * pmin(abcp$b, abcp$c) /
        (abcp$a + 2 * pmin(abcp$b, abcp$c))
    }
    if ("Sorensen" %in% metric) {
      res$Sorensen <- 1 - (abcp$b + abcp$c) / (2 * abcp$a + abcp$b + abcp$c)
    }
    if ("Simpson" %in% metric) {
      res$Simpson <- 1 - pmin(abcp$b, abcp$c) / (abcp$a + pmin(abcp$b, abcp$c))
    }
    
    # Attach abcp if abc in formula
    if (abcinformula) {
      a <- abcp$a
      b <- abcp$b
      c <- abcp$c
    }
  }
  
  # abca: compute ABC for abundance data
  testABC <- ((length(intersect(lsmetricABC, metric)) > 0) |
                ABCinformula) # check if ABC should be computed
  if (testABC) {
    # Use abc Rcpp function in src (three loops)
    abca <- abc(comat)
    abca <- data.frame(Site1 = siteid[abca[, 1]], Site2 = siteid[abca[, 2]],
                       A = abca[, 3], B = abca[, 4], C = abca[, 5])
    
    # Update res if NULL
    if (is.null(res)) {
      res <- abca[, 1:2]
    }
    
    # Compute metrics based on ABC
    if ("Bray" %in% metric) {
      res$Bray <- 1 - (abca$B + abca$C) / (2 * abca$A + abca$B + abca$C)
    }
    if ("Brayturn" %in% metric) {
      res$Brayturn <- 1 - pmin(abca$B, abca$C) /
        (abca$A + pmin(abca$B, abca$C))
    }
    
    # Attach abca if ABC in formula
    if (ABCinformula) {
      A <- abca$A
      B <- abca$B
      C <- abca$C
    }
  }
  
  # Compute equation in formula
  if (!is.null(formula)) {
    for (k in 1:length(formula)) {
      res <- cbind(res, eval(parse(text = formula[k])))
      colnames(res)[dim(res)[2]] <- formula[k]
    }
  }
  
  # Compute Euclidean similarity between site using dist()
  if ("Euclidean" %in% metric) {
    eucl <- as.matrix(stats::dist(comat))
    rownames(eucl) <- siteid
    colnames(eucl) <- siteid
    eucl[eucl == 0] <- -1
    eucl[lower.tri(eucl, diag = TRUE)] <- 0
    eucl <- mat_to_net(eucl, weight = TRUE, remove_zeroes = TRUE)
    colnames(eucl) <- c("Site1", "Site2", "Euclidean")
    eucl[eucl[, 3] == -1, 3] <- 0
    
    # Update res if NULL
    if (is.null(res)) {
      res <- eucl[, 1:2]
    }
    
    res$Euclidean <- 1 / (1 + eucl$Euclidean)
  }
  
  # Compute metrics based on abc
  if ("abc" %in% metric) {
    res$a <- abcp$a
    res$b <- abcp$b
    res$c <- abcp$c
  }
  
  if ("ABC" %in% metric) {
    res$A <- abca$A
    res$B <- abca$B
    res$C <- abca$C
  }

  # Create output class
  class(res) <- append("bioregion.pairwise.metric", class(res))
  
  # Inform nature of the output
  attr(res, "type") <- "similarity"
  
  # Return the output
  return(res)
}
