#' Compute similarity metrics between sites based on species composition
#'
#' This function generates a `data.frame` where each row provides one or
#' several similarity metrics between pairs of sites, based on a co-occurrence
#' `matrix` with sites as rows and species as columns.
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
#' `bioregion.pairwise`, containing one or several similarity
#' metrics between pairs of sites. The first two columns represent the pairs of 
#' sites. There is one column per similarity metric provided in `metric` and
#' `formula`, except for the `abc` and `ABC` metrics, which are stored in three 
#' separate columns (one for each letter).
#' 
#' @details
#' With `a` the number of species shared by a pair of sites, `b`
#' species only present in the first site and `c` species only present in
#' the second site.
#'
#' Jaccard = 1 - (b + c) / (a + b + c)
#'
#' Jaccardturn = 1 - 2min(b, c) / (a + 2min(b, c)) (Baselga, 2012)
#'
#' Sorensen = 1 - (b + c) / (2a + b + c)
#'
#' Simpson = 1 - min(b, c) / (a + min(b, c))
#'
#' If abundances data are available, Bray-Curtis and its turnover component can
#' also be computed with the following equation:
#'
#' Bray = 1 - (B + C) / (2A + B + C)
#'
#' Brayturn = 1 - min(B, C) / (A + min(B, C)) (Baselga, 2013)
#'
#' with `A` the sum of the lesser values for common species shared by a pair of
#' sites. `B` and `C` are the total number of specimens counted at both sites 
#' minus `A`.
#'
#' `formula` can be used to compute customized metrics with the terms
#' `a`, `b`, `c`, `A`, `B`, and `C`. For example
#' `formula = c("1 - pmin(b,c) / (a + pmin(b,c))", "1 - (B + C) / (2*A + B + C)")`
#' will compute the Simpson and Bray-Curtis similarity metrics, respectively. 
#' Note that `pmin` is used in the Simpson formula because `a`, `b`, `c`, `A`, 
#' `B` and `C` are `numeric` vectors.
#'
#' Euclidean computes the Euclidean similarity between each pair of sites
#' following this equation:
#'
#' Euclidean = 1 / (1 + d_ij)
#'
#' Where d_ij is the Euclidean distance between site i and 
#' site j in terms of species composition.
#' 
#' @references
#' Baselga A (2012) The Relationship between Species Replacement,
#' Dissimilarity Derived from Nestedness, and Nestedness.
#' \emph{Global Ecology and Biogeography} 21, 1223--1232.
#' 
#' Baselga A (2013) Separating the two components of abundance-based
#' dissimilarity: balanced changes in abundance vs. abundance gradients.
#' \emph{Methods in Ecology and Evolution} 4, 552--557.
#' 
#' @seealso 
#' For more details illustrated with a practical example, 
#' see the vignette: 
#' \url{https://biorgeo.github.io/bioregion/articles/a3_pairwise_metrics.html}.
#' 
#' Associated functions: 
#' [dissimilarity] [similarity_to_dissimilarity]
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
#' sim <- similarity(comat, metric = c("abc", "ABC", "Simpson", "Brayturn"))
#'
#' sim <- similarity(comat, metric = "all",
#' formula = "1 - (b + c) / (a + b + c)")
#' 
#' @export
similarity <- function(comat, 
                       metric = "Simpson", 
                       formula = NULL,
                       method = "prodmat"){
  
  # List of metrics based on abc
  lsmetricabc <- c("abc", "Jaccard", "Jaccardturn", "Sorensen", "Simpson")
  
  # List of metrics based on ABC
  lsmetricABC <- c("ABC", "Bray", "Brayturn")
  
  # List of metrics based on other features
  lsmetrico <- c("Euclidean")
  
  # Control inputs
  controls(args = method, data = NULL, type = "character")
  if (!(method %in% c("prodmat", "loops"))) {
    stop(paste0("Please choose method from the following:\n",
                "prodmat or loops."), 
         call. = FALSE)
  }
  
  if (is.null(metric) & is.null(formula)) {
    stop("metric or formula should be used.", call. = FALSE)
  }
  
  if(!is.null(metric)){
    controls(args = metric, data = NULL, type = "character_vector")
    if ("all" %in% metric) {
      metric <- c(lsmetricabc, lsmetricABC, lsmetrico)
    }
    if (length(intersect(c(lsmetricabc, lsmetricABC, lsmetrico), metric)) !=
        length(metric)) {
      stop(paste0("One or several metric(s) chosen are not", 
                  " available.\n",
                  "Please choose from the following:\n",
                  "abc, Jaccard, Jaccardturn, Sorensen, Simpson, ABC, ",
                  "Bray, Brayturn and Euclidean."),
           call. = FALSE)
    }
  }

  if(!is.null(formula)){
    controls(args = formula, data = NULL, type = "character_vector")
  }
  abcinformula <- (sum(c(grepl("a", formula), grepl("b", formula),
                         grepl("c", formula))) > 0)
  ABCinformula <- (sum(c(grepl("A", formula), grepl("B", formula),
                         grepl("C", formula))) > 0)  
  
  controls(args = NULL, data = comat, type = "input_matrix")
  minco <- min(comat)
  if (minco < 0) {
    if("Euclidean" %in% metric){
      message("Negative value(s) detected in comat!")
    }else{
      stop("Negative value detected in comat.", 
           call. = FALSE)
    }
  }
  
  # Extract site id
  siteid <- rownames(comat)
  if(is.null(siteid)){
    siteid <- as.character(1:dim(comat)[1])
    message("No rownames detected, they have been assigned automatically.")
  }
  
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
  
  # Compute equation in formula
  if (!is.null(formula)) {
    for (k in 1:length(formula)) {
      res <- cbind(res, eval(parse(text = formula[k])))
      colnames(res)[dim(res)[2]] <- formula[k]
    }
  }

  # Create output class
  class(res) <- append("bioregion.pairwise", class(res))
  
  # Inform nature of the output
  attr(res, "type") <- "similarity"
  attr(res, "nb_sites") <- nrow(comat)
  attr(res, "nb_species") <- ncol(comat)
  
  # Return the output
  return(res)
}
