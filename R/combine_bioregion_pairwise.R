#' Combine and enrich bioregion (dis)similarity object(s)
#'
#' Combine two `bioregion.pairwise.metric` objects and/or compute new pairwise
#' metrics based on the columns of the object(s).
#'
#' @param primary_metrics A `bioregion.pairwise.metric` object. This is the 
#' main set of pairwise metrics that will be used as a base for the combination.
#'   
#' @param secondary_metrics A second `bioregion.pairwise.metric` 
#' object to be combined with `primary_metrics`. Must have the same sites 
#' identifiers and pairwise structure. Can be set to `NULL` if `new_metrics` is 
#' specified.
#'   
#' @param new_metrics  A `character` vector or a single `character` string 
#' specifying custom formula(s) based on the column names of `primary_metrics` 
#' and `secondary_metrics` (see Details). The default is `NULL`.
#' 
#' @details
#' When both `primary_metrics` and `secondary_metrics` are provided and if the
#' pairwise structure is identical the function combine the two objects. If 
#' `new_metrics` is provided, each formula is evaluated based on the column 
#' names of `primary_metrics` (and `secondary_metrics` if provided). 
#' 
#' @return 
#' A new `bioregion.pairwise.metric` object containing the combined and/or
#' enriched data. It includes all original metrics from the inputs, as well as 
#' any newly computed metrics.
#' 
#' @seealso 
#' For more details illustrated with a practical example, 
#' see the vignette: 
#' \url{https://biorgeo.github.io/bioregion/articles/a3_pairwise_metrics.html}.
#' 
#' Associated functions: 
#' [dissimilarity] [similarity] [as_bioregion_pairwise]
#' 
#' @author
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) \cr
#' Boris Leroy (\email{leroy.boris@gmail.com}) \cr
#' Pierre Denelle (\email{pierre.denelle@gmail.com})
#' 
#' @examples
#' comat <- matrix(sample(0:1000, size = 50, replace = TRUE,
#' prob = 1 / 1:1001), 5, 10)
#' rownames(comat) <- paste0("s", 1:5)
#' colnames(comat) <- paste0("sp", 1:10)
#'
#' sim <- combine_bioregion_pairwise(primary_metrics = similarity(comat, 
#'                                                                metric = "abc"),
#'                                   secondary_metrics = similarity(comat, 
#'                                                                  metric = "Simpson"),
#'                                   new_metrics = "1 - (b + c) / (a + b + c)")
#' 
#' @export
combine_bioregion_pairwise <- function(primary_metrics,
                                       secondary_metrics,
                                       new_metrics = NULL) {
  
  # Control primary_metrics
  controls(args = NULL, data = primary_metrics, type = "input_pairwise")
  type1 <- attr(primary_metrics, "type")
  nbs1 <- attr(primary_metrics, "nb_sites")
  nbsp1 <- attr(primary_metrics, "nb_species")
  if(type1 != "similarity" & type1 != "dissimilarity"){
    stop(paste0("primary_metrics",
                " is a bioregion.pairwise.metric object but it has not ",
                "been possible to identify the object's type (similarity or ",
                " dissimilarity) probably because the ",
                "bioregion.pairwise.metric object has been altered."),
         call. = FALSE)
  }
  
  # Control secondary_metrics and combine
  if(!is.null(secondary_metrics)){
    controls(args = NULL, data = secondary_metrics, type = "input_pairwise")
    type2 <- attr(secondary_metrics, "type") 
    nbs2 <- attr(secondary_metrics, "nb_sites")
    nbsp2 <- attr(secondary_metrics, "nb_species")
    if(type2 != "similarity" & type2 != "dissimilarity"){
      stop(paste0("secondary_metrics",
                  " is a bioregion.pairwise.metric object but it has not ",
                  "been possible to identify the object's type (similarity or ",
                  " dissimilarity) probably because the ",
                  "bioregion.pairwise.metric object has been altered."),
           call. = FALSE)
    }
    if(type1 != type2){
      stop(paste0("primary_metrics and secondary_metrics should have the same ",
                  "type (similarity or dissimilarity)."),
           call. = FALSE)
    }
    if(nbs1 != nbs2){
      stop(paste0("primary_metrics and secondary_metrics should have the same ",
                  "number of sites."),
      call. = FALSE)
    }
    if(is.na(nbsp1) | is.na(nbsp2) | nbsp1 != nbsp2){
      message(paste0("primary_metrics and secondary_metrics are based ",
                     "on different number of species."))
      nbsp <- NA
    }else{
      nbsp <- nbsp1
    }
    site1 <- paste0(primary_metrics[,1],
                    "_",
                    primary_metrics[,2])
    site2 <- paste0(secondary_metrics[,1],
                    "_",
                    secondary_metrics[,2])
    if(sum(site1 == site2) != length(site1)){
      stop(paste0("primary_metrics and secondary_metrics should have the same ",
                  "sites identifiers and pairwise structure."),
      call. = FALSE)
    }else{
      mat <- cbind(primary_metrics, secondary_metrics[, -c(1,2)])
      if(dim(secondary_metrics)[2]==3){
        colnames(mat)[dim(mat)[2]] <- colnames(secondary_metrics)[3]
      }
    }
  }else{
    mat <- primary_metrics
    nbsp <- nbsp1
  }

  # Control new_metrics and compute
  if(!is.null(new_metrics)){
    controls(args = new_metrics, data = NULL, type = "character_vector")
    
    # Compute equation in new_metrics
    for (k in 1:length(new_metrics)) {
      mat <- cbind(mat, with(mat, eval(parse(text = new_metrics[k]))))
      colnames(mat)[dim(mat)[2]] <- new_metrics[k]
    }


  }

  # Export output
  attr(mat, "type") <- type1
  attr(mat, "nb_sites") <- nbs1
  attr(mat, "nb_species") <- nbsp
  
  class(mat) <- append("bioregion.pairwise.metric", class(mat))
  
  return(mat)
  
}
