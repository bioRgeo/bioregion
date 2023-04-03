#' Search for an optimal number of clusters in a list of partitions 
#'
#' This function aims at optimizing one or several criteria on a set of 
#' ordered partitions. It is usually applied to find one (or several) optimal
#' number(s) of clusters on, for example, a hierarchical tree to cut, or a 
#' range of partitions obtained from k-means or PAM. Users are advised to be 
#' careful if applied in other cases (e.g., partitions which are not ordered in 
#' an increasing or decreasing sequence, or partitions which are not related
#' to each other).
#' 
#' @param partitions a `bioregion.partition.metrics` object (output from 
#' [partition_metrics()] or a `data.frame` with the first two 
#' columns named "K" (partition name) and "n_clusters" (number of clusters) and
#' the following columns containing evaluation metrics (numeric values)
#' 
#' @param metrics_to_use character string or vector of character strings
#' indicating upon which metric(s) in `partitions` the optimal number of
#' clusters should be calculated. Defaults to `"all"` which means all 
#' metrics available in `partitions` will be used
#' 
#' @param criterion character string indicating the criterion to be used to
#' identify optimal number(s) of clusters. Available methods currently include
#' `"elbow"`,
#' `"increasing_step"`, `"decreasing_step"`, `"cutoff"`, 
#' `"breakpoints"`, `"min"` or
#' `"max"`. Default is `"elbow"`. See details.
#' 
#' @param step_quantile if `"increasing_step"` or `"decreasing_step"`,
#' specify here the quantile
#' of differences between two consecutive k to be used as the cutoff to
#' identify the most important steps in `eval_metric`
#' 
#' @param step_levels if `"increasing_step"` or `"decreasing_step"`, specify
#' here the number of largest steps to keep as cutoffs.
#' 
#' @param step_round_above a `boolean` indicating if the optimal number of 
#' clusters should be picked above or below the identified steps. Indeed, each
#' step will correspond to a sudden increase or decrease between partition X &
#' partition X+1: should the optimal partition be X+1 
#' (`step_round_above = TRUE`) or X (`step_round_above = FALSE`? 
#' Defaults to `TRUE` 
#' 
#' @param metric_cutoffs if `criterion = "cutoff"`, specify here the cutoffs
#' of `eval_metric` at which the number of clusters should be extracted
#' 
#' @param n_breakpoints specify here the number of breakpoints to look for in
#' the curve. Defaults to 1 
#' 
#' @param plot a boolean indicating if a plot of the first `eval_metric`
#' should be drawn with the identified optimal numbers of cutoffs
#'
#' @details
#' \loadmathjax
#'
#' This function explores the relationship evaluation metric ~ number of
#' clusters, and a criterion is applied to search an optimal number of
#' clusters.
#'
#' **Please read the note section about the following criteria.**
#'
#' Foreword: 
#' 
#' Here we implemented a set of criteria commonly found in the literature or
#' recommended in the bioregionalisation literature. Nevertheless, we also
#' advocate to move
#' beyond the "Search one optimal number of clusters" paradigm, and consider
#' investigating "multiple optimal numbers of clusters". Indeed, using only one
#' optimal number of clusters may simplify the natural complexity of biological
#' datasets, and, for example, ignore the often hierarchical / nested nature of
#' bioregionalisations. Using multiple partitions likely avoids this
#' oversimplification bias and may convey more information.
#' See, for example, the reanalysis of Holt et al. (2013)
#' by \insertCite{Ficetola2017}{bioregion}, where they used deep, intermediate
#' and shallow cuts. 
#' 
#' Following this rationale, several of the criteria implemented here can/will
#' return multiple "optimal" numbers of clusters, depending on user choices.
#'
#' **Criteria to find optimal number(s) of clusters**
#' \itemize{
#' \item{`elbow`:
#' This method consists in finding one elbow in the evaluation metric curve, as
#' is commonly done in clustering analyses. The idea is to approximate the
#' number of clusters at which the evaluation metric no longer increments.It is
#' based on a fast method finding the maximum distance between the curve and a
#' straight line linking the minimum and maximum number of points. The code we
#' use here is based on code written by Esben Eickhardt available here
#' <https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve/42810075#42810075>.
#' The code has been modified to work on both increasing and decreasing
#' evaluation metrics.}
#' \item{`increasing_step` or `decreasing_step`:
#' This method consists in identifying clusters at the most important changes,
#' or steps, in the evaluation metric. The objective can be to either look for
#' largest increases (`increasing_step`) or largest decreases
#' `decreasing_step`. Steps are calculated based on the pairwise differences
#' between partitions. Therefore, this is relative to the distribution of
#' differences in the evaluation metric over the tested partitions. Specify
#' `step_quantile` as the quantile cutoff above which steps will be selected as
#' most important (by default, 0.99, i.e. the largest 1\% steps will be
#' selected).Alternatively, you can also choose to specify the number of top
#' steps to keep, e.g. to keep the largest three steps, specify
#' `step_level = 3`. Basically this method will emphasize the most important
#' changes in the evaluation metric as a first approximation of where important
#' cuts can be chosen.
#' 
#' **Please note that you should choose between `increasing_step` and
#' `decreasing_step` depending on the nature of your evaluation metrics. For
#' example, for metrics that are monotonously decreasing (e.g., endemism
#' metrics `"avg_endemism" & "tot_endemism"`) with the number of clusters
#' should n_clusters, you should choose `decreasing_step`. On the contrary, for
#' metrics that are monotonously increasing with the number of clusters (e.g.,
#' `"pc_distance"`), you should choose `increasing_step`. **
#' }
#' \item{`cutoffs`:
#' This method consists in specifying the cutoff value(s) in the evaluation
#' metric from which the number(s) of clusters should be derived. This is the
#' method used by \insertCite{Holt2013}{bioregion}. Note, however, that the
#' cut-offs suggested by Holt et al. (0.9, 0.95, 0.99, 0.999) may be only
#' relevant at very large spatial scales, and lower cut-offs should be
#' considered at finer spatial scales.
#' }
#' \item{`breakpoints`:
#' This method consists in finding break points in the curve using a segmented
#' regression. Users have to specify the number of expected break points in
#' `n_breakpoints` (defaults to 1). Note that since this method relies on a 
#' regression model, it should probably not be applied with a low number of
#' partitions.}
#'
#' \item{`min` & `max`:
#' Picks the optimal partition(s) respectively at the minimum or maximum value
#' of the evaluation metric.}
#' }
#' @return
#' a `list` of class `bioregion.optimal.n` with three elements:
#' \itemize{
#' \item{`args`: input arguments
#' }
#' \item{`evaluation_df`: the input evaluation data.frame appended with
#' `boolean` columns identifying the optimal numbers of clusters
#' }
#' \item{`optimal_nb_clusters`: a list containing the optimal number(s)
#' of cluster(s) for each metric specified in `"metrics_to_use"`, based on
#' the chosen `criterion`
#' }
#' \item{`plot`: if requested, the plot will be stored in this slot}}
#' @note Please note that finding the optimal number of clusters is a procedure
#' which normally requires decisions from the users, and as such can hardly be
#' fully automatized. Users are strongly advised to read the references
#' indicated below to look for guidance on how to choose their optimal
#' number(s) of clusters. Consider the "optimal" numbers of clusters returned
#' by this function as first approximation of the best numbers for your 
#' bioregionalisation.
#' 
#' @references
#' \insertRef{Castro-Insua2018}{bioregion}
#'
#' \insertRef{Ficetola2017}{bioregion}
#'
#' \insertRef{Holt2013}{bioregion}
#'
#' \insertRef{Kreft2010}{bioregion}
#'
#' \insertRef{Langfelder2008}{bioregion}
#' 
#' @importFrom rlang .data
#' 
#' @author
#' Boris Leroy (\email{leroy.boris@gmail.com}),
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) and
#' Pierre Denelle (\email{pierre.denelle@gmail.com})
#' 
#' @examples
#' \donttest{
#' dissim <- dissimilarity(fishmat, metric = "all")
#'
#' # User-defined number of clusters
#' tree1 <- hclu_hierarclust(dissim,
#'                           n_clust = 2:50,
#'                           index = "Simpson")
#' tree1
#' 
#' a <- partition_metrics(tree1,
#'                    dissimilarity = dissim,
#'                    net = fishdf,
#'                    species_col = "Species",
#'                    site_col = "Site",
#'                    eval_metric = c("tot_endemism",
#'                                    "avg_endemism",
#'                                    "pc_distance",
#'                                    "anosim"))
#'                                    
#' find_optimal_n(a)
#' find_optimal_n(a, criterion = "increasing_step")
#' find_optimal_n(a, criterion = "decreasing_step")
#' find_optimal_n(a, criterion = "decreasing_step",
#'                step_levels = 3) 
#' find_optimal_n(a, criterion = "decreasing_step",
#'                step_quantile = .9) 
#' find_optimal_n(a, criterion = "decreasing_step",
#'                step_levels = 3) 
#' find_optimal_n(a, criterion = "decreasing_step",
#'                step_levels = 3)                 
#' find_optimal_n(a, criterion = "breakpoints")             
#' }
#'
#' @importFrom stats predict quantile
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes_string geom_line facet_wrap geom_vline
#' @importFrom ggplot2 theme_bw
#' 
#' @export

find_optimal_n <- function(
    partitions,
    metrics_to_use = "all",
    criterion = "elbow", 
    step_quantile = .99,
    step_levels = NULL,
    step_round_above = TRUE,
    metric_cutoffs = c(.5, .75, .9, .95, .99, .999),
    n_breakpoints = 1,
    plot = TRUE){
  
  if(!inherits(partitions, "bioregion.partition.metrics")){
    if(!inherits(partitions, "data.frame")){
      stop("partitions should be the output object from partition_metrics()",
           "or a data.frame")
    } else {
      if(all(colnames(partitions)[1:2] == c("K", "n_clusters")) &
         ncol(partitions) > 2) {
        if(all(sapply(partitions[, 3:ncol(partitions)],
                      is.numeric))) {
          partitions <- list(
            args = list(eval_metric =
                          colnames(partitions)[3: ncol(partitions)]),
            evaluation_df = partitions)
        } else {
          stop("Your partition data.frame contains non numeric columns. ", 
               "Only numeric columns are expected after the first two columns")
        }
        
      } else {
        stop(
          "partitions should be the output object from partition_metrics() ",
          "or a data.frame with the first two columns named as 'K' & ",
          "'n_clusters' and the following columns being the evaluation ",
          "metrics")
      }
    }
  }
  
  if(length(criterion) > 1) {
    stop("Please choose only one criterion")
  }
  
  if("all" %in% metrics_to_use) {
    metrics_to_use <- partitions$args$eval_metric
  } else if(any(!(metrics_to_use %in% colnames(partitions$evaluation_df)))) {
    stop("metrics_to_use should exist in the evaluation table")
  }
  print(metrics_to_use)
  if(nrow(partitions$evaluation_df) <= 4){
    stop("The number of partitions is too low (<=4) for this function to work
         properly")
  } else{
    message(paste0("Number of partitions: ",
                   nrow(partitions$evaluation_df), "\n"))
    
    if(criterion %in% c("elbow", "increasing_step", "decreasing_step",
                        "breakpoints")) {
      if(nrow(partitions$evaluation_df) <= 8) {
        message(paste0(
          "...Caveat: be cautious with the interpretation of metric analyses
          with such a low number of partitions"))
      }
    } else if(!(criterion %in% c("min", "max", "cutoff"))) {
      stop("criterion must be one of elbow, increasing_step, decreasing_step,
           min, max, cutoff or breakpoints")
    }
    
    if(criterion %in% c("increasing_step", "decreasing_step")){
      if(!is.numeric(step_quantile) || step_quantile <= 0 ||
         step_quantile >= 1){
        stop(
          "step_quantile must be a numeric in the ]0,1[ interval. See help of
           the function.")
      }
    }
    
    if(criterion %in% c("increasing_step", "decreasing_step")){
      if(!is.logical(step_round_above)){
        stop("step_round_above must be a boolean. See help of the function.")
      }
    }
    
    if(criterion %in% c("increasing_step", "decreasing_step")){
      if(!is.null(step_levels) && (!is.numeric(step_levels) ||
                                   step_levels < 0)){
        stop("step_levels must be a positive integer. See help of the
             function.")
      }
    }
    
    if(!is.logical(plot)){
      stop("plot should be a Boolean.")
    }
    
    message(paste0(
      "Searching for potential optimal number(s) of clusters based on the ",
      criterion, " method"))
    
    partitions$evaluation_df <- data.frame(
      partitions$evaluation_df,
      array(FALSE,
            dim = c(nrow(partitions$evaluation_df),
                    length(metrics_to_use)),
            dimnames = list(NULL,
                            paste0("optimal_n_", metrics_to_use))))
    if(criterion == "breakpoints"){
      controls(args = n_breakpoints, type = "positive_integer" )
      
      seg_res <- lapply(
        metrics_to_use,
        function(x, eval_df, n_breaks) {
          if(any(is.na(eval_df[, x]))){
            # NA_vals <- eval_df[which(is.na(eval_df[, x])), ]
            eval_df <- eval_df[-which(is.na(eval_df[, x])), ]
          } # else {
          # NA_vals <- NULL
          # }
          
          n_cl <- eval_df$n_clusters
          y <- eval_df[, x]
          fit <- stats::lm(y ~ n_cl)
          fit_seg <- segmented::segmented(fit, 
                                          npsi = n_breaks)
          preds <- data.frame(n_clusters = n_cl,
                              preds = segmented::broken.line(fit_seg)$fit,
                              metric = x)
          
          optim_cutoffs <- fit_seg$psi[, 2]
          
          
          # if(length(NA_vals)){
          #   preds <- rbind(data.frame(NA_vals, 
          #                             preds = NA),
          #                  preds)
          # }
          
          return(list(optim_cutoffs,
                      preds))
        }, eval_df = partitions$evaluation_df, n_breaks = n_breakpoints)
      
      
      optim_n <- lapply(seg_res, function(x) x[[1]])
      names(optim_n) <- metrics_to_use
      
      #Rounding to get the closest partition
      optim_n <- lapply(optim_n, round)
      if(any(!(unlist(optim_n) %in% partitions$evaluation_df$n_clusters))) {
        message("Exact break point not in the list of partitions: finding the",
                " closest partition...\n")
        for(m in names(optim_n)) {
          if(any(!(optim_n[[m]] %in% partitions$evaluation_df$n_clusters))) {
            optim_n[[m]][which(!(optim_n[[m]] %in% 
                                   partitions$evaluation_df$n_clusters))] <- 
              partitions$evaluation_df$n_clusters[
                which.min(abs(
                  partitions$evaluation_df$n_clusters - 
                    optim_n[[m]][which(!(optim_n[[m]] %in% 
                                           partitions$evaluation_df$n_clusters))
                                 ]
                  ))
                ]
          }
        }
      }

      
      seg_preds <- data.frame(data.table::rbindlist(lapply(seg_res,
                                                           function(x) x[[2]])))
   
      for(metric in names(optim_n)) {
        partitions$evaluation_df[which(partitions$evaluation_df$n_clusters %in%
                                         optim_n[[metric]]),
                                 paste0("optimal_n_", metric)] <- TRUE
      }
    }
    
    if(criterion == "elbow"){
      optim_n <- lapply(metrics_to_use,
                        function(x, eval_df) {
                          .elbow_finder(eval_df$n_clusters,
                                        eval_df[, x],
                                        correct_decrease = TRUE)[1]
                        }, eval_df = partitions$evaluation_df)
      names(optim_n) <- metrics_to_use
      
      for(metric in names(optim_n)) {
        partitions$evaluation_df[which(partitions$evaluation_df$n_clusters == 
                                         optim_n[[metric]]), 
                                 paste0("optimal_n_", metric)] <- TRUE
      }
      
      message(paste0("   * elbow found at:"))
      message(paste(paste(names(optim_n), optim_n, sep = " "),
                    collapse = "\n"))
      
      if("anosim" %in% metrics_to_use){
        warning(paste0(
          "The elbow method is likely not suitable for the ANOSIM",
          " metric. You should rather look for leaps in the curve",
          " (see criterion = 'increasing_step' or ",
          "decreasing_step)"))
      }
    }
    
    if(criterion %in% c("increasing_step", "decreasing_step")) {
      message(" - Step method")
      
      if(criterion == "increasing_step" & any(c("tot_endemism",
                                                "avg_endemism") %in% 
                                              metrics_to_use)) {
        warning(paste0(
          "Criterion 'increasing_step' cannot work properly with ",
          "metric 'tot_endemism', because this metric is usually ",
          "monotonously decreasing. Consider using ",
          "criterion = 'decreasing_step' instead.\n"))
      } else if(criterion == "decreasing_step" & any(c("pc_distance") %in% 
                                                     metrics_to_use)) {
        warning(paste0(
          "Criterion 'decreasing_step' cannot work properly with",
          " metrics 'pc_distance' or 'avg_endemism', because these",
          " metrics are usually monotonously decreasing. Consider ",
          "using criterion = 'increasing_step' instead.\n"))
      }
      
      optim_index <- lapply(
        metrics_to_use,
        function (x, eval_df, crit, s_lvl, s_qt, cl) {
          # Compute difference between each consecutive nb of clusters
          diffs <- eval_df[2:nrow(eval_df), x] -
            eval_df[1:(nrow(eval_df) - 1), x]
          if(crit == "decreasing_step"){
            diffs <- - diffs
          }
          
          if(!is.null(s_lvl)){
            level_diffs <- diffs[order(diffs, decreasing = TRUE)][1:s_lvl]
            optim_n <- which(diffs %in% level_diffs) + cl
            
            if(length(optim_n) >= s_lvl + 2) {
              warning(paste0("The number of optimal N for method '",
                             x, "' is suspiciously high, consider ",
                             "switching between 'increasing_step'",
                             " and 'decreasing_step'\n"))
            }
          } else if(!is.null(step_quantile)){
            qt <- stats::quantile(diffs,
                                  s_qt)
            optim_n <- which(diffs > qt) + cl
          }
          return(optim_n)
        }, eval_df = partitions$evaluation_df, crit = criterion,
        s_lvl = step_levels, s_qt = step_quantile, cl = step_round_above)
      
      names(optim_index) <- metrics_to_use
      
      optim_n <- lapply(
        optim_index,
        function(x) unique(partitions$evaluation_df$n_clusters[x]))
      
      
      for(metric in names(optim_n)) {
        partitions$evaluation_df[optim_index[[metric]],
                                 paste0("optimal_n_", metric)] <- TRUE
      }
      
    }
    if(criterion == "cutoff") {
      message(" - Cutoff method")
      
      if(length(metrics_to_use) > 1) {
        stop("Criterion 'cutoff' should probably be used with only one ",
             "evaluation metric (you have ",
             length(metrics_to_use),
             " evaluation metrics in 'partitions'). Indeed, metrics have ",
             "distinct orders of magnitude, and so the 'metric_cutoffs' you ",
             " chose are likely to be",
             " appropriate for only one of the metrics, but no the others.")
      }
      
      optim_index <- sapply(metric_cutoffs,
                            function(cutoff, vals) which(vals >= cutoff)[1],
                            vals = partitions$evaluation_df[, metrics_to_use])
      
      partitions$evaluation_df[, paste0("optimal_n_", metrics_to_use)] <- FALSE
      partitions$evaluation_df[optim_index, paste0("optimal_n_", 
                                                   metrics_to_use)] <- TRUE
      
      optim_n <- list(partitions$evaluation_df$n_clusters[optim_index])
      names(optim_n) <- metrics_to_use
    }
    
    if(criterion == "max"){
      message(" - Max value method")
      
      optim_index <- lapply(metrics_to_use,
                            function(x, eval_df) {
                              which(eval_df[, x] == max(eval_df[, x]))
                            }, eval_df = partitions$evaluation_df)
      names(optim_index) <- metrics_to_use
      
      optim_n <- lapply(
        optim_index,
        function(x) unique(partitions$evaluation_df$n_clusters[x]))
      
      for(metric in names(optim_n)) {
        partitions$evaluation_df[optim_index,
                                 paste0("optimal_n_", metric)] <- TRUE
      }
    }
    
    if(criterion == "min"){
      message(" - Min value method")
      
      optim_index <- lapply(metrics_to_use,
                            function(x, eval_df) {
                              which(eval_df[, x] == min(eval_df[, x]))
                            }, eval_df = partitions$evaluation_df)
      names(optim_index) <- metrics_to_use
      
      optim_n <- lapply(
        optim_index,
        function(x) unique(partitions$evaluation_df$n_clusters[x]))
      
      for(metric in names(optim_n)) {
        partitions$evaluation_df[optim_index, paste0("optimal_n_", metric)] <- 
          TRUE
      }
    }
    
    if(plot){
      ggdf2 <- tidyr::pivot_longer(
        data = as.data.frame(partitions$evaluation_df),
        cols = grep("optimal_n_", colnames(partitions$evaluation_df)),
        names_to = "variable")
      ggdf2 <- as.data.frame(ggdf2)
      ggdf2$variable <- gsub("optimal_n_", "", ggdf2$variable)
      ggdf2 <- ggdf2[ggdf2$value, ]
      ggdf <- as.data.frame(tidyr::pivot_longer(
        data = partitions$evaluation_df, cols = metrics_to_use,
        names_to = "variable"))
      
      message("Plotting results...")

      if(criterion == "breakpoints"){
        
        names(seg_preds) <- c("n_clusters", "value", "variable")
        
        message("   (the red line is the prediction from the segmented ",
                "regression)")
        
        p <- ggplot2::ggplot(ggdf, ggplot2::aes_string(x = "n_clusters",
                                                       y = "value")) +
          ggplot2::geom_line(col = "darkgrey") +
          ggplot2::facet_wrap(~ variable, scales = "free_y") +
          # ggplot2::geom_hline(yintercept = partitions$evaluation_df[
          #   partitions$evaluation_df$optimal_nclust, eval_metric[1]],
          #                     linetype = 2) +
          ggplot2::geom_vline(data = ggdf2,
                              ggplot2::aes_string(xintercept = "n_clusters"),
                              linetype = 2) +
          ggplot2::theme_bw() +
          ggplot2::geom_line(data = seg_preds,
                             col = "red")
      } else {
        p <- ggplot2::ggplot(ggdf, ggplot2::aes_string(x = "n_clusters",
                                                       y = "value")) +
          ggplot2::geom_line(col = "darkgrey") +
          ggplot2::facet_wrap(~ variable, scales = "free_y") +
          # ggplot2::geom_hline(yintercept = partitions$evaluation_df[
          #   partitions$evaluation_df$optimal_nclust, eval_metric[1]],
          #                     linetype = 2) +
          ggplot2::geom_vline(data = ggdf2,
                              ggplot2::aes_string(xintercept = "n_clusters"),
                              linetype = 2) +
          ggplot2::theme_bw()
      }
      print(p)
    } else {
      p <- FALSE
    }
  }
  
  outputs <- list(args = list(metrics_to_use = metrics_to_use,
                              criterion = criterion, 
                              step_quantile = step_quantile,
                              step_levels = step_levels,
                              step_round_above = step_round_above,
                              metric_cutoffs = metric_cutoffs,
                              n_breakpoints = n_breakpoints,
                              plot = plot
  ),
  evaluation_df = partitions$evaluation_df,
  optimal_nb_clusters = optim_n,
  plot = p)
  
  class(outputs) <- append("bioregion.optimal.n", class(outputs))
  return(outputs)
}
