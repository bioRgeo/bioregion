#' Search for an optimal number of bioregions in a bioregionalization 
#'
#' This function aims to optimize one or several criteria on a set of 
#' ordered partitions from a bioregionalization. It is typically used to find 
#' one or more optimal bioregion counts on hierarchical trees to cut or ranges 
#' of partitions from k-means or PAM. Users should exercise caution in 
#' other cases (e.g., unordered partitions or unrelated 
#' partitions).
#' 
#' @param evaluation_df a `data.frame` output from 
#' [bioregionalization_metrics()]) or a `data.frame` with the first two 
#' columns with partition name and number of bioregions, 
#' followed by columns with numeric evaluation metrics.
#' 
#' @param metrics_to_use A `character` vector or single string specifying 
#' metrics in `evaluation_df` for calculating optimal bioregions. Defaults 
#' to `"all"` (uses all metrics).
#' 
#' @param criterion A `character` string specifying the criterion to identify 
#' optimal clusters. Options include `"elbow"`, `"increasing_step"`, 
#' `"decreasing_step"`, `"cutoff"`, `"breakpoints"`, `"min"`, or `"max"`. 
#' Defaults to `"elbow"`. See Details.
#' 
#' @param step_quantile For `"increasing_step"` or `"decreasing_step"`, 
#' specifies the quantile of differences between consecutive bioregionalizations as 
#' the cutoff to identify significant steps in `eval_metric`.
#' 
#' @param step_levels For `"increasing_step"` or `"decreasing_step"`, specifies 
#' the number of largest steps to retain as cutoffs.
#' 
#' @param step_round_above A `boolean` indicating whether the optimal clusters 
#' are above (`TRUE`) or below (`FALSE`) identified steps. Defaults to `TRUE`.
#' 
#' @param metric_cutoffs For `criterion = "cutoff"`, specifies the cutoffs 
#' of `eval_metric` to extract cluster counts.
#' 
#' @param n_breakpoints Specifies the number of breakpoints to find in the 
#' curve. Defaults to 1.
#' 
#' @param plot A `boolean` indicating if a plot of the first `eval_metric` 
#' with identified optimal clusters should be drawn.
#' 
#' @param verbose A `boolean` indicating whether to 
#' display progress messages. Set to `FALSE` to suppress these messages.
#' 
#' @param bioregionalizations Deprecated.
#'
#' @return
#' A `list`containing these elements:
#' \itemize{
#' \item{`args`: Input arguments.}
#' \item{`evaluation_df`: The input evaluation `data.frame`, appended with 
#' `boolean` columns for optimal cluster counts.}
#' \item{`optimal_n`: A `list` with optimal cluster counts for each 
#' metric in `"metrics_to_use"`, based on the chosen `criterion`.}
#' \item{`plot`: The plot (if requested).}}
#'
#' @details
#' This function explores evaluation metric ~ cluster relationships, applying 
#' criteria to find optimal cluster counts.
#'
#' **Note on criteria:** Several criteria can return multiple optimal cluster 
#' counts, emphasizing hierarchical or nested bioregionalizations. This 
#' approach aligns with modern recommendations for biological datasets, as seen 
#' in Ficetola et al. (2017)'s reanalysis of Holt et al. (2013).
#' 
#' **Criteria for optimal clusters:** 
#' \itemize{
#' \item{`elbow`: Identifies the "elbow" point in the evaluation metric curve, 
#' where incremental improvements diminish. Based on a method to find the 
#' maximum distance from a straight line linking curve endpoints.}
#' \item{`increasing_step` or `decreasing_step`: Highlights significant 
#' increases or decreases in metrics by analyzing pairwise differences between 
#' bioregionalizations. Users specify `step_quantile` or `step_levels`.}
#' \item{`cutoffs`: Derives clusters from specified metric cutoffs, e.g., as in 
#' Holt et al. (2013). Adjust cutoffs based on spatial scale.}
#' \item{`breakpoints`: Uses segmented regression to find breakpoints. Requires 
#' specifying `n_breakpoints`.}
#' \item{`min` & `max`: Selects clusters at minimum or maximum metric values.}}
#' 
#' @note 
#' Please note that finding the optimal number of clusters is a procedure
#' which normally requires decisions from the users, and as such can hardly be
#' fully automatized. Users are strongly advised to read the references
#' indicated below to look for guidance on how to choose their optimal
#' number(s) of clusters. Consider the "optimal" numbers of clusters returned
#' by this function as first approximation of the best numbers for your 
#' bioregionalization.
#' 
#' @references
#' Holt BG, Lessard J, Borregaard MK, Fritz SA, Araújo MB, Dimitrov D, Fabre P, 
#' Graham CH, Graves GR, Jønsson Ka, Nogués-Bravo D, Wang Z, Whittaker RJ, 
#' Fjeldså J & Rahbek C (2013) An update of Wallace's zoogeographic regions of 
#' the world. \emph{Science} 339, 74-78.
#' 
#' Ficetola GF, Mazel F & Thuiller W (2017) Global determinants of 
#' zoogeographical boundaries. \emph{Nature Ecology & Evolution} 1, 0089.
#' 
#' @seealso 
#' For more details illustrated with a practical example, 
#' see the vignette: 
#' \url{https://biorgeo.github.io/bioregion/articles/a4_1_hierarchical_clustering.html#optimaln}.
#' 
#' Associated functions: 
#' [hclu_hierarclust]
#' 
#' @author
#' Boris Leroy (\email{leroy.boris@gmail.com}) \cr
#' Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) \cr
#' Pierre Denelle (\email{pierre.denelle@gmail.com})
#' 
#' @examples
#' comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
#' 20, 25)
#' rownames(comat) <- paste0("Site",1:20)
#' colnames(comat) <- paste0("Species",1:25)
#' 
#' dissim <- dissimilarity(comat, metric = "all")
#'
#' # User-defined number of clusters
#' bioreg <- hclu_hierarclust(dissim,
#'                            optimal_tree_method = "best",
#'                            n_clust = 5:10,
#'                            verbose = FALSE)
#' bioreg
#' 
#' evalmet <- bioregionalization_metrics(bioreg,
#'                                       eval_metrics = "anosim",
#'                                       dissimilarity = dissim)
#'                                    
#' find_optimal_n(evalmet, criterion = 'increasing_step', plot = FALSE)
#' 
#' @importFrom rlang .data
#' 
#' @export
find_optimal_n <- function(evaluation_df,
                           metrics_to_use = "all",
                           criterion = "elbow", 
                           step_quantile = .99,
                           step_levels = NULL,
                           step_round_above = TRUE,
                           metric_cutoffs = c(.5, .75, .9, .95, .99, .999),
                           n_breakpoints = 1,
                           plot = TRUE,
                           verbose = TRUE,
                           bioregionalizations = NULL){
  
  # Control deprecated
  if (!is.null(bioregionalizations)) {
    warning("bioregionalizations is deprecated.", 
            call. = FALSE)
  }
  
  # Control evaluation_df
  if(!inherits(evaluation_df, "bioregion.bioregionalization.metrics")){
    if(!inherits(evaluation_df, "data.frame")){
      stop(paste0("evaluation_df should be the output object from ", 
                  "bioregionalization_metrics() ",
                  "or a data.frame"), 
           call. = FALSE) 
    } else {
      if(ncol(evaluation_df) > 2) {
        if(all(sapply(evaluation_df[, 3:ncol(evaluation_df)],
                      is.numeric))) {
          #evaluation_df <- list(
          #  args = list(eval_metric = colnames(evaluation_df)[3: ncol(evaluation_df)]),
          #              evaluation_df = evaluation_df)
        } else {
          stop(paste0("evaluation_df contains non numeric ",
                      "columns. Only numeric columns are expected after the ",
                      "first two columns"), 
               call. = FALSE) 
        }
        
      } else {
        stop(paste0("evaluation_df should be the output object from ",
                    "bioregionalization_metrics() ",
                    "or a data.frame with the first two columns ",
                    "with partition name and number of bioregions ",
                    "followed by columns with numeric evaluation metrics"), 
             call. = FALSE) 
      }
    }
  }
  
  # Update colnmaes
  colnames(evaluation_df)[1:2] <- c("partition", "n_bioregions")
  
  # Control metrics_to_use
  controls(args = metrics_to_use, data = NULL, type = "character_vector")
  if("all" %in% metrics_to_use) {
    metrics_to_use <- colnames(evaluation_df)[-c(1,2)]
  } else if(any(!(metrics_to_use %in% colnames(evaluation_df)))) {
    stop(paste0("metrics_to_use should exist in evaluation_df."),
         call. = FALSE)
  }
  
  # Control criterion
  controls(args = criterion, data = NULL, type = "character")
  if (!(criterion %in% c("elbow", "increasing_step", "decreasing_step", 
                                "cutoff", "breakpoints", "min", "max"))) {
    stop(paste0("Please choose criterion from the following:\n",
                "elbow, increasing_step, decreasing_step, cutoff, breakpoints,",
                " min or max."), 
         call. = FALSE)   
  }
  
  # Control plot and verbose
  controls(args = plot, type = "boolean")
  controls(args = verbose, data = NULL, type = "boolean")
  
  # Verifying that metrics vary, otherwise we remove them
  nvals_metrics <- sapply(lapply(evaluation_df[, metrics_to_use,
                                                          drop = FALSE],
                               unique), length)
  if (any(nvals_metrics == 1)) {
    exclude_metrics <- names(nvals_metrics[which(nvals_metrics == 1)])
    warning(paste0("Metrics ", 
                   paste0(exclude_metrics, collapse = ", "),
                   " did not vary in evaluation_df, so they were removed."))
    metrics_to_use <- metrics_to_use[-which(metrics_to_use %in% 
                                              exclude_metrics)]
  } 
  if (any(nvals_metrics == 2) &
             criterion == "breakpoints") {
    exclude_metrics <- names(nvals_metrics[which(nvals_metrics == 2)])
    warning(paste0("Metrics ", 
                   paste0(exclude_metrics, collapse = ", "),
                   " do not vary sufficiently to use breakpoints",
                   ", so they were removed."))
    metrics_to_use <- metrics_to_use[-which(metrics_to_use %in% 
                                              exclude_metrics)]
  }
  if(!(length(metrics_to_use))) {
    stop(paste0("The selected evaluation metrics did not vary sufficiently ",
                "in input. Please check your metrics or increase ",
                "your range of partitions when computing ",
                "bioregionalization_metrics()"), 
         call. = FALSE)
  }
  
  #print(metrics_to_use)
  
  if(nrow(evaluation_df) <= 4){
    stop(paste0("The number of partitions is too low (<=4) ",
                "for this function to work properly"),
         call. = FALSE) 
         
  } else{
    if(verbose){
      message(paste0("Number of partitions: ",
                     nrow(evaluation_df), "\n"))
    }
    
    if(criterion %in% c("elbow", 
                        "increasing_step", 
                        "decreasing_step",
                        "breakpoints")) {
      if(nrow(evaluation_df) <= 8 & verbose) {
        message(paste0("...Caveat: be cautious with the interpretation of ",
                       "metric analyses with such a low number of partitions"))
      }
    } 
    #else if(!(criterion %in% c("min", "max", "cutoff"))) {
    #  stop("criterion must be one of elbow, increasing_step, decreasing_step,
    #       min, max, cutoff or breakpoints")
    #}
    
    if(criterion %in% c("increasing_step", "decreasing_step")){
      controls(args=step_quantile, type="strict_positive_numeric")
      if(step_quantile >= 1){
        stop(paste0("step_quantile must be in the ]0,1[ interval."), 
                    call. = FALSE)
      }
    }
    
    if(criterion %in% c("increasing_step", "decreasing_step")){
      controls(args=step_round_above, type="boolean")
    }
    
    if(criterion %in% c("increasing_step", "decreasing_step")){
      if(!is.null(step_levels)){
        controls(args = step_levels, type = "positive_integer") 
      }
    }
    
    if(verbose){
      message(paste0("Searching for potential optimal number(s) of clusters ",
                     "based on the ",
                     criterion, 
                     " method"))
    }
    
    evaluation_df <- data.frame(
      evaluation_df,
      array(FALSE,
            dim = c(nrow(evaluation_df),
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
          
          n_cl <- eval_df$n_bioregions
          y <- eval_df[, x]
          fit <- stats::lm(y ~ n_cl)
          if(fit$coefficients[2] != 0) {
            fit_seg <- segmented::segmented(fit, 
                                            npsi = n_breaks)
            preds <- data.frame(n_clusters = n_cl,
                                preds = segmented::broken.line(fit_seg)$fit,
                                metric = x)
            
            optim_cutoffs <- fit_seg$psi[, 2]
          } else {
            preds <- data.frame(n_clusters = n_cl,
                                preds = NA,
                                metric = x)
            optim_cutoffs <- NA
          }

          
          
          # if(length(NA_vals)){
          #   preds <- rbind(data.frame(NA_vals, 
          #                             preds = NA),
          #                  preds)
          # }
          
          return(list(optim_cutoffs,
                      preds))
        }, eval_df = evaluation_df, n_breaks = n_breakpoints)
      
      
      optim_n <- lapply(seg_res, function(x) x[[1]])
      names(optim_n) <- metrics_to_use
      
      #Rounding to get the closest bioregionalization
      optim_n <- lapply(optim_n, round)
      if(any(!(stats::na.omit(unlist(optim_n)) %in% evaluation_df$n_bioregions))) {
        if(verbose){
          message("Exact break point not in the list of partitions: finding the",
                  " closest partition...\n")
        }
        for(m in names(optim_n)) {
          if(any(!(optim_n[[m]] %in% evaluation_df$n_bioregions))) {
            for(cutoff in  optim_n[[m]][
              which(!(optim_n[[m]] %in% 
                      evaluation_df$n_bioregions))]) {
              optim_n[[m]][which(optim_n[[m]] == cutoff)] <- 
                evaluation_df$n_bioregions[
                  which.min(abs(
                    evaluation_df$n_bioregions - 
                      cutoff
                  ))
                ]
            }
          }
        }
      }

      
      seg_preds <- data.frame(data.table::rbindlist(lapply(seg_res,
                                                           function(x) x[[2]])))
   
      for(metric in names(optim_n)) {
        evaluation_df[which(evaluation_df$n_bioregions %in%
                                         optim_n[[metric]]),
                                 paste0("optimal_n_", metric)] <- TRUE
      }
    }
    
    if(criterion == "elbow"){
      optim_n <- lapply(metrics_to_use,
                        function(x, eval_df) {
                          elbow_finder(eval_df$n_bioregions,
                                       eval_df[, x],
                                       correct_decrease = TRUE)[1]
                        }, eval_df = evaluation_df)
      names(optim_n) <- metrics_to_use
      
      for(metric in names(optim_n)) {
        evaluation_df[which(evaluation_df$n_bioregions == 
                                         optim_n[[metric]]), 
                                 paste0("optimal_n_", metric)] <- TRUE
      }
      
      if(verbose){
        message(paste0("   * elbow found at:"))
        message(paste(paste(names(optim_n), optim_n, sep = " "),
                      collapse = "\n"))
      }
      
      if("anosim" %in% metrics_to_use & verbose){
        warning(paste0(
          "The elbow method is likely not suitable for the anosim",
          " metric. You should rather look for leaps in the curve",
          " (see criterion = 'increasing_step' or ",
          "decreasing_step)"))
      }
    }
    
    if(criterion %in% c("increasing_step", "decreasing_step")) {
      
      if(verbose){
        message(" - Step method")
      }  

      if(criterion == "increasing_step" & 
         verbose &
         any(c("tot_endemics", "mean_endemics") %in% metrics_to_use)) {
        warning(paste0(
          "Criterion 'increasing_step' cannot work properly with ",
          "metric 'tot_endemics', because this metric is usually ",
          "monotonously decreasing. Consider using ",
          "criterion = 'decreasing_step' instead.\n"))
      } else if(criterion == "decreasing_step" & 
                verbose &
                any(c("prop_between_dissim", "mean_endemics") %in% metrics_to_use)) {
        warning(paste0(
          "Criterion 'decreasing_step' cannot work properly with",
          " metrics 'prop_between_dissim' or 'mean_endemics', because these",
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
            
            if(length(optim_n) >= s_lvl + 2 & verbose) {
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
        }, eval_df = evaluation_df, crit = criterion,
        s_lvl = step_levels, s_qt = step_quantile, cl = step_round_above)
      
      names(optim_index) <- metrics_to_use
      
      optim_n <- lapply(
        optim_index,
        function(x) unique(evaluation_df$n_bioregions[x]))
      
      
      for(metric in names(optim_n)) {
        evaluation_df[optim_index[[metric]],
                                 paste0("optimal_n_", metric)] <- TRUE
      }
      
    }
    if(criterion == "cutoff") {
      
      if(verbose){
        message(" - Cutoff method")
      }
      
      if(length(metrics_to_use) > 1) {
        stop(paste0("Criterion 'cutoff' should probably be used with only one ",
             "evaluation metric (you have ",
             length(metrics_to_use),
             " evaluation metrics in 'evaluation_df'). Indeed, metrics have ",
             "distinct orders of magnitude, and so the 'metric_cutoffs' you ",
             " chose are likely to be",
             " appropriate for only one of the metrics, but no the others."), 
             call. = FALSE)
      }
      
      optim_index <- sapply(metric_cutoffs,
                            function(cutoff, vals) which(vals >= cutoff)[1],
                            vals = evaluation_df[, metrics_to_use])
      
      evaluation_df[, paste0("optimal_n_", metrics_to_use)] <- FALSE
      evaluation_df[optim_index, paste0("optimal_n_", 
                                                   metrics_to_use)] <- TRUE
      
      optim_n <- list(evaluation_df$n_bioregions[optim_index])
      names(optim_n) <- metrics_to_use
    }
    
    if(criterion == "max"){
      
      if(verbose){
        message(" - Max value method")
      }
      
      optim_index <- lapply(metrics_to_use,
                            function(x, eval_df) {
                              which(eval_df[, x] == max(eval_df[, x]))
                            }, eval_df = evaluation_df)
      names(optim_index) <- metrics_to_use
      
      optim_n <- lapply(
        optim_index,
        function(x) unique(evaluation_df$n_bioregions[x]))
      
      print(optim_index)
      
      for(metric in names(optim_n)) {
        evaluation_df[optim_index,
                                 paste0("optimal_n_", metric)] <- TRUE
      }
    }
    
    if(criterion == "min"){
      
      if(verbose){
        message(" - Min value method")
      }
      
      optim_index <- lapply(metrics_to_use,
                            function(x, eval_df) {
                              which(eval_df[, x] == min(eval_df[, x]))
                            }, eval_df = evaluation_df)
      names(optim_index) <- metrics_to_use
      
      optim_n <- lapply(
        optim_index,
        function(x) unique(evaluation_df$n_bioregions[x]))
      
      for(metric in names(optim_n)) {
        evaluation_df[optim_index, paste0("optimal_n_", metric)] <- 
          TRUE
      }
    }
    
    if(plot){
      ggdf2 <- tidyr::pivot_longer(
        data = as.data.frame(evaluation_df),
        cols = grep("optimal_n_", colnames(evaluation_df)),
        names_to = "variable")
      ggdf2 <- as.data.frame(ggdf2)
      ggdf2$variable <- gsub("optimal_n_", "", ggdf2$variable)
      ggdf2 <- ggdf2[ggdf2$value, ]
      ggdf <- as.data.frame(tidyr::pivot_longer(
        data = evaluation_df, cols = metrics_to_use,
        names_to = "variable"))
      
      if(verbose){
        message("Plotting results...")
      }
      
      if(criterion == "breakpoints"){
        
        names(seg_preds) <- c("n_bioregions", "value", "variable")
        
        if(verbose){
          message("   (the red line is the prediction from the segmented ",
                  "regression)")
        }
        
        p <- ggplot2::ggplot(ggdf, ggplot2::aes(x = .data$n_bioregions,
                                                       y = .data$value)) +
          ggplot2::geom_line(col = "darkgrey") +
          ggplot2::facet_wrap(~ variable, scales = "free_y") +
          # ggplot2::geom_hline(yintercept = evaluation_df[
          #   evaluation_df$optimal_nclust, eval_metric[1]],
          #                     linetype = 2) +
          ggplot2::geom_vline(data = ggdf2,
                              ggplot2::aes(xintercept = .data$n_bioregions),
                              linetype = 2) +
          ggplot2::theme_bw() +
          ggplot2::geom_line(data = seg_preds,
                             col = "red")
      } else {
        p <- ggplot2::ggplot(ggdf, ggplot2::aes(x = .data$n_bioregions,
                                                       y = .data$value)) +
          ggplot2::geom_line(col = "darkgrey") +
          ggplot2::facet_wrap(~ variable, scales = "free_y") +
          # ggplot2::geom_hline(yintercept = evaluation_df[
          #   evaluation_df$optimal_nclust, eval_metric[1]],
          #                     linetype = 2) +
          ggplot2::geom_vline(data = ggdf2,
                              ggplot2::aes(xintercept = .data$n_bioregions),
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
  evaluation_df = evaluation_df,
  optimal_n = optim_n,
  plot = p)
  
  class(outputs) <- append("bioregion.optimal.n", class(outputs))
  return(outputs)
}
