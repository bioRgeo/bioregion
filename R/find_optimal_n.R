#' Search for an optimal number of clusters in a list of partitions 
#'
#' This function aims to optimize one or several criteria on a set of 
#' ordered partitions. It is typically used to find one or more optimal 
#' cluster counts on hierarchical trees to cut or ranges of partitions 
#' from k-means or PAM. Users should exercise caution in other cases 
#' (e.g., unordered partitions or unrelated partitions).
#' 
#' @param partitions A `bioregion.partition.metrics` object (output from 
#' [bioregionalization_metrics()]) or a `data.frame` with the first two 
#' columns named `K` (partition name) and `n_clusters` (number of clusters), 
#' followed by columns with numeric evaluation metrics.
#' 
#' @param metrics_to_use A `character` vector or single string specifying 
#' metrics in `partitions` for calculating optimal clusters. Defaults to 
#' `"all"` (uses all metrics).
#' 
#' @param criterion A `character` string specifying the criterion to identify 
#' optimal clusters. Options include `"elbow"`, `"increasing_step"`, 
#' `"decreasing_step"`, `"cutoff"`, `"breakpoints"`, `"min"`, or `"max"`. 
#' Defaults to `"elbow"`. See Details.
#' 
#' @param step_quantile For `"increasing_step"` or `"decreasing_step"`, 
#' specifies the quantile of differences between consecutive partitions as 
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
#' @return
#' A `list` of class `bioregion.optimal.n` with these elements:
#' \itemize{
#' \item{`args`: Input arguments.}
#' \item{`evaluation_df`: The input evaluation `data.frame`, appended with 
#' `boolean` columns for optimal cluster counts.}
#' \item{`optimal_nb_clusters`: A `list` with optimal cluster counts for each 
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
#' partitions. Users specify `step_quantile` or `step_levels`.}
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
#' comnet <- mat_to_net(comat)
#' 
#' dissim <- dissimilarity(comat, metric = "all")
#'
#' # User-defined number of clusters
#' tree1 <- hclu_hierarclust(dissim,
#'                           n_clust = 2:15)
#' tree1
#' 
#' a <- bioregionalization_metrics(tree1,
#'                    dissimilarity = dissim,
#'                    net = comnet,
#'                    species_col = "Node2",
#'                    site_col = "Node1",
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
#'
#' @importFrom stats predict quantile
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes_string geom_line facet_wrap geom_vline
#' @importFrom ggplot2 theme_bw
#' @importFrom rlang .data
#' 
#' @export
find_optimal_n <- function(partitions,
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
      stop(paste0("partitions should be the output object from ", 
                  "bioregionalization_metrics() ",
                  "or a data.frame"), 
           call. = FALSE) 
    } else {
      if(all(colnames(partitions)[1:2] == c("K", "n_clusters")) &
         ncol(partitions) > 2) {
        if(all(sapply(partitions[, 3:ncol(partitions)],
                      is.numeric))) {
          partitions <- list(
            args = list(eval_metric = colnames(partitions)[3: ncol(partitions)]),
                        evaluation_df = partitions)
        } else {
          stop(paste0("Your partition data.frame contains non numeric ",
                      "columns. Only numeric columns are expected after the ",
                      "first two columns"), 
               call. = FALSE) 
        }
        
      } else {
        stop(paste0("partitions should be the output object from ",
                    "bioregionalization_metrics() ",
                    "or a data.frame with the first two columns named as ",
                    "'K' & 'n_clusters' and the following columns being the ",
                    "evaluation metrics"), 
             call. = FALSE) 
        
      }
    }
  }
  
  controls(args = metrics_to_use, data = NULL, type = "character_vector")
  if("all" %in% metrics_to_use) {
    metrics_to_use <- partitions$args$eval_metric
  } else if(any(!(metrics_to_use %in% colnames(partitions$evaluation_df)))) {
    stop(paste0("metrics_to_use should exist in the evaluation table."),
         call. = FALSE)
  }
  
  controls(args = criterion, data = NULL, type = "character")
  if (!(criterion %in% c("elbow", "increasing_step", "decreasing_step", 
                                "cutoff", "breakpoints", "min", "max"))) {
    stop(paste0("Please choose criterion from the following:\n",
                "elbow, increasing_step, decreasing_step, cutoff, breakpoints,",
                " min or max"), 
         call. = FALSE)   
  }
  
  # Verifying that metrics vary, otherwise we remove them
  nvals_metrics <- sapply(lapply(partitions$evaluation_df[, metrics_to_use,
                                                          drop = FALSE],
                               unique), length)
  if (any(nvals_metrics == 1)) {
    exclude_metrics <- names(nvals_metrics[which(nvals_metrics == 1)])
    warning(paste0("Metrics ", 
                   paste0(exclude_metrics, collapse = ", "),
                   " did not vary in partitions, so they were removed."))
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
    stop(paste0("The selected partition metrics did not vary sufficiently ",
                "in input. Please check your partition metrics or increase ",
                "your range of partitions when computing ",
                "bioregionalization_metrics()"), 
         call. = FALSE)
  }
  
  #print(metrics_to_use)
  
  if(nrow(partitions$evaluation_df) <= 4){
    stop(paste0("The number of partitions is too low (<=4) ",
                "for this function to work properly"),
         call. = FALSE) 
         
  } else{
    message(paste0("Number of partitions: ",
                   nrow(partitions$evaluation_df), "\n"))
    
    if(criterion %in% c("elbow", "increasing_step", "decreasing_step",
                        "breakpoints")) {
      if(nrow(partitions$evaluation_df) <= 8) {
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
    
    controls(args = plot, type = "boolean")
    
    message(paste0("Searching for potential optimal number(s) of clusters ",
                   "based on the ",
                   criterion, 
                   " method"))
    
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
        }, eval_df = partitions$evaluation_df, n_breaks = n_breakpoints)
      
      
      optim_n <- lapply(seg_res, function(x) x[[1]])
      names(optim_n) <- metrics_to_use
      
      #Rounding to get the closest partition
      optim_n <- lapply(optim_n, round)
      if(any(!(na.omit(unlist(optim_n)) %in% partitions$evaluation_df$n_clusters))) {
        message("Exact break point not in the list of partitions: finding the",
                " closest partition...\n")
        for(m in names(optim_n)) {
          if(any(!(optim_n[[m]] %in% partitions$evaluation_df$n_clusters))) {
            for(cutoff in  optim_n[[m]][
              which(!(optim_n[[m]] %in% 
                      partitions$evaluation_df$n_clusters))]) {
              optim_n[[m]][which(optim_n[[m]] == cutoff)] <- 
                partitions$evaluation_df$n_clusters[
                  which.min(abs(
                    partitions$evaluation_df$n_clusters - 
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
        stop(paste0("Criterion 'cutoff' should probably be used with only one ",
             "evaluation metric (you have ",
             length(metrics_to_use),
             " evaluation metrics in 'partitions'). Indeed, metrics have ",
             "distinct orders of magnitude, and so the 'metric_cutoffs' you ",
             " chose are likely to be",
             " appropriate for only one of the metrics, but no the others."), 
             call. = FALSE)
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
