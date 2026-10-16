# Helpers for the workflow tests (test-workflow.R), also used by some unit
# tests.
#
# The workflow tests run the analyses the way users run them, from the
# site-species matrix to the final outputs, and check every hand-off between
# two functions. The checks below therefore look at the content of the objects
# and of their printed outputs, and build the expected values from the objects
# themselves: a renamed slot then makes a test fail, while a legitimate change
# in the values does not require editing the tests.


# Shared data ------------------------------------------------------------------
# 80 sites of fishmat, spread over the whole dataset, without the species that
# do not occur in any of them. One of these sites (Kifissos.Beo) shares no
# species with the others: it is absent from the similarity network, so
# network clustering leaves it without a cluster (NA), a case that every
# function using a bioregionalization must handle.
wf_comat <- fishmat[floor(seq(1, nrow(fishmat), length.out = 80)), ]
wf_comat <- wf_comat[, colSums(wf_comat) > 0]
wf_dissim <- dissimilarity(wf_comat, metric = "Simpson")
wf_sim <- similarity(wf_comat, metric = "Simpson")


# expect_output_contains --------------------------------------------------------
# Checks that the printed output of expr contains each of strings. Runs of
# spaces and line breaks are treated as a single space, in both the output and
# strings, so that the checks do not depend on the exact spacing of cat().
# With whole_line = TRUE, each string must be a complete line of the output.
expect_output_contains <- function(expr, strings, whole_line = FALSE) {
  ex <- substitute(expr)
  env <- parent.frame()
  out <- utils::capture.output(eval(ex, env))

  squish <- function(s) trimws(gsub("[[:space:]]+", " ", s))
  if (whole_line) {
    found <- squish(strings) %in% squish(out)
  } else {
    text <- squish(paste(out, collapse = "\n"))
    found <- vapply(squish(strings),
                    function(s) grepl(s, text, fixed = TRUE),
                    logical(1))
  }

  testthat::expect(
    all(found),
    paste0("The printed output of `", deparse1(ex), "` does not contain ",
           paste0("\"", strings[!found], "\"", collapse = ", "),
           ".\nPrinted output:\n", paste(out, collapse = "\n")))

  invisible(out)
}


# format_n_print -----------------------------------------------------------------
# Numbers as print.bioregion.clusters() shows them: the first 10 values, then
# the number of values left out.
format_n_print <- function(v) {
  if (length(v) > 10) {
    paste0(paste(v[1:10], collapse = " "),
           " ... (with ", length(v) - 10, " more values)")
  } else {
    paste(v, collapse = " ")
  }
}


# optimal_n_lines -----------------------------------------------------------------
# The lines that print() of a find_optimal_n() output must show: each metric
# followed by its optimal numbers of clusters, or "none found".
optimal_n_lines <- function(opt) {
  vapply(names(opt$optimal_n), function(m) {
    n <- opt$optimal_n[[m]]
    n <- n[!is.na(n)]
    paste(m, "-", if (length(n) == 0) "none found" else paste(n, collapse = " "))
  }, character(1), USE.NAMES = FALSE)
}


# expect_valid_bioregion_clusters ------------------------------------------------
# Checks the properties that every bioregion.clusters object must have,
# whatever the algorithm that produced it.
expect_valid_bioregion_clusters <- function(x) {
  testthat::expect_s3_class(x, "bioregion.clusters")

  # clusters: first column with the site (or node) names, then one column per
  # partition, and a node_type attribute with one value per row
  cl <- x$clusters
  testthat::expect_true(is.data.frame(cl))
  testthat::expect_gte(ncol(cl), 2)
  testthat::expect_type(cl[[1]], "character")
  testthat::expect_false(anyNA(cl[[1]]))
  testthat::expect_false(anyDuplicated(cl[[1]]) > 0)
  node_type <- attr(cl, "node_type")
  testthat::expect_length(node_type, nrow(cl))
  testthat::expect_true(all(node_type %in% c("site", "species")))

  # cluster_info: one row per partition, matching the actual partitions
  info <- x$cluster_info
  testthat::expect_true(is.data.frame(info))
  testthat::expect_equal(nrow(info), ncol(cl) - 1)
  testthat::expect_equal(info$partition_name, names(cl)[-1])
  n_observed <- vapply(cl[, -1, drop = FALSE],
                       function(y) length(unique(y[!is.na(y)])),
                       numeric(1))
  testthat::expect_equal(unname(info$n_clust), unname(n_observed))

  # print() shows the numbers of clusters of all partitions (and, for
  # hclu_hierarclust(), all the requested numbers of clusters)
  expected <- paste("Number of clusters:", format_n_print(info$n_clust))
  if (x$name == "hclu_hierarclust" && !is.null(x$args$n_clust)) {
    expected <- c(expected,
                  paste("Number of clusters requested by the user:",
                        format_n_print(x$args$n_clust)))
  }
  expect_output_contains(print(x), expected)

  # summary() shows the number of clusters of the partitions it describes
  # (the first 3 by default)
  n_shown <- min(3, nrow(info))
  expect_output_contains(summary(x),
                         paste("Total clusters:", info$n_clust[seq_len(n_shown)]))

  testthat::expect_no_error(utils::capture.output(str(x)))

  invisible(x)
}


# ssm_metric_names ---------------------------------------------------------------
# The metrics stored in a site_species_metrics() output, recovered from the
# column names of its data.frames: without the identifier columns, without the
# _occ / _abd suffixes, and with the raw counts (n_sb, n_s, ...) grouped under
# "core_terms". With columns = TRUE, the column names themselves.
ssm_metric_names <- function(x, columns = FALSE) {
  parts <- if (is.data.frame(x[[1]])) list(x) else x
  cols <- unique(unlist(lapply(parts, function(p) lapply(p, names))))
  cols <- setdiff(cols, c("species", "site", "bioregion", "chorotype"))
  if (columns) {
    return(cols)
  }
  metrics <- sub("_(occ|abd)$", "", cols)
  metrics[grepl("^[nw]_(sb|s|b|gc|g|c)$", metrics)] <- "core_terms"
  unique(metrics)
}
