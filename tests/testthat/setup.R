
quietly <- function(expr) {
  invisible(capture.output(
    suppressMessages(
      suppressWarnings(expr)
    )
  ))
}

skip_if_not_installed_quiet <- function(pkg, version = NULL) {

  ip <- rownames(installed.packages())
  if (!pkg %in% ip) {
    testthat::skip(paste("Package", pkg, "not installed"))
  }
  
  if (!is.null(version)) {
    ver <- tryCatch(utils::packageVersion(pkg), error = function(e) NULL)
    if (is.null(ver) || ver < version) {
      testthat::skip(paste("Package", pkg, "too old"))
    }
  }
  
  invisible(TRUE)
}

