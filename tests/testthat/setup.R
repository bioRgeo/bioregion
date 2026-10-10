# quietly
quietly <- function(expr) {
  invisible(capture.output(
    suppressMessages(
      suppressWarnings(expr)
    )
  ))
}

# quiet install
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

# expect no error for plots
expect_no_error_plotless <- function(expr) {
  
  old_dev <- dev.cur()
  
  pdf(NULL) 
  on.exit({
    while(dev.cur() > 1) dev.off()
    if(old_dev > 1) dev.set(old_dev)
  }, add = TRUE)
  
  setHook("before.new.device", 
                     function(...) TRUE, action = "replace")
  
  testthat::expect_no_error(eval(expr))
  
}