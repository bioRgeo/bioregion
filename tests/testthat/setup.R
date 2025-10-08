
quietly <- function(expr) {
  invisible(capture.output(
    suppressMessages(
      suppressWarnings(expr)
    )
  ))
}

