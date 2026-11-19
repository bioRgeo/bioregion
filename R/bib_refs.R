#' Format references from REFERENCES.bib for the help files
#'
#' Internal helper used when building the documentation. It formats one or
#' several references from `inst/REFERENCES.bib` with the citation style of
#' the vignettes (`vignettes/style_citation.csl`), so that help files and
#' vignettes display references identically, and references are maintained
#' in a single place (the .bib file).
#'
#' It is meant to be called with inline R code in roxygen blocks, which is
#' evaluated when the documentation is generated with
#' `roxygen2::roxygenise()`:
#'
#' ```
#' #' @references
#' #' `r bib_refs(c("Kreft2010", "Holt2013"))`
#' ```
#'
#' The formatted references are then written in the `.Rd` files, so neither
#' pandoc nor the .bib file are needed when the package is installed or used.
#'
#' @param keys A `character` vector of citation keys from `REFERENCES.bib`.
#'
#' @param root A `character` string, the path to the root of the package
#' source. By default, the folder from which the package was loaded, which
#' is the source folder when the documentation is generated with
#' `roxygen2::roxygenise()` or `devtools::document()`.
#'
#' @return A `character` string in markdown, with one paragraph per
#' reference, sorted as in the bibliography of the vignettes (by author,
#' then year).
#'
#' @details
#' Requires pandoc (shipped with RStudio), and the package source (the
#' citation style is in `vignettes/`, which is not part of the installed
#' package).
#'
#' @author
#' Boris Leroy (\email{leroy.boris@@gmail.com})
#'
#' @noRd
bib_refs <- function(keys,
                     root = find.package("bioregion")) {

  if (!is.character(keys) || length(keys) == 0 || anyNA(keys)) {
    stop("keys must be a character vector of citation keys.", call. = FALSE)
  }

  bib <- file.path(root, "inst", "REFERENCES.bib")
  csl <- file.path(root, "vignettes", "style_citation.csl")
  for (f in c(bib, csl)) {
    if (!file.exists(f)) {
      stop("File not found: ", f, ". bib_refs() needs the package source ",
           "(run it from roxygen2::roxygenise() or devtools::document()).",
           call. = FALSE)
    }
  }

  if (!rmarkdown::pandoc_available()) {
    stop("pandoc is required to format references.", call. = FALSE)
  }

  # Check that all keys exist in the .bib file
  bib_lines <- readLines(bib, encoding = "UTF-8", warn = FALSE)
  bib_keys <- sub("^@[A-Za-z]+[{]([^,]+),.*$", "\\1",
                  grep("^@[A-Za-z]+[{]", bib_lines, value = TRUE))
  missing_keys <- setdiff(keys, bib_keys)
  if (length(missing_keys) > 0) {
    stop("Citation key(s) not found in REFERENCES.bib: ",
         paste(missing_keys, collapse = ", "), call. = FALSE)
  }

  # pandoc formats the keys listed in 'nocite', with the citation style
  input <- tempfile(fileext = ".md")
  output <- tempfile(fileext = ".md")
  on.exit(unlink(c(input, output)), add = TRUE)
  writeLines(c("---",
               paste0("nocite: \"", paste0("@", unique(keys),
                                           collapse = ", "), "\""),
               "---"),
             input, useBytes = TRUE)

  rmarkdown::pandoc_convert(input,
                            to = "commonmark-raw_html",
                            output = output,
                            citeproc = TRUE,
                            options = c("--bibliography", bib,
                                        "--csl", csl,
                                        "--wrap=none"))

  refs <- readLines(output, encoding = "UTF-8", warn = FALSE)
  refs <- refs[nzchar(trimws(refs))]

  paste(refs, collapse = "\n\n")
}
