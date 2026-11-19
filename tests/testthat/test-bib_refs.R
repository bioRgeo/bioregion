# bib_refs() needs the package source (the citation style is in vignettes/)
# and pandoc: these tests are skipped on the installed package (R CMD check)
root <- find.package("bioregion")
has_source <- file.exists(file.path(root, "vignettes", "style_citation.csl"))

test_that("bib_refs() formats references with the citation style", {
  skip_if_not(has_source, "Package source not available")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc not available")

  refs <- bib_refs(c("Kreft2010", "Holt2013"))
  expect_type(refs, "character")
  expect_length(refs, 1)

  # One paragraph per reference, sorted by author
  paragraphs <- strsplit(refs, "\n\n", fixed = TRUE)[[1]]
  expect_length(paragraphs, 2)
  expect_match(paragraphs[1], "^Holt BG")
  expect_match(paragraphs[2], "^Kreft H & Jetz W \\(2010\\)")

  # Journal in italics, DOI displayed as a link to doi.org
  expect_match(paragraphs[2], "*Journal of Biogeography*", fixed = TRUE)
  expect_match(paragraphs[2],
               paste0("[10.1111/j.1365-2699.2010.02375.x]",
                      "(https://doi.org/10.1111/j.1365-2699.2010.02375.x)"),
               fixed = TRUE)

  # Duplicated keys are only formatted once
  expect_identical(bib_refs(c("Kreft2010", "Kreft2010")),
                   bib_refs("Kreft2010"))
})

test_that("bib_refs() returns errors for invalid keys", {
  expect_error(bib_refs(character(0)), "character vector",
               fixed = TRUE)
  expect_error(bib_refs(NA_character_), "character vector",
               fixed = TRUE)

  skip_if_not(has_source, "Package source not available")
  expect_error(bib_refs(c("Kreft2010", "Unknown2000")),
               "Citation key(s) not found in REFERENCES.bib: Unknown2000",
               fixed = TRUE)
})
