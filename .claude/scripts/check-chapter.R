# Standalone single-chapter validation: knit the chunks of one chapter and run
# pandoc --citeproc against the book's .bib files, so an unresolved citation
# key surfaces as "[WARNING] Citeproc: citation <key> not found".
# Cross-chapter \@ref() targets are checked statically by check-xrefs.R.
main <- function(chapter) {
  stopifnot(file.exists(chapter))
  base <- tools::file_path_sans_ext(basename(chapter))
  stem <- paste0("_check-", base)
  tmp_rmd <- paste0(stem, ".Rmd")
  out_html <- paste0(stem, ".html")
  clean <- function() {
    unlink(c(tmp_rmd, out_html, paste0(stem, c(".knit.md", ".utf8.md"))))
    unlink(paste0(stem, "_files"), recursive = TRUE)
  }
  clean()
  on.exit(clean(), add = TRUE)
  bibs <- c("references.bib", "packages.bib", character(0))
  bibs <- bibs[file.exists(bibs)]
  yaml <- c(
    "---",
    sprintf('title: "check: %s"', base),
    sprintf("bibliography: [%s]", paste(bibs, collapse = ", ")),
    "link-citations: true",
    "---", "",
    "```{r _chk_common, include=FALSE}",
    'source("_common.R")',
    "```", ""
  )
  writeLines(c(yaml, readLines(chapter, warn = FALSE, encoding = "UTF-8")),
             tmp_rmd, useBytes = TRUE)
  rmarkdown::render(
    tmp_rmd,
    output_format = bookdown::html_document2(
      number_sections = FALSE, self_contained = FALSE, theme = NULL
    ),
    output_file = out_html, quiet = FALSE, envir = globalenv()
  )
}
main(commandArgs(trailingOnly = TRUE)[1])
