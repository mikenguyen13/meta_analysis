# Sourced before every chapter (before_chapter_script in _bookdown.yml).
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  cache = TRUE,
  warning = FALSE,
  message = FALSE,
  echo = TRUE,
  dpi = 150,
  cache.lazy = FALSE,
  autodep = TRUE,
  out.width = "85%",
  fig.align = "center",
  fig.width = 7,
  fig.height = 5
)

# Prefix auto-generated chunk labels with the chapter name so unnamed chunks in
# different chapters never collide.
options(knitr.duplicate.label = "allow")
local({
  ci <- tryCatch(knitr::current_input(), error = function(e) NULL)
  if (!is.null(ci) && nzchar(ci)) {
    prefix <- gsub("[^A-Za-z0-9]+", "-", tools::file_path_sans_ext(basename(ci)))
    prefix <- gsub("^-+|-+$", "", prefix)
    if (nzchar(prefix)) knitr::opts_knit$set(unnamed.chunk.label = prefix)
  }
})

options(
  crayon.enabled = FALSE,
  width = 80,
  digits = 4,
  mc.cores = 2
)

set.seed(2026)

suppressPackageStartupMessages({
  library(metafor)
  library(ggplot2)
})
theme_set(theme_minimal(base_size = 12))

# Chunks that need JAGS (RoBMA, gemtc) carry the chunk option `jags = TRUE`.
# The published build (CI, Ubuntu with JAGS installed) runs them like any other
# chunk and fails loudly if JAGS is missing. Only the local single-chapter check
# (.claude/scripts/check-chapter.sh sets META_BOOK_LOCAL_CHECK=1) skips them,
# because this Windows machine has no JAGS.
knitr::opts_hooks$set(jags = function(options) {
  if (isTRUE(options$jags) && Sys.getenv("META_BOOK_LOCAL_CHECK") == "1" &&
      !isTRUE(tryCatch(suppressMessages(requireNamespace("rjags", quietly = TRUE)),
                       error = function(e) FALSE))) {
    options$eval <- FALSE
  }
  options
})
