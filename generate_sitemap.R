#!/usr/bin/env Rscript
# generate_sitemap.R - writes sitemap.xml into the rendered book output dir and
# gives every page a <link rel="canonical"> to its URL on the site (bs4_book sets
# og:url but no canonical link).
# Run after bookdown::render_book():  Rscript generate_sitemap.R
# The site URL is read from the `url:` field in index.Rmd.

read_field <- function(file, field) {
  lines <- readLines(file, warn = FALSE)
  m <- grep(paste0("^", field, ":"), lines, value = TRUE)
  if (!length(m)) return(NA_character_)
  v <- sub(paste0("^", field, ":[[:space:]]*"), "", m[1])
  v <- gsub("\\:", ":", v, fixed = TRUE)   # bookdown escapes the colon in url:
  gsub('^"|"$', "", trimws(v))
}

out_dir <- "_book"
if (file.exists("_bookdown.yml")) {
  v <- read_field("_bookdown.yml", "output_dir")
  if (!is.na(v)) out_dir <- v
}
base_url <- read_field("index.Rmd", "url")
if (is.na(base_url)) stop("Add a url: field to the index.Rmd YAML header first.")
if (!grepl("/$", base_url)) base_url <- paste0(base_url, "/")
if (!dir.exists(out_dir)) stop("Output dir '", out_dir, "' not found. Render the book first.")

pages <- setdiff(list.files(out_dir, pattern = "[.]html$"), "404.html")
lastmod <- format(file.mtime(file.path(out_dir, pages)), "%Y-%m-%d")
loc <- paste0(base_url, ifelse(pages == "index.html", "", pages))
xml <- c('<?xml version="1.0" encoding="UTF-8"?>',
         '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">',
         paste0("  <url><loc>", loc, "</loc><lastmod>", lastmod, "</lastmod></url>"),
         "</urlset>")
writeLines(xml, file.path(out_dir, "sitemap.xml"))
cat("Wrote", file.path(out_dir, "sitemap.xml"), "-", length(pages), "URLs\n")

# Canonical link on every page, pointing at the same URL the sitemap lists.
canonical <- 0L
for (i in seq_along(pages)) {
  f <- file.path(out_dir, pages[i])
  html <- readLines(f, warn = FALSE, encoding = "UTF-8")
  if (any(grepl('rel="canonical"', html, fixed = TRUE))) next
  head_end <- grep("</head>", html, fixed = TRUE)[1]
  if (is.na(head_end)) next
  html <- append(html, paste0('<link rel="canonical" href="', loc[i], '">'), head_end - 1)
  con <- file(f, open = "wb")
  writeLines(enc2utf8(html), con, useBytes = TRUE)
  close(con)
  canonical <- canonical + 1L
}
cat("Added a canonical link to", canonical, "pages\n")
