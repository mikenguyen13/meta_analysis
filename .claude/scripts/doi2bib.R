# Verified citations from Crossref.
#   Rscript .claude/scripts/doi2bib.R doi <doi> <key>      -> BibTeX entry (fails if DOI unknown)
#   Rscript .claude/scripts/doi2bib.R search "<free text>"  -> top 5 matches with DOIs
args <- commandArgs(trailingOnly = TRUE)
ua <- "meta-analysis-book/1.0 (bibliography check)"
get <- function(url) {
  h <- curl::new_handle(useragent = ua, followlocation = TRUE)
  r <- curl::curl_fetch_memory(url, handle = h)
  if (r$status_code != 200) return(NULL)
  rawToChar(r$content)
}
if (args[1] == "search") {
  q <- utils::URLencode(args[2], reserved = TRUE)
  j <- get(sprintf("https://api.crossref.org/works?query.bibliographic=%s&rows=5&select=DOI,title,author,issued,container-title,type", q))
  it <- jsonlite::fromJSON(j, simplifyVector = FALSE)$message$items
  for (x in it) cat(sprintf("%s | %s | %s | %s | %s | %s\n", x$DOI,
    tryCatch(x$issued$`date-parts`[[1]][[1]], error = function(e) NA),
    tryCatch(x$author[[1]]$family, error = function(e) NA),
    tryCatch(x$title[[1]], error = function(e) NA),
    tryCatch(x$`container-title`[[1]], error = function(e) NA), x$type))
} else if (args[1] == "doi") {
  b <- get(sprintf("https://api.crossref.org/works/%s/transform/application/x-bibtex", args[2]))
  if (is.null(b) || !grepl("@", b)) { message("DOI NOT FOUND: ", args[2]); quit(status = 1) }
  b <- gsub("<[^>]+>", "", b)
  b <- gsub("[[:space:]]+", " ", trimws(b))
  if (length(args) >= 3) b <- sub("^@([A-Za-z]+)[{][^,]*,", sprintf("@\\1{%s,", args[3]), b)
  # one field per line
  b <- gsub(", ([A-Za-z]+)=[{]", ",\n  \\1={", b)
  b <- sub(" [}]$", "\n}", b)
  cat(b, "\n\n", sep = "")
}
