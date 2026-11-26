# Merge bib/*.bib into references.bib, dropping duplicate keys (first wins),
# then report keys cited in the built chapters that are missing, and entries
# in references.bib that no chapter cites.
split_entries <- function(lines) {
  start <- grep("^[[:space:]]*@", lines)
  if (!length(start)) return(list())
  end <- c(start[-1] - 1, length(lines))
  Map(function(s, e) lines[s:e], start, end)
}
key_of <- function(entry) sub("^[[:space:]]*@[A-Za-z]+[{]([^,]+),.*$", "\\1", entry[1])

ref <- split_entries(readLines("references.bib", encoding = "UTF-8", warn = FALSE))
extra <- unlist(lapply(list.files("bib", "[.]bib$", full.names = TRUE), function(f)
  split_entries(readLines(f, encoding = "UTF-8", warn = FALSE))), recursive = FALSE)
all <- c(ref, extra)
keys <- vapply(all, key_of, "")
dups <- unique(keys[duplicated(tolower(keys))])
keep <- all[!duplicated(tolower(keys))]
keys <- vapply(keep, key_of, "")
o <- order(tolower(keys))
drop_trailing_blank <- function(e) {
  e <- trimws(e, "right")
  while (length(e) && e[length(e)] == "") e <- e[-length(e)]
  e
}
out <- unlist(lapply(keep[o], function(e) c(drop_trailing_blank(e), "")))
writeLines(out, "references.bib", useBytes = TRUE)
cat("references.bib:", length(keep), "entries;", length(dups), "duplicate keys dropped\n")
if (length(dups)) cat("  duplicates:", paste(dups, collapse = ", "), "\n")

files <- yaml::read_yaml("_bookdown.yml")$rmd_files
txt <- unlist(lapply(files, readLines, encoding = "UTF-8", warn = FALSE))
txt <- txt[!grepl("^```", txt)]
cited <- unique(unlist(regmatches(txt, gregexpr("(?<![A-Za-z0-9_.])@[A-Za-z][A-Za-z0-9_:-]*[A-Za-z0-9]", txt, perl = TRUE))))
cited <- sub("^@", "", cited)
cited <- cited[!grepl("^R-", cited)]
missing <- setdiff(cited, keys)
cat("cited keys:", length(cited), "; missing from references.bib:", length(missing), "\n")
if (length(missing)) cat("  ", paste(missing, collapse = ", "), "\n")
unused <- setdiff(keys, cited)
cat("uncited entries:", length(unused), "\n")
