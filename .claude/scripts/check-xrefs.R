# Book-wide cross-reference integrity check.
#
# A single-chapter render cannot resolve \@ref() targets that live in other
# chapters, and a full render_book is too heavy for a routine check. This does
# the resolution statically: it collects every label the book defines and every
# reference the book uses, then reports the ones that point at nothing.
#
# Two independent link systems are checked, because they fail differently.
#
#   1. bookdown \@ref(id) references, resolved against explicit {#id} anchors,
#      (\#eq:id) labels, and chunk-derived fig:/tab: labels.
#
#   2. Markdown [text](#id) links, resolved against explicit {#id} anchors plus
#      the identifiers pandoc auto-generates from heading text. These are the
#      dangerous ones: an underscore-prefixed chapter is excluded from the build
#      but still defines anchors on disk, so a link into one looks fine to grep
#      and silently 404s in the rendered book. Only non-underscore files are
#      scanned on either side, so such a link is reported.
#
# Regexes use character classes ([{], [\\]) rather than backslash escapes so the
# file survives being written through shells that mangle backslashes.

files <- setdiff(list.files(".", pattern = "[.]Rmd$"),
                 list.files(".", pattern = "^_"))
txt <- lapply(files, readLines, warn = FALSE, encoding = "UTF-8")
names(txt) <- files

grab <- function(x, re) {
  out <- unlist(regmatches(x, gregexpr(re, x, perl = TRUE)))
  if (!length(out)) return(character(0))
  sub(re, "\\1", out, perl = TRUE)
}

hdr_re  <- "[{]#([A-Za-z0-9:_.-]+)[^}]*[}]"
# One backslash in prose, two inside an R string. 40-report's stargazer titles
# are the R-string case, and its report-stargazer-html-anchors chunk makes those
# anchors resolve in HTML as well as LaTeX.
lab_re  <- "[(][\\\\]{1,2}#([A-Za-z]+:[A-Za-z0-9_.-]+)[)]"
ref_re  <- "[\\\\]@ref[(]([^)]+)[)]"
link_re <- "[]][(]#([A-Za-z0-9:_.-]+)[)]"

anchors <- unlist(lapply(txt, function(x) c(grab(x, hdr_re), grab(x, lab_re))))

# fig: and tab: references resolve against chunk labels, not explicit labels.
chunks <- unlist(lapply(txt, function(x) {
  hdr <- grep("^```[{]r[ ,]", x, value = TRUE)
  lab <- sub("^```[{]r[ ,]+([A-Za-z0-9_.-]+).*$", "\\1", hdr)
  keep <- lab != hdr
  # A chunk whose fig.cap is a vector c("..", "..") draws one numbered figure per
  # caption, labelled <label>-1, <label>-2, and so on.
  multi <- keep & grepl("fig[.]cap[ ]*=[ ]*c[(]", hdr)
  ncap <- vapply(hdr[multi], function(h) {
    cap <- sub("^.*fig[.]cap[ ]*=[ ]*c[(]", "", h)
    cap <- sub("[)].*$", "", gsub('"(?:[^"]|[\\\\]")*"', "S", cap, perl = TRUE))
    length(gregexpr("S", cap)[[1]])
  }, 1L)
  c(lab[keep], unlist(Map(function(l, k) paste0(l, "-", seq_len(k)), lab[multi], ncap)))
}))
defined <- unique(c(anchors, paste0("fig:", chunks), paste0("tab:", chunks)))

used <- lapply(txt, function(x) unique(grab(x, ref_re)))

bad <- 0L
for (f in names(used)) {
  miss <- setdiff(used[[f]], defined)
  if (length(miss)) {
    bad <- bad + length(miss)
    cat(f, ":\n", sep = "")
    for (m in miss) cat("   DANGLING  ref(", m, ")\n", sep = "")
  }
}
cat("\nXREF: ", length(unlist(used)), " references, ", length(defined),
    " labels defined, ", bad, " dangling\n", sep = "")
cat(if (bad == 0) "XREF: PASS\n" else "XREF: FAIL\n")

# --- Markdown [text](#id) links -------------------------------------------
#
# Approximates pandoc's auto_identifiers: strip inline markup, drop everything
# up to the first letter, lowercase, keep alphanumerics and -_. , turn runs of
# whitespace into single hyphens. A heading carrying an explicit {#id} uses that
# id instead and contributes no auto-id.

auto_id <- function(heading) {
  h <- sub("^#+[[:space:]]*", "", heading)
  h <- sub("[[:space:]]*[{][^}]*[}][[:space:]]*$", "", h)   # trailing attributes
  h <- gsub("[[]([^]]*)[]][(][^)]*[)]", "\\1", h)           # [text](target)
  h <- gsub("[$][^$]*[$]", "", h)                            # inline math
  h <- gsub("[`*_]", "", h)                                  # code and emphasis
  h <- sub("^[^[:alpha:]]*", "", h)                          # up to first letter
  h <- tolower(h)
  h <- gsub("[^[:alnum:][:space:]._-]", "", h)
  h <- gsub("[[:space:]]+", "-", trimws(h))
  h
}

headings <- lapply(txt, function(x) grep("^#{1,6}[[:space:]]", x, value = TRUE))
auto <- unlist(lapply(headings, function(hs) {
  hs <- hs[!grepl("[{]#", hs)]        # explicit anchors already collected
  vapply(hs, auto_id, character(1), USE.NAMES = FALSE)
}))
targets <- unique(c(anchors, auto[nzchar(auto)]))

# pandoc disambiguates repeated identifiers with -1, -2, ... suffixes.
resolves <- function(id) {
  id %in% targets || grepl("-[0-9]+$", id) && sub("-[0-9]+$", "", id) %in% targets
}

links <- lapply(txt, function(x) unique(grab(x, link_re)))
badlink <- 0L
for (f in names(links)) {
  miss <- Filter(function(id) !resolves(id), links[[f]])
  if (length(miss)) {
    badlink <- badlink + length(miss)
    cat(f, ":\n", sep = "")
    for (m in miss) cat("   DEAD LINK  [...](#", m, ")\n", sep = "")
  }
}
cat("\nLINKS: ", length(unlist(links)), " markdown links, ", length(targets),
    " anchors reachable, ", badlink, " dead\n", sep = "")
cat(if (badlink == 0) "LINKS: PASS\n" else "LINKS: FAIL\n")

if (bad > 0 || badlink > 0) quit(status = 1)
