# meta_analysis — working notes

Bookdown book "Meta-Analysis: A Practical and Modern Guide with R". Published
format is `bookdown::bs4_book` (listed first in `_output.yml`). Chapter order
is the explicit `rmd_files` list in `_bookdown.yml`; `_`-prefixed files are not
built. Site: https://mikenguyen13.github.io/meta_analysis/ , built by
`.github/workflows/deploy.yml` on push to `master` and pushed to `gh-pages`.

## Checking a chapter

```bash
bash .claude/scripts/check-chapter.sh 04-models.Rmd     # knit + citeproc, prints RENDER: PASS/FAIL
"/c/Program Files/R/R-4.4.3/bin/Rscript.exe" .claude/scripts/check-xrefs.R   # book-wide \@ref() check
```

The local machine has no JAGS, so chunks that need it (`RoBMA`, `gemtc`) only
run in CI (Ubuntu, `apt install jags`). Everything else must knit locally.

Full local build (release step only):

```bash
export R_LIBS=C:/Users/miken/AppData/Local/R/win-library/4.4 RSTUDIO_PANDOC="C:/Program Files/Quarto/bin/tools"
"/c/Program Files/R/R-4.4.3/bin/Rscript.exe" -e 'bookdown::render_book("index.Rmd", "bookdown::bs4_book")'
```

## Writing conventions

- Every heading carries a stable id: `# Title {#sec-<slug>}`, `## Sub {#sec-<chapterslug>-<slug>}`.
  Ids are lowercase, hyphenated, unique book-wide, and never renamed once published
  (external tools link to them). Unnumbered headings: `{#sec-x .unnumbered}`.
- Cross-references: `\@ref(sec-...)`, figures `\@ref(fig:chunk-label)`.
- All code runs during the render. No `eval = FALSE`, no `try()` wrappers that hide
  failures. Use public data (`metadat`/`metafor` `dat.*`, `meta`, `netmeta`,
  `metaSEM`, `mada`, `multinma`, `baggr` datasets) or simulated data with `set.seed()`.
- Keep each chunk fast (MCMC: few chains/iterations, still converged; GOSH: subsets).
- Name chunks that make figures and give them `fig.cap`.
- Chunks that need JAGS (`RoBMA`, `gemtc`) carry the chunk option `jags = TRUE`.
  The local check skips them (hook in `_common.R`); CI runs them. Nothing else may
  depend on objects they create.
- Callouts: `::: {.rmdnote}` / `::: {.rmdtip}` / `::: {.rmdwarning}`.

## Citations

Every citation has a DOI verified against Crossref. Get the BibTeX from Crossref,
never type it by hand:

```bash
"/c/Program Files/R/R-4.4.3/bin/Rscript.exe" .claude/scripts/doi2bib.R search "free text"
"/c/Program Files/R/R-4.4.3/bin/Rscript.exe" .claude/scripts/doi2bib.R doi 10.1136/bmj.n71 page2021prisma
```

Keys are `firstauthorYEARfirstword` in lowercase. Books without DOIs (Cochrane
Handbook, Borenstein et al.) may be entered by hand with an ISBN or URL. All entries
live in `references.bib`; `packages.bib` is generated.

## Commits

Work on `master` is fine for this repo. `reading/` (copyrighted PDFs) and
`packages/` (third-party zips) are gitignored and must never be committed.
