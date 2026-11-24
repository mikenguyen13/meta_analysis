# Meta-Analysis: A Practical and Modern Guide with R

[![Build and deploy book](https://github.com/mikenguyen13/meta_analysis/actions/workflows/deploy.yml/badge.svg)](https://github.com/mikenguyen13/meta_analysis/actions/workflows/deploy.yml)
[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-sa/4.0/)

An open textbook on systematic reviews and meta-analysis, with R code that runs on public data
in every chapter.

**Read online:** <https://mikenguyen13.github.io/meta_analysis/>

## Contents

**Part I. Foundations**

1. Introduction to Research Synthesis
2. The Systematic Review Workflow (question, protocol and preregistration, search, screening, PRISMA 2020, data extraction, RoB 2 and ROBINS-I)
3. Effect Sizes (SMD and Hedges' g, log OR/RR/RD, correlations and Fisher's z, conversions, computing from reported statistics, artifact corrections)

**Part II. Core Methods**

4. Common-Effect and Random-Effects Models (tau² estimators, Hartung-Knapp-Sidik-Jonkman, prediction intervals, Q, I², H², GLMMs for binary data)
5. Heterogeneity, Subgroups, and Meta-Regression (permutation tests, multimodel inference, location-scale models)
6. Dependent Effect Sizes (three-level models, robust variance estimation, correlated-and-hierarchical effects)
7. Publication Bias and Small-Study Effects (funnel plots, Egger, trim-and-fill, PET-PEESE, selection models, p-curve, p-uniform, RoBMA)
8. Sensitivity Analysis and Diagnostics (influence, leave-one-out, Baujat, GOSH, cumulative meta-analysis, multiverse)

**Part III. Advanced Methods**

9. Bayesian Meta-Analysis (bayesmeta, brms, heterogeneity priors, model averaging)
10. Network Meta-Analysis (netmeta, multinma, gemtc)
11. Meta-Analytic Structural Equation Modeling (TSSEM, OSMASEM with metaSEM)
12. Multivariate, Dose-Response, and IPD Meta-Analysis
13. Diagnostic Test Accuracy (bivariate and HSROC models)
14. Meta-Analysis of Causal Estimates (RCTs and quasi-experiments, bias adjustment, hierarchical aggregation)

**Part IV. Reporting**

15. Certainty of Evidence (GRADE) and Reporting
16. A Workflow Checklist

## Building the book

The book is a [bookdown](https://bookdown.org/) project rendered with `bookdown::bs4_book`.
Every push to `master` triggers `.github/workflows/deploy.yml`, which installs the packages
listed in `DESCRIPTION` (plus JAGS), renders the book with all code chunks evaluated, and
publishes `_book/` to the `gh-pages` branch. Each deploy commit on `gh-pages` names the source
commit it was built from, and the deployed site carries the same SHA in `SOURCE_COMMIT`.

To build locally:

```r
# install.packages("remotes")
remotes::install_deps(dependencies = TRUE)
bookdown::render_book("index.Rmd", output_format = "bookdown::bs4_book")
```

The Bayesian chapters need a C++ toolchain for Stan; the `RoBMA` and `gemtc` examples need
[JAGS](https://mcmc-jags.sourceforge.io/).

## License

Text and code are licensed under
[CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/). See [LICENSE](LICENSE).

## Citation

See [CITATION.cff](CITATION.cff). Suggested citation:

> Nguyen, M. (2026). *Meta-Analysis: A Practical and Modern Guide with R*. https://mikenguyen13.github.io/meta_analysis/
