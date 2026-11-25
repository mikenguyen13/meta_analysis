# A Workflow Checklist {#sec-conclusion}

The preceding chapters each took one part of a meta-analysis apart. This closing chapter puts
the parts back together as a sequence of decisions, flags the mistakes that most often undermine
published meta-analyses, and ends with a compact analysis template that runs the core steps on
one dataset. It is meant to be used as a checklist while planning, analyzing, and reviewing.

## Before the analysis {#sec-conclusion-before}

A meta-analysis inherits every weakness of the review that feeds it. The decisions below should
be fixed in a protocol and registered before the data are seen (Chapter \@ref(sec-systematic-review)).

1. **Question.** State the population, intervention or exposure, comparator, outcomes, and
   study designs. Decide whether the target is a single common effect or the mean and spread of a
   distribution of effects, because that choice drives the model (Chapter \@ref(sec-models)).
2. **Protocol and registration.** Register on PROSPERO or OSF. Specify the primary outcome, the
   effect-size metric, the model and heterogeneity estimator, planned moderators, the handling of
   dependent effects, and the sensitivity analyses. Deviations are allowed; undisclosed deviations
   are not.
3. **Search and screening.** Search several databases plus registries and grey literature, screen
   in duplicate, and record counts for the PRISMA 2020 flow diagram.
4. **Extraction.** Pilot the form, extract in duplicate, and keep the raw statistics (means, SDs,
   counts, test statistics) rather than only computed effect sizes, so that every effect size can
   be recomputed (Chapter \@ref(sec-effect-sizes)).
5. **Risk of bias.** Rate each study with a design-appropriate tool (RoB 2, ROBINS-I, QUADAS-2) and
   plan how the ratings will enter the synthesis.

## Choosing the model {#sec-conclusion-model}

Table \@ref(tab:conclusion-decisions) maps common data situations to the approach this book
recommends and the chapter that develops it.



Table: (\#tab:conclusion-decisions)Data situations, recommended approaches, and where they are covered.

|Situation                                                    |Approach                                                                                                |Chapter |
|:------------------------------------------------------------|:-------------------------------------------------------------------------------------------------------|:-------|
|One effect per study, studies differ in design or population |Random-effects model, REML or Paule-Mandel tau^2, Hartung-Knapp CI, prediction interval                 |4       |
|Very few studies (k < 5)                                     |Bayesian model with an informative or weakly informative prior on tau; report sensitivity               |9       |
|Binary outcome with rare events or zero cells                |Mantel-Haenszel, Peto, or binomial GLMM instead of continuity corrections                               |4       |
|Several effects per study, sampling covariances unknown      |CHE working model with robust variance estimation (CR2, Satterthwaite df)                               |6       |
|Effects nested in studies nested in labs or papers           |Three-level model; RVE on top if the dependence structure is uncertain                                  |6       |
|Moderators of interest                                       |Mixed-effects meta-regression with Knapp-Hartung tests; permutation tests; few pre-specified moderators |5       |
|Suspected selective publication                              |Funnel plot plus a triangulation of PET-PEESE, selection models, and sensitivity analysis               |7       |
|Influential or outlying studies                              |Influence diagnostics, leave-one-out, GOSH; report results with and without                             |8       |
|Three or more competing interventions                        |Network meta-analysis, checking transitivity and consistency                                            |10      |
|Correlation matrices and a theoretical path model            |Two-stage or one-stage MASEM                                                                            |11      |
|Correlated outcomes with known within-study correlations     |Multivariate random-effects model                                                                       |12      |
|Participant-level data available                             |One-stage or two-stage IPD meta-analysis separating within- and across-study interactions               |12      |
|Sensitivity and specificity of a test                        |Bivariate random-effects (Reitsma) or HSROC model                                                       |13      |
|Estimates from quasi-experiments with different estimands    |Pool only comparable estimands; design as moderator; bias-adjusted or hierarchical models               |14      |



Two defaults deserve emphasis. First, report the prediction interval alongside the confidence
interval for the mean whenever a random-effects model is used: the confidence interval describes
uncertainty about the average effect, while the prediction interval describes where the true
effect of a new comparable study is likely to fall, and the two can lead to very different
practical conclusions. Second, prefer methods with honest small-sample behavior: the
Hartung-Knapp-Sidik-Jonkman adjustment for the mean, and small-sample corrected robust variance
estimation for dependent effects.

## Common pitfalls {#sec-conclusion-pitfalls}

The following errors recur in published meta-analyses and in manuscripts sent for review.

- **Double counting.** Using a shared control group several times, or entering several outcomes
  from the same participants as independent effects, understates standard errors
  (Chapter \@ref(sec-dependent-effects)).
- **Mixing metrics without conversion.** Pooling standardized mean differences computed from
  change scores with those from post-test scores, or correlations with partial correlations from
  regressions with different control sets, combines quantities that estimate different things.
- **Reading I² as the amount of heterogeneity.** I² is a ratio that depends on study precision.
  Report tau and the prediction interval to describe how much effects vary.
- **Fishing for moderators.** Testing many moderators on few studies produces false positives.
  Pre-specify, adjust for multiplicity, and treat post hoc findings as hypotheses.
- **Ecological inference.** A study-level association between mean age and effect size says
  nothing reliable about whether older individuals benefit more.
- **Treating trim-and-fill as a correction.** It is a sensitivity analysis resting on a strong
  assumption about why studies are missing. No single publication-bias method should be read
  as the answer.
- **Fixed-effect inference with heterogeneous studies.** A common-effect model applied to
  heterogeneous studies produces confidence intervals that are too narrow for any claim beyond
  the included studies.
- **Ignoring risk of bias.** Averaging high- and low-quality studies without examining whether
  results differ by risk of bias can let the weakest evidence dominate.

## A compact analysis template {#sec-conclusion-template}

The code below runs the core steps of a univariate random-effects meta-analysis from raw data to
diagnostics, using the BCG vaccine trials. Each step corresponds to a chapter, and each would be
expanded in a real project.


``` r
library(metafor)

# 1. Effect sizes: log risk ratios and sampling variances (Chapter 3)
dat <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg,
              data = dat.bcg, slab = paste(author, year, sep = ", "))

# 2. Random-effects model with REML and Hartung-Knapp inference (Chapter 4)
res <- rma(yi, vi, data = dat, method = "REML", test = "knha")
res
#> 
#> Random-Effects Model (k = 13; tau^2 estimator: REML)
#> 
#> tau^2 (estimated amount of total heterogeneity): 0.3132 (SE = 0.1664)
#> tau (square root of estimated tau^2 value):      0.5597
#> I^2 (total heterogeneity / total variability):   92.22%
#> H^2 (total variability / sampling variability):  12.86
#> 
#> Test for Heterogeneity:
#> Q(df = 12) = 152.2330, p-val < .0001
#> 
#> Model Results:
#> 
#> estimate      se     tval  df    pval    ci.lb    ci.ub     
#>  -0.7145  0.1808  -3.9522  12  0.0019  -1.1084  -0.3206  ** 
#> 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# 3. Average effect, confidence interval, and prediction interval on the RR scale
pred <- predict(res, transf = exp, digits = 2)
pred
#> 
#>  pred ci.lb ci.ub pi.lb pi.ub 
#>  0.49  0.33  0.73  0.14  1.76

# 4. Confidence interval for tau^2 and tau
confint(res)
#> 
#>        estimate   ci.lb   ci.ub 
#> tau^2    0.3132  0.1197  1.1115 
#> tau      0.5597  0.3460  1.0543 
#> I^2(%)  92.2214 81.9177 97.6781 
#> H^2     12.8558  5.5303 43.0680

# 5. Moderator: absolute latitude (Chapter 5)
res_mod <- rma(yi, vi, mods = ~ ablat, data = dat, method = "REML", test = "knha")
res_mod
#> 
#> Mixed-Effects Model (k = 13; tau^2 estimator: REML)
#> 
#> tau^2 (estimated amount of residual heterogeneity):     0.0764 (SE = 0.0591)
#> tau (square root of estimated tau^2 value):             0.2763
#> I^2 (residual heterogeneity / unaccounted variability): 68.39%
#> H^2 (unaccounted variability / sampling variability):   3.16
#> R^2 (amount of heterogeneity accounted for):            75.62%
#> 
#> Test for Residual Heterogeneity:
#> QE(df = 11) = 30.7331, p-val = 0.0012
#> 
#> Test of Moderators (coefficient 2):
#> F(df1 = 1, df2 = 11) = 12.5905, p-val = 0.0046
#> 
#> Model Results:
#> 
#>          estimate      se     tval  df    pval    ci.lb    ci.ub     
#> intrcpt    0.2515  0.2839   0.8857  11  0.3948  -0.3735   0.8764     
#> ablat     -0.0291  0.0082  -3.5483  11  0.0046  -0.0472  -0.0111  ** 
#> 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# 6. Influence diagnostics (Chapter 8)
inf <- influence(res)
which(inf$is.infl)
#> integer(0)

# 7. Small-study effects: Egger-type test and trim-and-fill (Chapter 7)
regtest(res)
#> 
#> Regression Test for Funnel Plot Asymmetry
#> 
#> Model:     mixed-effects meta-regression model
#> Predictor: standard error
#> 
#> Test for Funnel Plot Asymmetry: t = -0.7812, df = 11, p = 0.4512
#> Limit Estimate (as sei -> 0):   b = -0.5104 (CI: -1.2123, 0.1915)
trimfill(rma(yi, vi, data = dat))
#> 
#> Estimated number of missing studies on the right side: 1 (SE = 2.4528)
#> 
#> Random-Effects Model (k = 14; tau^2 estimator: REML)
#> 
#> tau^2 (estimated amount of total heterogeneity): 0.3313 (SE = 0.1701)
#> tau (square root of estimated tau^2 value):      0.5756
#> I^2 (total heterogeneity / total variability):   92.14%
#> H^2 (total variability / sampling variability):  12.72
#> 
#> Test for Heterogeneity:
#> Q(df = 13) = 154.6750, p-val < .0001
#> 
#> Model Results:
#> 
#> estimate      se     zval    pval    ci.lb    ci.ub      
#>  -0.6571  0.1785  -3.6805  0.0002  -1.0070  -0.3072  *** 
#> 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

The pooled risk ratio printed in step 3 is about 0.49, a substantial
average protective effect of BCG, but the 95% prediction interval runs from about
0.14 to 1.76: in a new setting comparable
to these trials, the vaccine could be highly protective or have no benefit at all. Step 5 shows
that absolute latitude accounts for about 76% of the between-trial variance
(the pseudo-R² printed by the model): trials farther from the equator report stronger protection. This is the typical shape of a careful
answer. The average effect matters less than the question of where and for whom the effect holds,
and the diagnostics in steps 6 and 7 show how much the conclusion depends on individual studies
and on assumptions about missing ones.

## Reporting checklist {#sec-conclusion-reporting}

Whatever the method, a report should let a reader reproduce the analysis and judge the certainty
of its conclusions (Chapter \@ref(sec-reporting)). At a minimum it states the registered protocol
and any deviations; presents the PRISMA 2020 flow diagram; tabulates study characteristics and
risk-of-bias ratings; names the effect-size metric, model, heterogeneity estimator, and
confidence-interval method; reports the pooled estimate with its confidence interval, tau, and
the prediction interval; describes moderator, sensitivity, and publication-bias analyses,
including those that did not change the conclusion; rates certainty with GRADE; and shares the
extracted data and analysis code.

## Further reading {#sec-conclusion-reading}

The *Cochrane Handbook for Systematic Reviews of Interventions* [@higgins2019cochrane] remains
the most complete practical reference for reviews of interventions, and
@borenstein2021introduction is a clear conceptual introduction to the statistics of meta-analysis.
The documentation and worked examples on the metafor website, and the CRAN Task View on
Meta-Analysis, track current software. Methodological developments appear in *Research Synthesis
Methods*, *Psychological Methods*, *Statistics in Medicine*, and the *Journal of Economic Surveys*.
