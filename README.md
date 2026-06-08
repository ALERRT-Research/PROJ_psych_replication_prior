# In the Meantime: A Bayesian replication prior for psychology

A Bayesian hierarchical prior that forecasts the **replication** effect size of
a new psychology study from its original (published) effect and the planned
replication sample size. The prior is trained on five replication corpora and
validated out-of-sample on SCORE.

This repository is the reproducible analysis and replication package for the
manuscript *"In the Meantime: A Bayesian stopgap for reading psychology while
reform takes hold."*

## The idea

Across five independent replication programs (OSC 2015, Many Labs 1, Many Labs
2, Camerer 2018, and SCORE), the field shows the same pattern: large originals
come back smaller, small originals come back near zero, and paper-level
replication rates cluster near half. Once that pattern has held across five
programs, it stops being news about any one paper and becomes information about
the field, which is what a prior is for. This package treats a new study as one
more datum in a reference class the field has already built.

## Headline result

The pooled forecast covers **95.8%** of held-out SCORE replications inside a
nominal 95% interval, against **45.8%** for the originals' own confidence
intervals and **55.2%** for an honest frequentist prediction interval that
propagates both standard errors. A standard random-effects meta-regression with
the same variance structure covers at 94.8%, so the gain comes from modeling the
replications' residual scatter around the reference class, not from the choice
of Bayesian machinery.

## The model

On the Fisher-z scale, for a study with original effect `z_orig` and a
replication of sample size `n_repli`:

```
z_repli ~ Normal(intercept + slope * z_orig + u, sqrt(sigma^2 + se_repli^2))
u        ~ Normal(0, sd_corpus)          # new-corpus offset
se_repli = 1 / sqrt(n_repli - 3)
```

r-scale effects are recovered with `r = tanh(z)`. The fitted post-SCORE prior
(n = 251 paired effects across five corpora):

```
intercept ~ Normal(-0.0449, 0.0321)
slope     ~ Normal( 0.6197, 0.0419)      # coefficient on z_orig
sigma     ~ Lognormal(-1.6484, 0.0543)
sd_corpus ~ Gamma(shape = 1.1313, rate = 27.4129)
```

## Contents

| File | What it is |
|------|------------|
| `analysis.R` | Self-contained script that reproduces every analysis in the manuscript and supplement, refitting all models from a single data file. No precomputed fits are read. |
| `analysis_data.csv` | 251 paired original/replication effects across the five corpora, on the Fisher-z scale with sampling-variance offsets. |
| `psych_prior.R` | The portable go-forward prior. Closed-form expressions plus paste-ready brms priors; runs in base R with no fitted object required. |
| `psych_posterior_draws_2026.rds` | 12,000 posterior draws for exact joint-posterior forecasting (preserves parameter correlations). |

## Using the prior in your own work

No fitted object needed. The closed-form path runs in base R:

```r
source("psych_prior.R")

# A new study reports r = 0.40; you plan a replication with n = 100.
forecast_replication(r_orig = 0.40, n_repli = 100)
#> predicted r and 95% prediction interval, on both z and r scales
```

For exact propagation of parameter correlations, use the draws-based path
(reads `psych_posterior_draws_2026.rds`, shipped beside the script):

```r
forecast_replication_draws(r_orig = 0.40, n_repli = 100)
```

To refit the full model with these priors, paste the `brms` block from the top
of `psych_prior.R`.

## Reproducing the analysis

Run from the repository folder:

```sh
Rscript analysis.R
```

Requires R (>= 4.6) with `brms` (>= 2.23), `dplyr`, `readr`, `ggplot2`,
`tidyr`, `posterior`, `metafor`, and a working Stan toolchain
(`rstan` / `cmdstanr`). Figures are written to `./figures/`. A seed is set for
reproducibility.

Set `RUN_SENSITIVITY <- FALSE` near the top of `analysis.R` for a fast
main-results-only run (about 20 minutes). `TRUE` also runs the sensitivity
suite (multilab ladder, leave-one-corpus-out, prior sensitivity), about
1 to 1.5 hours total.

## Citation

If you use this software or the replication prior, please cite it using the
metadata in `CITATION.cff`, or the archived release via its Zenodo DOI. DOI: 10.5281/zenodo.20597680

## License

Released under CC-BY-4.0. See `LICENSE`.
