# =============================================================================
# analysis.R  —  "In the Meantime: A Bayesian stopgap for reading psychology
#                 while reform takes hold."
#
# One self-contained script that reproduces every analysis reported in the
# manuscript and supplement, refitting all models from a single data file
# (analysis_data.csv). No precomputed fits are read.
#
# Run from this folder:
#   "C:/Program Files/R/R-4.6.0/bin/Rscript.exe" analysis.R
#
# Requires: brms (>= 2.23), dplyr, readr, ggplot2, tidyr, posterior, metafor, and a
# working Stan toolchain (rstan/cmdstanr). Figures are written to ./figures/.
#
# Sections:
#   1. Load & split the combined data
#   2. OSC-only Layer-1a fit (main_cog)  -> held-out SCORE coverage ~ 91.7%
#   3. Pooled hierarchical M1 fits (pre/post-SCORE) + M_pooled baseline
#   4. Forecast & coverage (pooled 95.8%; frequentist naive 45.8% / honest 55.2%)
#   5. Figures (Figure 1 + supplement S1-S9)
#   6. Sensitivity suite (multilab ladder M0/M2/M3, LOCO, prior sensitivity)
#   7. Export post-SCORE M1 posterior (draws + summaries for psych_prior.R)
#
# RUN_SENSITIVITY <- FALSE runs main results + figures only (~20 min).
# RUN_SENSITIVITY <- TRUE  also runs Section 6 (~1-1.5 h total).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(brms)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(posterior)
  library(metafor)
})

# ---- Resolve script directory so paths work regardless of cwd ----
.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grep("^--file=", .args)])
HERE  <- if (length(.file)) dirname(normalizePath(.file)) else getwd()
setwd(HERE)

FIG_DIR <- file.path(HERE, "figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, showWarnings = FALSE)

RUN_SENSITIVITY <- TRUE   # set FALSE for a fast main-results-only run

set.seed(20260504)
safe_z   <- function(r) atanh(pmin(pmax(r, -0.9999), 0.9999))
save_fig <- function(p, fname, w = 7, h = 4.5) {
  ggsave(file.path(FIG_DIR, fname), p, width = w, height = h, dpi = 150)
  cat("  saved figures/", fname, "\n", sep = "")
}

# =============================================================================
# 1. LOAD & SPLIT
# =============================================================================
d <- read_csv(file.path(HERE, "analysis_data.csv"), show_col_types = FALSE)

corpus_levels <- c("osc", "ml1", "ml2", "camerer2018", "score")

train <- d %>% filter(set == "train") %>%
  mutate(corpus = factor(corpus, levels = corpus_levels))           # n = 155
full  <- d %>%
  mutate(corpus = factor(corpus, levels = corpus_levels))           # n = 251
osc   <- d %>% filter(corpus == "osc")                              # n = 97
score <- d %>% filter(set == "score") %>% mutate(row_id = row_number())  # n = 96
psych_only_ids <- score$row_id[score$pub_expanded == "psychology"]  # n = 69

cat(sprintf("Loaded %d rows: train=%d, score=%d (psych-only=%d), osc=%d\n",
            nrow(d), nrow(train), nrow(score), length(psych_only_ids), nrow(osc)))

# Posterior-predictive forecast helper (z-scale summaries)
forecast_pp <- function(fit, newdata, new_levels = FALSE) {
  pp <- posterior_predict(fit, newdata = newdata, allow_new_levels = TRUE,
                          sample_new_levels = if (new_levels) "gaussian" else "uncertainty")
  data.frame(
    zhat   = apply(pp, 2, mean),
    z_lo80 = apply(pp, 2, quantile, 0.10),
    z_hi80 = apply(pp, 2, quantile, 0.90),
    z_lo95 = apply(pp, 2, quantile, 0.025),
    z_hi95 = apply(pp, 2, quantile, 0.975)
  )
}
cov95 <- function(obs, lo, hi) mean(obs >= lo & obs <= hi)

# =============================================================================
# 2. OSC-ONLY LAYER-1a FIT  (model label: main_cog)
#    z_repli | se(se_repli, sigma=TRUE) ~ z_orig + social   on OSC (n=97).
#    Forecast SCORE with social = 0 (cognitive baseline). Seed 20260409.
#    This is the "OSC original fit" — reproduced from the data, not carried.
# =============================================================================
cat("\n========== 2. OSC-only Layer-1a fit (main_cog) ==========\n")
priors_osc <- c(
  prior(normal(0, 0.5),   class = "Intercept"),
  prior(normal(0.5, 0.5), class = "b", coef = "z_orig"),
  prior(normal(0, 0.3),   class = "b", coef = "social"),
  prior(normal(0, 0.3),   class = "sigma")
)
fit_osc <- brm(
  z_repli | se(se_repli, sigma = TRUE) ~ z_orig + social,
  data = osc, family = gaussian(), prior = priors_osc,
  chains = 4, iter = 2000, warmup = 1000, cores = 4,
  seed = 20260409, silent = 2, refresh = 0
)
score_osc <- score %>% mutate(social = 0)   # SCORE has no discipline tag
fc_osc <- forecast_pp(fit_osc, score_osc)
osc_cov95 <- cov95(score$z_repli, fc_osc$z_lo95, fc_osc$z_hi95)
cat(sprintf("OSC-only held-out SCORE 95%% coverage: %.3f (%d/%d)\n",
            osc_cov95, sum(score$z_repli >= fc_osc$z_lo95 & score$z_repli <= fc_osc$z_hi95),
            nrow(score)))

# =============================================================================
# 3. POOLED HIERARCHICAL M1 FITS + M_pooled BASELINE
#    M1: z_repli | se(se_repli, sigma=TRUE) ~ z_orig + (1 | corpus). Seed 20260504.
# =============================================================================
cat("\n========== 3. Pooled M1 fits (pre/post) + baseline ==========\n")
priors_M1 <- c(
  prior(normal(0, 0.5),   class = "Intercept"),
  prior(normal(0.8, 0.3), class = "b", coef = "z_orig"),
  prior(normal(0.2, 0.1), class = "sigma", lb = 0),
  prior(normal(0, 0.2),   class = "sd", lb = 0)
)
formula_M1 <- bf(z_repli | se(se_repli, sigma = TRUE) ~ z_orig + (1 | corpus))

fit_M1 <- function(data, label) {
  cat("  fitting", label, "(n =", nrow(data), ")\n")
  brm(formula_M1, data = data, prior = priors_M1,
      chains = 4, warmup = 2000, iter = 5000, cores = 4,
      control = list(adapt_delta = 0.999, max_treedepth = 15),
      seed = 20260504, refresh = 0)
}
fit_pre  <- fit_M1(train, "M1 pre-SCORE")
fit_post <- fit_M1(full,  "M1 post-SCORE")

# Level-0 complete-pooling baseline (supplement S5.3)
priors_pooled <- c(
  prior(normal(0, 0.5),   class = "Intercept"),
  prior(normal(0.8, 0.3), class = "b", coef = "z_orig"),
  prior(normal(0.2, 0.1), class = "sigma", lb = 0)
)
fit_baseline <- brm(
  bf(z_repli | se(se_repli, sigma = TRUE) ~ z_orig),
  data = train, prior = priors_pooled,
  chains = 4, warmup = 2000, iter = 5000, cores = 4,
  control = list(adapt_delta = 0.999, max_treedepth = 15),
  seed = 20260504, refresh = 0
)

diag1 <- function(fit, label) {
  div  <- sum(subset(nuts_params(fit), Parameter == "divergent__")$Value)
  rmax <- max(rhat(fit), na.rm = TRUE)
  fx   <- fixef(fit); sp <- summary(fit)$spec_pars
  sdc  <- tryCatch(summary(fit)$random$corpus[1, "Estimate"], error = function(e) NA)
  cat(sprintf("  %-16s div=%d rhat=%.4f  intercept=%.3f slope=%.3f sigma=%.3f sd_corpus=%s\n",
              label, div, rmax, fx["Intercept","Estimate"], fx["z_orig","Estimate"],
              sp[1,"Estimate"], ifelse(is.na(sdc), "NA", sprintf("%.3f", sdc))))
}
cat("\n--- M1 / baseline diagnostics ---\n")
diag1(fit_pre, "M1 pre");  diag1(fit_post, "M1 post"); diag1(fit_baseline, "M_pooled")

# =============================================================================
# 4. FORECAST & COVERAGE
# =============================================================================
cat("\n========== 4. Forecast & coverage ==========\n")
score_eval <- score %>% mutate(corpus = "score")
fc_pooled  <- forecast_pp(fit_pre, score_eval, new_levels = TRUE)

pooled_cov95 <- cov95(score$z_repli, fc_pooled$z_lo95, fc_pooled$z_hi95)
pooled_cov95_psych <- cov95(score$z_repli[psych_only_ids],
                            fc_pooled$z_lo95[psych_only_ids], fc_pooled$z_hi95[psych_only_ids])

# Frequentist baselines (closed-form from the SEs)
se_orig  <- 1 / sqrt(pmax(score$n_orig  - 3, 1))
se_repli <- 1 / sqrt(pmax(score$n_repli - 3, 1))
se_pi    <- sqrt(se_orig^2 + se_repli^2)
naive_cov  <- cov95(score$z_repli, score$z_orig - 1.96 * se_orig, score$z_orig + 1.96 * se_orig)
honest_cov <- cov95(score$z_repli, score$z_orig - 1.96 * se_pi,   score$z_orig + 1.96 * se_pi)

# Frequentist parallel to M1: random-effects meta-regression with the same variance
# structure (known replication sampling variance V = se_repli^2, residual heterogeneity
# via ~1|.rid, corpus random intercept via ~1|corpus), REML. The prediction interval for
# a held-out study in a new corpus propagates parameter uncertainty + both variance
# components + the new study's sampling variance. Confirms the held-out calibration is a
# property of the reference-class model, not of Bayesian machinery. Matches brms M1.
train_f  <- train %>% mutate(.rid = row_number())
fit_freq <- rma.mv(yi = z_repli, V = se_repli^2, mods = ~ z_orig,
                   random = list(~ 1 | corpus, ~ 1 | .rid), data = train_f, method = "REML")
Xs_f     <- cbind(1, score$z_orig)
mu_f     <- as.vector(Xs_f %*% coef(fit_freq))
pv_f     <- rowSums((Xs_f %*% vcov(fit_freq)) * Xs_f) + sum(fit_freq$sigma2) + score$se_repli^2
freq_cov <- cov95(score$z_repli, mu_f - 1.96 * sqrt(pv_f), mu_f + 1.96 * sqrt(pv_f))
cat(sprintf("Frequentist RE meta-regression PI (M1 structure) coverage: %.3f (slope %.3f, tau_sd %.3f, sigma_sd %.3f, mean 95%% width %.3f)\n",
            freq_cov, coef(fit_freq)[["z_orig"]], sqrt(fit_freq$sigma2[1]), sqrt(fit_freq$sigma2[2]),
            mean(2 * 1.96 * sqrt(pv_f))))

comparison <- data.frame(
  method = c("Pooled Bayesian M1 (psych+health)", "Pooled Bayesian M1 (psych-only)",
             "Frequentist RE meta-regression (M1 structure)",
             "OSC-only (main_cog)", "Frequentist naive CI", "Frequentist honest PI"),
  n = c(nrow(score), length(psych_only_ids), nrow(score), nrow(score), nrow(score), nrow(score)),
  coverage95 = round(c(pooled_cov95, pooled_cov95_psych, freq_cov, osc_cov95, naive_cov, honest_cov), 3)
)
write_csv(comparison, file.path(HERE, "coverage_comparison.csv"))
cat("Coverage comparison:\n"); print(comparison, row.names = FALSE)

# =============================================================================
# 5. FIGURES
# =============================================================================
cat("\n========== 5. Figures ==========\n")
theme_set(theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank()))

# ---- Figure 1: three-panel per-study intervals (manuscript) ----
fig1 <- score %>%
  mutate(
    bayes_lo = fc_pooled$z_lo95, bayes_hi = fc_pooled$z_hi95,
    se_orig = se_orig, se_pi = se_pi,
    naive_lo  = z_orig - 1.96 * se_orig, naive_hi  = z_orig + 1.96 * se_orig,
    honest_lo = z_orig - 1.96 * se_pi,   honest_hi = z_orig + 1.96 * se_pi,
    inside_bayes  = z_repli >= bayes_lo  & z_repli <= bayes_hi,
    inside_honest = z_repli >= honest_lo & z_repli <= honest_hi,
    inside_naive  = z_repli >= naive_lo  & z_repli <= naive_hi
  ) %>% arrange(z_orig) %>% mutate(plot_idx = row_number())

cb <- mean(fig1$inside_bayes); ch <- mean(fig1$inside_honest); cn <- mean(fig1$inside_naive)
long1 <- fig1 %>%
  transmute(plot_idx, z_repli,
            bayes_lo, bayes_hi, bayes_inside = inside_bayes,
            honest_lo, honest_hi, honest_inside = inside_honest,
            naive_lo, naive_hi, naive_inside = inside_naive) %>%
  pivot_longer(-c(plot_idx, z_repli), names_to = c("interval_type", ".value"), names_sep = "_") %>%
  mutate(interval_type = factor(interval_type, levels = c("bayes", "honest", "naive"),
            labels = c(sprintf("Bayesian 95%% CrI (coverage %.2f)", cb),
                       sprintf("Honest 95%% PI (coverage %.2f)", ch),
                       sprintf("Naive 95%% CI (coverage %.2f)", cn))))
icols <- setNames(c("#2c7fb8", "#737373", "#bdbdbd"), levels(long1$interval_type))
p1 <- ggplot(long1, aes(x = plot_idx)) +
  geom_hline(yintercept = 0, color = "grey70", linetype = "dashed", linewidth = 0.3) +
  geom_linerange(aes(ymin = lo, ymax = hi, color = interval_type), linewidth = 0.45, alpha = 0.9) +
  geom_point(data = function(d) d[d$inside, ],  aes(y = z_repli), shape = 16, color = "black", size = 0.9) +
  geom_point(data = function(d) d[!d$inside, ], aes(y = z_repli), shape = 4,  color = "red3", size = 0.9, stroke = 0.5) +
  scale_color_manual(values = icols) +
  facet_wrap(~ interval_type, ncol = 1) +
  labs(x = "Study index (sorted by original Fisher-z)",
       y = "Fisher-z replication effect") +
  coord_cartesian(ylim = c(-2.1, 1.8)) +
  theme(legend.position = "none",
        axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        strip.text = element_text(size = 12, hjust = 0))
# Sized to the page text width (6.5in) so Word does not downscale it: text renders at true 12pt.
ggsave(file.path(FIG_DIR, "figure_intervals_per_study.png"), p1, width = 6.5, height = 7.5, dpi = 300)
cat("  saved figures/figure_intervals_per_study.png\n")

# ---- Figure S1: pre/post posterior 95% CrI overlap ----
params <- c("b_Intercept" = "Intercept", "b_z_orig" = "Slope on z_orig",
            "sigma" = "Residual SD", "sd_corpus__Intercept" = "Corpus-level SD")
sp_rows <- function(fit, variant) {
  dd <- as_draws_df(fit)
  do.call(rbind, lapply(names(params), function(p) data.frame(
    param = params[p], variant = variant, median = median(dd[[p]]),
    q025 = quantile(dd[[p]], 0.025), q975 = quantile(dd[[p]], 0.975), row.names = NULL)))
}
ov <- rbind(sp_rows(fit_pre, "Pre-SCORE (155 rows, 4 corpora)"),
            sp_rows(fit_post, "Post-SCORE (251 rows, 5 corpora)"))
ov$param   <- factor(ov$param, levels = rev(c("Slope on z_orig","Intercept","Residual SD","Corpus-level SD")))
ov$variant <- factor(ov$variant, levels = c("Pre-SCORE (155 rows, 4 corpora)","Post-SCORE (251 rows, 5 corpora)"))
pS1 <- ggplot(ov, aes(x = median, y = param, color = variant)) +
  geom_vline(xintercept = 0, color = "grey70", linetype = "dashed", linewidth = 0.3) +
  geom_linerange(aes(xmin = q025, xmax = q975), position = position_dodge(width = 0.5), linewidth = 1.2) +
  geom_point(position = position_dodge(width = 0.5), size = 2.4) +
  scale_color_manual(values = c("Pre-SCORE (155 rows, 4 corpora)" = "#2c7fb8",
                                "Post-SCORE (251 rows, 5 corpora)" = "#c0504d")) +
  labs(x = "Posterior value (Fisher-z scale)", y = NULL) +
  theme(legend.position = "bottom", legend.title = element_blank())
save_fig(pS1, "figure_slope_overlap.png", w = 9, h = 5)

# ---- Figures S2-S5: posterior predictive checks (titles live in the supplement text) ----
save_fig(pp_check(fit_pre,  ndraws = 100) + xlab("z_repli (Fisher-z)"), "M1_pp_density_pre.png")
save_fig(pp_check(fit_post, ndraws = 100) + xlab("z_repli (Fisher-z)"), "M1_pp_density_post.png")
save_fig(pp_check(fit_pre,  type = "intervals_grouped", group = "corpus", prob = 0.5, prob_outer = 0.95),
         "M1_pp_intervals_pre.png", w = 8, h = 5)
save_fig(pp_check(fit_post, type = "intervals_grouped", group = "corpus", prob = 0.5, prob_outer = 0.95),
         "M1_pp_intervals_post.png", w = 8, h = 5)

# ---- Figures S6-S7: trace plots ----
trace_vars <- c("b_Intercept", "b_z_orig", "sigma", "sd_corpus__Intercept")
save_fig(mcmc_plot(fit_pre,  type = "trace", variable = trace_vars), "M1_trace_pre.png", w = 8, h = 5)
save_fig(mcmc_plot(fit_post, type = "trace", variable = trace_vars), "M1_trace_post.png", w = 8, h = 5)

# ---- Figure S8: per-study held-out forecast ----
fdf <- fig1 %>% transmute(plot_idx, z_obs = z_repli, zhat = fc_pooled$zhat[order(score$z_orig)],
                          lo95 = bayes_lo, hi95 = bayes_hi,
                          lo80 = fc_pooled$z_lo80[order(score$z_orig)],
                          hi80 = fc_pooled$z_hi80[order(score$z_orig)], inside95 = inside_bayes)
pS8 <- ggplot(fdf, aes(x = plot_idx)) +
  geom_linerange(aes(ymin = lo95, ymax = hi95), color = "grey70", linewidth = 0.6) +
  geom_linerange(aes(ymin = lo80, ymax = hi80), color = "grey40", linewidth = 0.9) +
  geom_point(aes(y = z_obs, shape = inside95, color = inside95), size = 1.6) +
  geom_point(aes(y = zhat), color = "steelblue", size = 0.6) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 4)) +
  scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "red3")) +
  labs(x = "SCORE replication (ranked by original Fisher-z)", y = "Fisher-z replication effect") +
  theme(legend.position = "none")
save_fig(pS8, "M1_score_forecast.png", w = 9, h = 4.5)

# ---- Figure S9: PIT histogram ----
pp_pre <- posterior_predict(fit_pre, newdata = score_eval, allow_new_levels = TRUE,
                            sample_new_levels = "gaussian")
pit <- vapply(seq_len(nrow(score)), function(i) mean(pp_pre[, i] <= score$z_repli[i]), numeric(1))
# Calibration as an effect size, not a significance verdict: PIT mean/SD against the
# Uniform(0,1) targets (0.5, 1/sqrt(12) = 0.289) and the 1-Wasserstein distance to
# uniform, plus the mean Bayesian 95% interval width (sharpness). No KS p-value.
.spit  <- sort(pit); .npit <- length(pit)
pit_w1 <- mean(abs(.spit - (seq_len(.npit) - 0.5) / .npit))
.plo   <- apply(pp_pre, 2, quantile, 0.025); .phi <- apply(pp_pre, 2, quantile, 0.975)
cat(sprintf("PIT calibration: mean %.3f (unif 0.500), sd %.3f (unif 0.289), 1-Wasserstein-to-uniform %.3f; mean Bayesian 95%% width %.3f\n",
            mean(pit), sd(pit), pit_w1, mean(.phi - .plo)))
pS9 <- ggplot(data.frame(pit = pit), aes(x = pit)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "white", boundary = 0) +
  geom_hline(yintercept = nrow(score) * 0.05, linetype = "dashed", color = "grey50") +
  labs(x = "PIT value", y = "count")
save_fig(pS9, "M1_score_pit.png", w = 6.5, h = 4)

# =============================================================================
# 6. SENSITIVITY SUITE  (supplement S5.4, S8)
#    Multilab ladder (M0/M2/M3), leave-one-corpus-out, prior sensitivity.
#    The multilab indicator (`multilab`) was dropped from the main spec (M1)
#    after this analysis; M0/M2/M3 are the rejected rungs of the ladder.
# =============================================================================
if (RUN_SENSITIVITY) {
  cat("\n========== 6. Sensitivity suite ==========\n")

  fit_quiet <- function(formula, data, prior, ad = 0.999) {
    brm(formula, data = data, prior = prior, chains = 4, warmup = 2000, iter = 5000,
        cores = 4, control = list(adapt_delta = ad, max_treedepth = 15),
        seed = 20260504, refresh = 0)
  }
  cov_of <- function(fit, new_levels = TRUE) {
    fc <- forecast_pp(fit, score_eval, new_levels = new_levels)
    cov95(score$z_repli, fc$z_lo95, fc$z_hi95)
  }
  # Forecast-based metrics on the held-out SCORE set (optionally a row subset by ids).
  metrics_fc <- function(fc, ids = NULL) {
    obs <- score$z_repli
    if (!is.null(ids)) { fc <- fc[ids, ]; obs <- obs[ids] }
    inside <- sum(obs >= fc$z_lo95 & obs <= fc$z_hi95)
    list(mae   = round(mean(abs(obs - fc$zhat)), 3),
         cov80 = round(mean(obs >= fc$z_lo80 & obs <= fc$z_hi80), 3),
         cov95 = round(mean(obs >= fc$z_lo95 & obs <= fc$z_hi95), 3),
         inside = inside)
  }
  # Convergence diagnostics for one fit: max Rhat, min bulk ESS, divergent transitions.
  diag_row <- function(fit, label) {
    div  <- sum(subset(nuts_params(fit), Parameter == "divergent__")$Value)
    rmax <- round(max(rhat(fit), na.rm = TRUE), 3)
    essb <- round(min(summarise_draws(as_draws_array(fit), "ess_bulk")$ess_bulk, na.rm = TRUE))
    data.frame(fit = label, rhat_max = rmax, ess_bulk_min = essb, divergences = div)
  }

  # ---- Multilab ladder on the pre-SCORE pool ----
  priors_M0 <- c(prior(normal(0,0.5),class="Intercept"), prior(normal(0.8,0.3),class="b",coef="z_orig"),
                 prior(normal(0,0.3),class="b",coef="multilab"), prior(normal(0.2,0.1),class="sigma",lb=0),
                 prior(normal(0,0.2),class="sd",lb=0))
  fit_M0 <- fit_quiet(bf(z_repli | se(se_repli, sigma=TRUE) ~ z_orig + multilab + (1|corpus)), train, priors_M0)
  fit_M2 <- fit_quiet(bf(z_repli | se(se_repli, sigma=TRUE) ~ z_orig + (z_orig|corpus)), train,
                      c(prior(normal(0,0.5),class="Intercept"), prior(normal(0.8,0.3),class="b",coef="z_orig"),
                        prior(normal(0.2,0.1),class="sigma",lb=0), prior(normal(0,0.2),class="sd",lb=0)))
  fit_M3 <- fit_quiet(bf(z_repli | se(se_repli, sigma=TRUE) ~ z_orig * multilab + (1|corpus)), train,
                      c(prior(normal(0,0.5),class="Intercept"), prior(normal(0.8,0.3),class="b",coef="z_orig"),
                        prior(normal(0,0.3),class="b",coef="multilab"),
                        prior(normal(0,0.3),class="b",coef="z_orig:multilab"),
                        prior(normal(0.2,0.1),class="sigma",lb=0), prior(normal(0,0.2),class="sd",lb=0)))
  ladder <- data.frame(
    spec = c("M1 + indicator", "M1 (adopted)", "M2 (z_orig|corpus)", "M3 (z_orig*multilab)"),
    formula = c("z_orig + multilab + (1|corpus)", "z_orig + (1|corpus)",
                "z_orig + (z_orig|corpus)", "z_orig * multilab + (1|corpus)"),
    score_cov95 = round(c(cov_of(fit_M0), pooled_cov95, cov_of(fit_M2), cov_of(fit_M3)), 3))
  write_csv(ladder, file.path(HERE, "sensitivity_multilab_ladder.csv"))
  cat("Multilab ladder:\n"); print(ladder, row.names = FALSE)
  ladder_diag <- rbind(diag_row(fit_M0, "M1 + indicator"),
                       diag_row(fit_M2, "M2 varying-slopes"),
                       diag_row(fit_M3, "M3 interaction"))
  write_csv(ladder_diag, file.path(HERE, "sensitivity_ladder_diag.csv"))
  cat("Ladder diagnostics (M1 settings, adapt_delta = 0.999):\n"); print(ladder_diag, row.names = FALSE)

  # ---- Leave-one-corpus-out (M1 spec): full pool + four drops ----
  drop_map <- c(osc = "Drop OSC", ml1 = "Drop Many Labs 1",
                ml2 = "Drop Many Labs 2", camerer2018 = "Drop Camerer 2018 psych")
  loco_fits <- list(); loco_diag <- list()
  m_full <- metrics_fc(fc_pooled)
  loco <- data.frame(variant = "Full pool", n = nrow(train), mae = m_full$mae,
                     cov80 = m_full$cov80, cov95 = m_full$cov95,
                     replications = sprintf("%d/96", m_full$inside),
                     slope = round(fixef(fit_pre)["z_orig","Estimate"], 3))
  for (drop in names(drop_map)) {
    sub <- train %>% filter(corpus != drop) %>% mutate(corpus = droplevels(corpus))
    fit <- fit_M1(sub, paste0("LOCO drop ", drop))
    loco_fits[[drop]] <- fit
    loco_diag[[drop]] <- diag_row(fit, sprintf("LOCO drop %s (n = %d)", drop, nrow(sub)))
    m <- metrics_fc(forecast_pp(fit, score_eval, new_levels = TRUE))
    loco <- rbind(loco, data.frame(variant = drop_map[[drop]], n = nrow(sub), mae = m$mae,
                  cov80 = m$cov80, cov95 = m$cov95, replications = sprintf("%d/96", m$inside),
                  slope = round(fixef(fit)["z_orig","Estimate"], 3)))
  }
  write_csv(loco, file.path(HERE, "sensitivity_loco.csv"))
  cat("Leave-one-corpus-out (M1):\n"); print(loco, row.names = FALSE)

  # ---- Prior sensitivity (M1 priors scaled tight/loose) ----
  priors_tight <- c(prior(normal(0,0.25),class="Intercept"), prior(normal(0.8,0.15),class="b",coef="z_orig"),
                    prior(normal(0.2,0.05),class="sigma",lb=0), prior(normal(0,0.1),class="sd",lb=0))
  priors_loose <- c(prior(normal(0,1.0),class="Intercept"), prior(normal(0.8,0.6),class="b",coef="z_orig"),
                    prior(normal(0.2,0.2),class="sigma",lb=0), prior(normal(0,0.4),class="sd",lb=0))
  fit_tight <- fit_quiet(formula_M1, train, priors_tight)
  fit_loose <- fit_quiet(formula_M1, train, priors_loose)
  # Robustness-fit convergence diagnostics (S5.5): tight/loose prior fits + LOCO fits.
  fit_diag <- rbind(diag_row(fit_tight, "Prior sensitivity tight priors"),
                    diag_row(fit_loose, "Prior sensitivity loose priors"),
                    do.call(rbind, loco_diag))
  write_csv(fit_diag, file.path(HERE, "sensitivity_fit_diag.csv"))
  cat("Robustness-fit diagnostics (M1):\n"); print(fit_diag, row.names = FALSE)
  # Per-spec posterior means (M1 has no multilab term).
  pm <- function(fit) c(round(fixef(fit)["Intercept","Estimate"], 3),
                        round(fixef(fit)["z_orig","Estimate"], 3),
                        round(summary(fit)$spec_pars[1,"Estimate"], 3),
                        round(summary(fit)$random$corpus[1,"Estimate"], 3))
  prior_means <- data.frame(
    parameter = c("Intercept","Slope on z_orig","Residual SD (sigma)","Corpus-level SD"),
    tight = pm(fit_tight), main = pm(fit_pre), loose = pm(fit_loose))
  write_csv(prior_means, file.path(HERE, "sensitivity_prior_means.csv"))
  cat("Prior-spec posterior means (M1):\n"); print(prior_means, row.names = FALSE)

  fc_tight <- forecast_pp(fit_tight, score_eval, new_levels = TRUE)
  fc_loose <- forecast_pp(fit_loose, score_eval, new_levels = TRUE)
  mt <- metrics_fc(fc_tight); mm <- metrics_fc(fc_pooled); ml <- metrics_fc(fc_loose)
  prior_sens <- data.frame(
    priors = c("tight (SDs halved)", "main", "loose (SDs doubled)"),
    mae   = c(mt$mae,   mm$mae,   ml$mae),
    cov80 = c(mt$cov80, mm$cov80, ml$cov80),
    cov95 = c(mt$cov95, mm$cov95, ml$cov95),
    sd_corpus = round(c(summary(fit_tight)$random$corpus[1,"Estimate"],
                        summary(fit_pre)$random$corpus[1,"Estimate"],
                        summary(fit_loose)$random$corpus[1,"Estimate"]), 3))
  write_csv(prior_sens, file.path(HERE, "sensitivity_priors.csv"))
  cat("Prior sensitivity (M1):\n"); print(prior_sens, row.names = FALSE)

  # ---- Psych-only vs psych-and-health subset coverage (S12) ----
  fc_noml1 <- forecast_pp(loco_fits[["ml1"]], score_eval, new_levels = TRUE)
  subset_cov <- data.frame(
    model = c("Pooled full, psych-and-health", "Pooled full, psych-only",
              "Pooled no-ML1, psych-and-health", "Pooled no-ML1, psych-only",
              "OSC-only, psych-and-health", "OSC-only, psych-only"),
    n = c(96, 69, 96, 69, 96, 69),
    cov95 = round(c(metrics_fc(fc_pooled)$cov95, metrics_fc(fc_pooled, psych_only_ids)$cov95,
                    metrics_fc(fc_noml1)$cov95,  metrics_fc(fc_noml1, psych_only_ids)$cov95,
                    cov95(score$z_repli, fc_osc$z_lo95, fc_osc$z_hi95),
                    cov95(score$z_repli[psych_only_ids], fc_osc$z_lo95[psych_only_ids],
                          fc_osc$z_hi95[psych_only_ids])), 3))
  write_csv(subset_cov, file.path(HERE, "sensitivity_psych_subset.csv"))
  cat("Psych subset coverage (S12):\n"); print(subset_cov, row.names = FALSE)
} else {
  cat("\n(Section 6 sensitivity suite skipped: RUN_SENSITIVITY = FALSE)\n")
}

# =============================================================================
# 7. EXPORT POST-SCORE M1 POSTERIOR
#    The go-forward prior. Writes the 12000-draw posterior and prints summaries
#    used to populate psych_prior.R. See psych_prior.R for the closed forms.
# =============================================================================
cat("\n========== 7. Export post-SCORE M1 posterior ==========\n")
draws <- as_draws_df(fit_post) %>%
  select(intercept = b_Intercept, slope = b_z_orig, sigma = sigma, sd_corpus = sd_corpus__Intercept)
draws <- as.data.frame(draws)
saveRDS(draws, file.path(HERE, "psych_posterior_draws_2026.rds"))

cat("Post-SCORE M1 posterior (n =", nrow(draws), " draws):\n")
cat("  parameter     mean      sd      2.5%     97.5%\n")
for (nm in names(draws)) {
  x <- draws[[nm]]
  cat(sprintf("  %-10s %8.4f %8.4f %8.4f %8.4f\n",
              nm, mean(x), sd(x), quantile(x, 0.025), quantile(x, 0.975)))
}

# brmsprior object (independent-Normal approximation; mirrors psych_prior.R)
psych_prior_2026 <- c(
  set_prior(sprintf("normal(%.4f, %.4f)", mean(draws$intercept), sd(draws$intercept)), class = "Intercept"),
  set_prior(sprintf("normal(%.4f, %.4f)", mean(draws$slope),     sd(draws$slope)),     class = "b", coef = "z_orig"),
  set_prior(sprintf("normal(%.4f, %.4f)", mean(draws$sigma),     sd(draws$sigma)),     class = "sigma", lb = "0"),
  set_prior(sprintf("normal(%.4f, %.4f)", mean(draws$sd_corpus), sd(draws$sd_corpus)), class = "sd", lb = "0")
)
saveRDS(psych_prior_2026, file.path(HERE, "psych_prior_2026.rds"))

cat("\nWrote: psych_posterior_draws_2026.rds, psych_prior_2026.rds, coverage_comparison.csv\n")
cat("Done.\n")
