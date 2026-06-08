# =============================================================================
# psych_prior.R  —  The psychology replication prior (post-SCORE M1, 2026).
#
# A portable, self-contained description of the go-forward prior for predicting
# the replication effect of a new general-interest psychology study. You do NOT
# need brms, the fitted model, or any serialized object to use it: the prior is
# given below as closed-form expressions you can paste straight into your own
# model. A 12000-draw option (exact joint posterior) is also provided.
#
# Model (Fisher-z scale): for a study with original effect z_orig and a
# replication of sample size n_repli,
#     z_repli ~ Normal(intercept + slope * z_orig + u, sqrt(sigma^2 + se_repli^2))
# where u ~ Normal(0, sd_corpus) is a new-corpus offset and
# se_repli = 1 / sqrt(n_repli - 3). r-scale effects are r = tanh(z).
#
# Runs in base R. No dependencies for the closed-form path.
# =============================================================================

# -----------------------------------------------------------------------------
# THE PRIOR  (post-SCORE M1 fit; n = 251 paired effects across 5 corpora)
#
# Each parameter's posterior was summarized by the closed-form family that best
# fit its 12000 draws (selected by AIC + Kolmogorov-Smirnov; see analysis.R / the
# replication package). Copy these straight into your own model.
#
#   intercept  ~  Normal(mean = -0.0449, sd = 0.0321)
#   slope      ~  Normal(mean =  0.6197, sd = 0.0419)        # coefficient on z_orig
#   sigma      ~  Lognormal(meanlog = -1.6484, sdlog = 0.0543)   # residual SD (> 0)
#   sd_corpus  ~  Gamma(shape = 1.1313, rate = 27.4129)          # between-corpus SD (> 0)
#
# Note: the closed forms are best-fit *marginal* summaries (they treat the four
# parameters as independent). This matches how a prior is normally specified and
# is what the brms prior below encodes. For exact propagation of the joint
# posterior (preserving parameter correlations), use the 12000-draw path further
# down with psych_posterior_draws_2026.rds.
# -----------------------------------------------------------------------------

PSYCH_PRIOR_2026 <- list(
  intercept = list(family = "normal",    mean = -0.0449, sd    = 0.0321),
  slope     = list(family = "normal",    mean =  0.6197, sd    = 0.0419),
  sigma     = list(family = "lognormal", meanlog = -1.6484, sdlog = 0.0543),
  sd_corpus = list(family = "gamma",     shape = 1.1313, rate  = 27.4129)
)

# Paste-ready brms priors (same closed forms), for refitting the M1 spec
# z_repli | se(se_repli, sigma = TRUE) ~ z_orig + (1 | corpus):
#
#   library(brms)
#   psych_prior_brms <- c(
#     set_prior("normal(-0.0449, 0.0321)",      class = "Intercept"),
#     set_prior("normal(0.6197, 0.0419)",       class = "b", coef = "z_orig"),
#     set_prior("lognormal(-1.6484, 0.0543)",   class = "sigma"),
#     set_prior("gamma(1.1313, 27.4129)",       class = "sd")
#   )
#   fit <- brm(formula, data = your_data, prior = psych_prior_brms)

print_psych_prior <- function() {
  cat("Psychology replication prior (post-SCORE M1, 2026)\n")
  cat("  intercept ~ Normal(-0.0449, 0.0321)\n")
  cat("  slope     ~ Normal( 0.6197, 0.0419)   [coefficient on z_orig]\n")
  cat("  sigma     ~ Lognormal(-1.6484, 0.0543)\n")
  cat("  sd_corpus ~ Gamma(shape = 1.1313, rate = 27.4129)\n")
}

# -----------------------------------------------------------------------------
# OPTION A — closed-form forecast (portable, base R, no fitted object)
#
# Forecast the replication effect of a new psychology study from its original
# effect size and the replication's sample size. Simulates parameters from the
# closed-form prior above.
#
#   r_orig   original (published) effect, on the correlation scale
#   n_repli  planned/observed replication sample size
#   n_draws  Monte Carlo draws
# Returns a list with the posterior-predictive mean and 95% interval, on both
# the Fisher-z and r scales.
# -----------------------------------------------------------------------------
forecast_replication <- function(r_orig, n_repli, n_draws = 20000, seed = 20260504) {
  set.seed(seed)
  p <- PSYCH_PRIOR_2026
  z_orig   <- atanh(min(max(r_orig, -0.9999), 0.9999))
  se_repli <- 1 / sqrt(max(n_repli - 3, 1))

  intercept <- rnorm(n_draws, p$intercept$mean, p$intercept$sd)
  slope     <- rnorm(n_draws, p$slope$mean,     p$slope$sd)
  sigma     <- rlnorm(n_draws, p$sigma$meanlog, p$sigma$sdlog)
  sd_corpus <- rgamma(n_draws, shape = p$sd_corpus$shape, rate = p$sd_corpus$rate)
  u         <- rnorm(n_draws, 0, sd_corpus)              # new-corpus offset

  mu     <- intercept + slope * z_orig + u
  z_pred <- rnorm(n_draws, mu, sqrt(sigma^2 + se_repli^2))

  z_q <- quantile(z_pred, c(0.025, 0.5, 0.975))
  list(z_mean = mean(z_pred),
       z_95 = unname(z_q[c(1, 3)]),
       r_median = tanh(unname(z_q[2])),
       r_95 = tanh(unname(z_q[c(1, 3)])))
}

# -----------------------------------------------------------------------------
# OPTION B — 12000-draw forecast (exact joint posterior)
#
# Uses the saved posterior draws (psych_posterior_draws_2026.rds, shipped beside
# this file). Preserves parameter correlations. Identical interface to Option A.
# -----------------------------------------------------------------------------
forecast_replication_draws <- function(r_orig, n_repli,
                                       draws_file = "psych_posterior_draws_2026.rds",
                                       seed = 20260504) {
  set.seed(seed)
  draws <- as.data.frame(readRDS(draws_file))
  z_orig   <- atanh(min(max(r_orig, -0.9999), 0.9999))
  se_repli <- 1 / sqrt(max(n_repli - 3, 1))

  u      <- rnorm(nrow(draws), 0, draws$sd_corpus)       # new-corpus offset
  mu     <- draws$intercept + draws$slope * z_orig + u
  z_pred <- rnorm(nrow(draws), mu, sqrt(draws$sigma^2 + se_repli^2))

  z_q <- quantile(z_pred, c(0.025, 0.5, 0.975))
  list(z_mean = mean(z_pred),
       z_95 = unname(z_q[c(1, 3)]),
       r_median = tanh(unname(z_q[2])),
       r_95 = tanh(unname(z_q[c(1, 3)])))
}

# -----------------------------------------------------------------------------
# Demonstration (runs when sourced interactively or via Rscript)
# -----------------------------------------------------------------------------
if (sys.nframe() == 0 || identical(environment(), globalenv())) {
  print_psych_prior()
  cat("\nExample: a new study reports r = 0.40, replication n = 100.\n")
  a <- forecast_replication(r_orig = 0.40, n_repli = 100)
  cat(sprintf("  Closed-form : predicted r = %.3f, 95%% PI [%.3f, %.3f]\n",
              a$r_median, a$r_95[1], a$r_95[2]))
  if (file.exists("psych_posterior_draws_2026.rds")) {
    b <- forecast_replication_draws(r_orig = 0.40, n_repli = 100)
    cat(sprintf("  12000-draw  : predicted r = %.3f, 95%% PI [%.3f, %.3f]\n",
                b$r_median, b$r_95[1], b$r_95[2]))
  }
}
