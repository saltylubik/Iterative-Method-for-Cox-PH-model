# =============================================================================
#  Iterative Method for Modelling the Distribution of Regression Effects
#  in the Cox Proportional Hazards Model
#
#  Authors : Lubomír Seif & Ivana Malá
#            Prague University of Economics and Business
#
#  Purpose : Full implementation, concentration-coefficient estimators,
#            simulation study, and validity checks on simulation.
#            Supplementary material for the SAfJR 2026 poster.
#
#  Usage   : source("cox_iterative_supplement.R")
#            All outputs are written to ./outputs/ (created automatically).
#
#  Sections:
#    0.  Packages & setup
#    1.  Configuration (edit here to change the simulation design)
#    2.  Data generator
#    3.  Core: bootstrap p-value grid        (bootstrap_pvals_cox)
#    4.  Probability measure                 (pi_from_pbar)
#    5.  Concentration-coefficient estimators
#        5a. CI matching                     (lambda_ci_match)
#        5b. Contrast ratio                  (lambda_contrast)
#    6.  Interval summaries                  (hdi_discrete, quantile_interval)
#    7.  Main wrapper                        (cox_iterative_method)
#    8.  Plot helpers
#    9.  Single-dataset illustration
#   10.  Simulation study
#        10a. Single-replicate runner        (run_one_rep)
#        10b. Full runner                    (run_simulation)
#        10c. Figures
#   11.  Validity checks & summary diagnostics
# =============================================================================


# -----------------------------------------------------------------------------
# 0.  Packages & setup
# -----------------------------------------------------------------------------
required_pkgs <- c("survival", "ggplot2", "dplyr", "tidyr",
                   "kableExtra", "scales", "patchwork")
invisible(lapply(required_pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))

# Output directory (relative to working directory)
OUT_DIR <- "outputs"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
out_path <- function(...) file.path(OUT_DIR, ...)


# -----------------------------------------------------------------------------
# 1.  Configuration
#     Edit this block to change the simulation design without touching the code.
# -----------------------------------------------------------------------------
CONFIG <- list(
  # Single-dataset illustration
  demo_n          = 300,
  demo_beta       = c(log(1.5), log(0.8), 0.3),  # true log-HRs for x1, x2, x3
  demo_cens_rate  = 0.02,
  demo_B          = 1000,   # bootstrap replicates for illustration
  demo_m          = 1000,   # grid points for illustration
  demo_seed       = 42,
  demo_boot_seed  = 123,

  # Simulation study
  sim_n_vec       = c(200, 500, 1000),
  sim_beta_vec    = c(log(1.1), log(1.5), log(2.0)),
  sim_cens_rates  = c(0.02, 0.08),   # targets ≈ 20 % and ≈ 50 % censoring
  sim_R_reps      = 100,             # replicates per condition (use 500 for final)
  sim_B           = 500,             # bootstrap replicates per replicate
  sim_m           = 500,             # grid points per replicate
  sim_base_seed   = 2025,

  # Method settings
  target_mass     = 0.95,            # CI-matching calibration target
  R_contrast      = 50               # contrast-ratio target R
)


# -----------------------------------------------------------------------------
# 2.  Data generator
#     Weibull (shape = 1, i.e. exponential) baseline, exponential censoring.
# -----------------------------------------------------------------------------
sim_cox_data <- function(n         = 300,
                         lambda0   = 0.05,
                         beta      = c(log(1.5), log(0.8), 0.3),
                         cens_rate = 0.02,
                         seed      = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x1  <- rnorm(n)
  x2  <- rbinom(n, 1, 0.5)
  x3  <- rnorm(n)
  X   <- cbind(x1 = x1, x2 = x2, x3 = x3)
  eta <- drop(X %*% beta)

  T_event <- -log(runif(n)) / (lambda0 * exp(eta))
  C       <- rexp(n, rate = cens_rate)
  time    <- pmin(T_event, C)
  status  <- as.integer(T_event <= C)

  data.frame(time = time, status = status, x1 = x1, x2 = x2, x3 = x3)
}


# -----------------------------------------------------------------------------
# 3.  Bootstrap p-value grid
#
#  For each bootstrap resample b = 1, …, B:
#    • Refit the Cox model on the resampled data (case resampling).
#    • Compute the Wald statistic  Z_k^(b) = (β̂^(b) − β_k) / SE^(b)
#      for every grid point β_k.
#    • Record the two-sided p-value  p_k^(b) = 2(1 − Φ(|Z_k^(b)|)).
#
#  Aggregate over B resamples:
#    p̄_k = (1/B) Σ_b p_k^(b)      (mean p-value at grid point k)
#    SE(p̄_k) = SD_b(p_k^(b)) / √B  (Monte Carlo standard error)
# -----------------------------------------------------------------------------
bootstrap_pvals_cox <- function(df,
                                grid,
                                j_name    = "x1",
                                formula   = Surv(time, status) ~ x1 + x2 + x3,
                                ties      = "efron",
                                B         = 1000,
                                min_valid = 200,
                                seed      = NULL,
                                verbose   = TRUE) {
  if (!is.null(seed)) set.seed(seed)

  n <- nrow(df)
  m <- length(grid)
  P <- matrix(NA_real_, nrow = m, ncol = B)

  boot_once <- function(dfb) {
    fitb <- try(coxph(formula, data = dfb, ties = ties), silent = TRUE)
    if (inherits(fitb, "try-error")) return(NULL)
    vcb  <- try(vcov(fitb), silent = TRUE)
    if (inherits(vcb, "try-error")) return(NULL)
    bh   <- try(unname(coef(fitb)[j_name]), silent = TRUE)
    se   <- try(sqrt(vcb[j_name, j_name]), silent = TRUE)
    if (any(sapply(list(bh, se), inherits, "try-error"))) return(NULL)
    if (!is.finite(bh) || !is.finite(se) || se <= 0)     return(NULL)
    Z <- (bh - grid) / se
    2 * (1 - pnorm(abs(Z)))
  }

  if (verbose) pb <- txtProgressBar(min = 0, max = B, style = 3)
  for (b in seq_len(B)) {
    idx       <- sample.int(n, replace = TRUE)
    pvec      <- boot_once(df[idx, , drop = FALSE])
    if (!is.null(pvec)) P[, b] <- pvec
    if (verbose) setTxtProgressBar(pb, b)
  }
  if (verbose) close(pb)

  valid_counts <- rowSums(is.finite(P))
  pbar <- rowMeans(P, na.rm = TRUE)
  sepb <- apply(P, 1, sd, na.rm = TRUE) / sqrt(pmax(valid_counts, 1))
  keep <- valid_counts >= min_valid & is.finite(pbar) & is.finite(sepb)

  list(grid = grid[keep], pbar = pbar[keep], sepb = sepb[keep],
       valid_counts = valid_counts[keep], B = B, j_name = j_name)
}


# -----------------------------------------------------------------------------
# 4.  Probability measure  π(β_k)
#
#   w_k   = exp(λ · p̄_k)    [equivalent to exp(−λ(1 − p̄_k)) up to const]
#   π(β_k) = w_k / Σ_ℓ w_ℓ
#
#  The formulations exp(λ p̄) and exp(−λ(1 − p̄)) are identical after
#  normalisation because exp(−λ) cancels.  The quantity (1 − p̄) measures
#  implausibility, playing the role of "energy" in the Boltzmann analogy
#  (1/λ is the "temperature").
# -----------------------------------------------------------------------------
pi_from_pbar <- function(pbar, lambda) {
  stopifnot(is.numeric(pbar), is.numeric(lambda),
            length(lambda) == 1, lambda >= 0)
  w <- exp(lambda * pbar)
  w[!is.finite(w)] <- 0
  total <- sum(w)
  if (total == 0) return(rep(1 / length(w), length(w)))
  w / total
}


# -----------------------------------------------------------------------------
# 5a. λ estimator: CI matching
#
#  Finds λ such that  Σ_{k : β_k ∈ Wald CI} π(β_k) = target_mass.
#  This aligns the new plausibility distribution with standard Wald-CI
#  coverage, making the result directly comparable to classical inference.
# -----------------------------------------------------------------------------
lambda_ci_match <- function(grid, pbar,
                            bhat, sehat,
                            target_mass = 0.95,
                            z_crit      = 1.96,
                            interval    = c(0, 1000),
                            tol         = 1e-8,
                            maxiter     = 200) {
  inside <- (grid >= bhat - z_crit * sehat) &
            (grid <= bhat + z_crit * sehat)

  if (!any(inside)) {
    message("[lambda_ci_match] No grid points inside Wald CI; returning lambda = 0")
    return(list(lambda = 0, mass = 0, note = "no_grid_in_CI"))
  }

  f <- function(lam) {
    w <- exp(lam * pbar); w <- w / sum(w)
    sum(w[inside]) - target_mass
  }

  f0 <- f(interval[1]); f1 <- f(interval[2])
  if (sign(f0) == sign(f1)) {
    for (upper in c(2000, 5000, 10000)) {
      if (sign(f0) != sign(f(upper))) { interval[2] <- upper; break }
    }
  }

  if (sign(f(interval[1])) == sign(f(interval[2]))) {
    best <- if (abs(f(interval[1])) < abs(f(interval[2]))) interval[1] else interval[2]
    w <- exp(best * pbar); w <- w / sum(w)
    return(list(lambda = best, mass = sum(w[inside]),
                note = "bracket_failed_returning_boundary"))
  }

  root <- uniroot(f, interval = interval, tol = tol, maxiter = maxiter)
  lam  <- root$root
  w    <- exp(lam * pbar); w <- w / sum(w)
  list(lambda = lam, mass = sum(w[inside]), note = "OK")
}


# -----------------------------------------------------------------------------
# 5b. λ estimator: Contrast ratio (closed form)
#
#  Chooses λ so that  π(β̂) / π(β_ref) = R,
#  where β_ref is the median of p̄ (default) or the second-highest p̄ value.
#
#  From  exp(λ · p̄_max) / exp(λ · p̄_ref) = R
#        → λ = log(R) / (p̄_max − p̄_ref)
#
#  Interpretation: the mode is always R times more plausible than the
#  chosen reference point, making λ a "signal-to-noise" scale.
# -----------------------------------------------------------------------------
lambda_contrast <- function(pbar,
                            R   = 50,
                            ref = c("median", "second_best"),
                            eps = 1e-9) {
  ref      <- match.arg(ref)
  pmax_val <- max(pbar, na.rm = TRUE)

  pref <- switch(ref,
    median = median(pbar, na.rm = TRUE),
    second_best = {
      ord <- order(pbar, decreasing = TRUE, na.last = NA)
      pbar[ord[pmin(2L, length(ord))]]
    }
  )

  denom <- pmax_val - pref
  if (!is.finite(denom) || denom < eps) {
    message("[lambda_contrast] Flat p̄ curve; returning lambda = 0")
    return(list(lambda = 0, ratio = 1, note = "flat_pbar"))
  }

  lam   <- log(R) / (denom + eps)
  ratio <- exp(lam * pmax_val) / exp(lam * pref)
  list(lambda = lam, ratio = ratio, note = "OK")
}


# -----------------------------------------------------------------------------
# 6.  Interval summaries
# -----------------------------------------------------------------------------

# Highest Density Interval — narrowest interval covering a given mass.
hdi_discrete <- function(grid, pi, mass = 0.95) {
  stopifnot(length(grid) == length(pi), abs(sum(pi) - 1) < 1e-6)
  ord   <- order(pi, decreasing = TRUE)
  cum   <- cumsum(pi[ord])
  thr   <- pi[ord][which(cum >= mass)[1]]
  in_hdi <- pi >= thr
  c(lo = min(grid[in_hdi]), hi = max(grid[in_hdi]))
}

# Equal-tailed quantile interval.
quantile_interval <- function(grid, pi, mass = 0.95) {
  stopifnot(length(grid) == length(pi))
  lo_tail <- (1 - mass) / 2
  hi_tail <- 1 - lo_tail
  cum <- cumsum(pi)
  c(lo = grid[which(cum >= lo_tail)[1]],
    hi = grid[which(cum >= hi_tail)[1]])
}


# -----------------------------------------------------------------------------
# 7.  Main wrapper: cox_iterative_method()
# -----------------------------------------------------------------------------
cox_iterative_method <- function(df,
                                 j_name        = "x1",
                                 formula       = Surv(time, status) ~ x1 + x2 + x3,
                                 ties          = "efron",
                                 B             = 1000,
                                 m             = 1000,
                                 width_se      = 3,
                                 min_valid     = 200,
                                 lambda_method = c("ci_match", "contrast", "both"),
                                 R_contrast    = 50,
                                 target_mass   = 0.95,
                                 seed          = NULL,
                                 verbose       = TRUE) {

  lambda_method <- match.arg(lambda_method)

  # Step 1 — fit on observed data
  fit0  <- coxph(formula, data = df, ties = ties)
  bhat  <- unname(coef(fit0)[j_name])
  sehat <- sqrt(vcov(fit0)[j_name, j_name])
  D     <- sum(df$status == 1)

  if (verbose) {
    cat(sprintf("\n[cox_iterative_method]\n"))
    cat(sprintf("  Coefficient : %s\n", j_name))
    cat(sprintf("  beta_hat = %.4f   SE = %.4f   Events = %d\n", bhat, sehat, D))
    cat(sprintf("  Search interval: [%.4f, %.4f]\n",
                bhat - width_se * sehat, bhat + width_se * sehat))
    cat(sprintf("  B = %d   m = %d\n\n", B, m))
  }

  # Step 2 — grid
  grid <- seq(bhat - width_se * sehat,
              bhat + width_se * sehat,
              length.out = m)

  # Step 3 — bootstrap p-values
  boot <- bootstrap_pvals_cox(df, grid, j_name = j_name,
                              formula = formula, ties = ties,
                              B = B, min_valid = min_valid,
                              seed = seed, verbose = verbose)

  # Step 4 — concentration coefficient & π(β)
  results <- list()

  if (lambda_method %in% c("ci_match", "both")) {
    ci_out <- lambda_ci_match(boot$grid, boot$pbar, bhat, sehat,
                              target_mass = target_mass)
    pi_ci  <- pi_from_pbar(boot$pbar, ci_out$lambda)
    results[["ci_match"]] <- list(
      lambda     = ci_out$lambda,
      mass_in_CI = ci_out$mass,
      pi         = pi_ci,
      hdi        = hdi_discrete(boot$grid, pi_ci),
      qi         = quantile_interval(boot$grid, pi_ci)
    )
    if (verbose)
      cat(sprintf("  [CI match]  lambda = %.4f   mass in CI = %.4f\n",
                  ci_out$lambda, ci_out$mass))
  }

  if (lambda_method %in% c("contrast", "both")) {
    cr_out <- lambda_contrast(boot$pbar, R = R_contrast)
    pi_cr  <- pi_from_pbar(boot$pbar, cr_out$lambda)
    results[["contrast"]] <- list(
      lambda = cr_out$lambda,
      ratio  = cr_out$ratio,
      pi     = pi_cr,
      hdi    = hdi_discrete(boot$grid, pi_cr),
      qi     = quantile_interval(boot$grid, pi_cr)
    )
    if (verbose)
      cat(sprintf("  [Contrast]  lambda = %.4f   ratio = %.2f\n",
                  cr_out$lambda, cr_out$ratio))
  }

  wald_ci <- c(lo = bhat - 1.96 * sehat, hi = bhat + 1.96 * sehat)

  list(j_name = j_name, bhat = bhat, sehat = sehat, D = D,
       wald_ci = wald_ci, grid = boot$grid, pbar = boot$pbar,
       sepb = boot$sepb, valid_counts = boot$valid_counts,
       B = B, results = results)
}


# -----------------------------------------------------------------------------
# 8.  Plot helpers
# -----------------------------------------------------------------------------

# 8a. Bootstrap p̄(β) curve with pointwise confidence band
plot_pbar_band <- function(out, alpha_band = 0.15,
                           title = expression(bar(p)(beta)~"bootstrap curve")) {
  z975 <- qnorm(0.975)
  dfp  <- data.frame(
    beta = out$grid,
    pbar = out$pbar,
    lo   = pmax(0, out$pbar - z975 * out$sepb),
    hi   = pmin(1, out$pbar + z975 * out$sepb)
  )

  ggplot(dfp, aes(beta, pbar)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#004e92", alpha = alpha_band) +
    geom_line(linewidth = 0.9, colour = "#004e92") +
    geom_vline(xintercept = out$bhat, linetype = "dashed",
               colour = "#e05c00", linewidth = 0.7) +
    geom_hline(yintercept = 0.05, linetype = "dotted",
               colour = "grey40", linewidth = 0.6) +
    annotate("text", x = out$bhat, y = max(dfp$pbar) * 1.04,
             label = expression(hat(beta)), colour = "#e05c00", size = 4) +
    scale_x_continuous(name = expression(beta~"(log-HR)")) +
    scale_y_continuous(name = expression(bar(p)(beta)),
                       limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.06))) +
    ggtitle(title) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 12, face = "bold"),
          panel.grid.minor = element_blank())
}

# 8b. π(β) overlay for multiple λ strategies
plot_pi_overlay <- function(out,
                            which_methods = names(out$results),
                            title = expression(pi(beta)~"under different"~lambda~"estimators")) {
  colours <- c(ci_match = "#004e92", contrast = "#e05c00")

  df_all <- do.call(rbind, lapply(which_methods, function(nm)
    data.frame(method = nm, beta = out$grid, pi = out$results[[nm]]$pi)
  ))

  ci_lo <- out$wald_ci["lo"]; ci_hi <- out$wald_ci["hi"]

  ggplot(df_all, aes(beta, pi, colour = method, group = method)) +
    annotate("rect", xmin = ci_lo, xmax = ci_hi,
             ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.4) +
    geom_line(linewidth = 1.0) +
    geom_vline(xintercept = out$bhat, linetype = "dashed",
               colour = "grey30", linewidth = 0.7) +
    scale_colour_manual(
      values = colours[which_methods],
      labels = c(ci_match = "CI matching", contrast = "Contrast ratio"),
      name   = expression(lambda~"estimator")
    ) +
    scale_x_continuous(name = expression(beta~"(log-HR)")) +
    scale_y_continuous(name = expression(pi(beta)),
                       expand = expansion(mult = c(0, 0.06))) +
    ggtitle(title) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 12, face = "bold"),
          legend.position = "bottom",
          panel.grid.minor = element_blank())
}

# 8c. Simulation summary plots
plot_sim_summary <- function(sim_df, metric, ylab,
                             title = NULL) {
  ggplot(sim_df, aes(x = factor(n), y = .data[[metric]],
                     colour = method, group = method)) +
    stat_summary(fun = mean, geom = "line",  linewidth = 0.9) +
    stat_summary(fun = mean, geom = "point", size = 2.5) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    facet_wrap(~ cens_pct, labeller = label_both, scales = "free_y") +
    scale_colour_manual(values = c(ci_match = "#004e92", contrast = "#e05c00"),
                        labels = c(ci_match = "CI matching",
                                   contrast = "Contrast ratio"),
                        name = NULL) +
    labs(x = "Sample size n", y = ylab,
         title = if (!is.null(title)) title else metric) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "#004e92"),
          strip.text = element_text(colour = "white", face = "bold"))
}


# =============================================================================
# 9.  Single-dataset illustration
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("  SECTION 9: Single-dataset illustration\n")
cat(strrep("=", 70), "\n\n", sep = "")

set.seed(CONFIG$demo_seed)
df_demo <- sim_cox_data(n         = CONFIG$demo_n,
                        beta      = CONFIG$demo_beta,
                        cens_rate = CONFIG$demo_cens_rate)
cat(sprintf("n = %d   Events: %d (%.1f%%)   Censored: %d (%.1f%%)\n",
            nrow(df_demo),
            sum(df_demo$status), 100 * mean(df_demo$status),
            sum(!df_demo$status), 100 * mean(df_demo$status == 0)))

out_demo <- cox_iterative_method(
  df            = df_demo,
  j_name        = "x1",
  B             = CONFIG$demo_B,
  m             = CONFIG$demo_m,
  lambda_method = "both",
  R_contrast    = CONFIG$R_contrast,
  target_mass   = CONFIG$target_mass,
  seed          = CONFIG$demo_boot_seed,
  verbose       = TRUE
)

cat("\n--- Interval summaries (illustration dataset) ---\n")
cat(sprintf("  True beta (x1): %.4f  (HR = %.2f)\n",
            CONFIG$demo_beta[1], exp(CONFIG$demo_beta[1])))
for (nm in names(out_demo$results)) {
  r <- out_demo$results[[nm]]
  cat(sprintf("  %-14s  lambda=%7.3f  HDI=[%6.4f, %6.4f]  QI=[%6.4f, %6.4f]\n",
              nm, r$lambda,
              r$hdi["lo"], r$hdi["hi"],
              r$qi["lo"],  r$qi["hi"]))
}
cat(sprintf("  %-14s             Wald CI=[%6.4f, %6.4f]\n",
            "Wald CI", out_demo$wald_ci["lo"], out_demo$wald_ci["hi"]))

# Save illustration plots
fig_pbar <- plot_pbar_band(out_demo)
fig_pi   <- plot_pi_overlay(out_demo)
fig_demo <- fig_pbar + fig_pi +
  plot_annotation(title = "Illustration: n = 300, HR = 1.5, ~20% censoring",
                  theme = theme(plot.title = element_text(face = "bold")))
ggsave(out_path("fig_illustration.pdf"), fig_demo,
       width = 12, height = 5, device = "pdf")
cat("\nIllustration figure saved to", out_path("fig_illustration.pdf"), "\n")


# =============================================================================
# 10.  Simulation study
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("  SECTION 10: Simulation study\n")
cat(strrep("=", 70), "\n\n", sep = "")

# -----------------------------------------------------------------------------
# 10a. Single-replicate runner
# -----------------------------------------------------------------------------
run_one_rep <- function(n, beta_true, cens_rate,
                        B = 500, m = 500,
                        j_name = "x1", seed = NULL) {
  df <- sim_cox_data(n         = n,
                     beta      = c(beta_true, log(0.8), 0.3),
                     cens_rate = cens_rate,
                     seed      = seed)
  cens_pct_actual <- 100 * mean(df$status == 0)

  out <- tryCatch(
    cox_iterative_method(df, j_name = j_name,
                         B = B, m = m,
                         lambda_method = "both",
                         seed = seed,
                         verbose = FALSE),
    error = function(e) NULL
  )
  if (is.null(out)) return(NULL)

  rows <- lapply(names(out$results), function(nm) {
    r <- out$results[[nm]]

    inside_wald <- (out$grid >= out$wald_ci["lo"]) &
                   (out$grid <= out$wald_ci["hi"])

    data.frame(
      method          = nm,
      n               = n,
      beta_true       = beta_true,
      HR_true         = exp(beta_true),
      cens_rate       = cens_rate,
      cens_pct_actual = cens_pct_actual,
      lambda          = r$lambda,
      hdi_lo          = unname(r$hdi["lo"]),
      hdi_hi          = unname(r$hdi["hi"]),
      hdi_covers      = unname((r$hdi["lo"] <= beta_true) &
                                (r$hdi["hi"] >= beta_true)),
      hdi_width       = unname(r$hdi["hi"] - r$hdi["lo"]),
      qi_lo           = unname(r$qi["lo"]),
      qi_hi           = unname(r$qi["hi"]),
      qi_covers       = unname((r$qi["lo"] <= beta_true) &
                                (r$qi["hi"] >= beta_true)),
      qi_width        = unname(r$qi["hi"] - r$qi["lo"]),
      wald_lo         = unname(out$wald_ci["lo"]),
      wald_hi         = unname(out$wald_ci["hi"]),
      wald_covers     = unname((out$wald_ci["lo"] <= beta_true) &
                                (out$wald_ci["hi"] >= beta_true)),
      wald_width      = unname(out$wald_ci["hi"] - out$wald_ci["lo"]),
      mass_in_wald_ci = sum(r$pi[inside_wald]),
      bias_mode       = out$grid[which.max(r$pi)] - beta_true,
      bias_mean       = sum(out$grid * r$pi) - beta_true,
      bhat            = out$bhat,
      D               = out$D,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


# -----------------------------------------------------------------------------
# 10b. Full runner
# -----------------------------------------------------------------------------
run_simulation <- function(n_vec       = c(200, 500, 1000),
                           beta_vec    = c(log(1.1), log(1.5), log(2.0)),
                           cens_rates  = c(0.02, 0.08),
                           R_reps      = 100,
                           B           = 500,
                           m           = 500,
                           base_seed   = 2025,
                           verbose_sim = TRUE) {
  conditions <- expand.grid(n         = n_vec,
                            beta_true = beta_vec,
                            cens_rate = cens_rates,
                            stringsAsFactors = FALSE)
  all_results <- vector("list", nrow(conditions) * R_reps)
  idx <- 1L

  for (ci in seq_len(nrow(conditions))) {
    cond <- conditions[ci, ]
    if (verbose_sim)
      cat(sprintf("\nCondition %d/%d:  n=%d  HR=%.2f  cens_rate=%.3f\n",
                  ci, nrow(conditions),
                  cond$n, exp(cond$beta_true), cond$cens_rate))
    if (verbose_sim) pb2 <- txtProgressBar(min = 0, max = R_reps, style = 3)

    for (r in seq_len(R_reps)) {
      seed_r <- base_seed + ci * 1000 + r
      res <- run_one_rep(n         = cond$n,
                         beta_true = cond$beta_true,
                         cens_rate = cond$cens_rate,
                         B = B, m = m, seed = seed_r)
      all_results[[idx]] <- res
      idx <- idx + 1L
      if (verbose_sim) setTxtProgressBar(pb2, r)
    }
    if (verbose_sim) close(pb2)
  }

  df_raw <- do.call(rbind, Filter(Negate(is.null), all_results))
  # Bin actual censoring % to nearest 10 for display
  df_raw$cens_pct <- round(df_raw$cens_pct_actual / 10) * 10
  df_raw
}

cat("Running simulation — this may take several minutes ...\n")
cat(sprintf("(Design: n = {%s}, HR = {%s}, R_reps = %d, B = %d)\n\n",
            paste(CONFIG$sim_n_vec, collapse = ", "),
            paste(round(exp(CONFIG$sim_beta_vec), 2), collapse = ", "),
            CONFIG$sim_R_reps, CONFIG$sim_B))

sim_df <- run_simulation(
  n_vec      = CONFIG$sim_n_vec,
  beta_vec   = CONFIG$sim_beta_vec,
  cens_rates = CONFIG$sim_cens_rates,
  R_reps     = CONFIG$sim_R_reps,
  B          = CONFIG$sim_B,
  m          = CONFIG$sim_m,
  base_seed  = CONFIG$sim_base_seed
)
saveRDS(sim_df, out_path("sim_results.rds"))
cat(sprintf("\nSimulation complete: %d rows saved to %s\n",
            nrow(sim_df), out_path("sim_results.rds")))


# -----------------------------------------------------------------------------
# 10c. Figures
# -----------------------------------------------------------------------------
cat("\n--- Building simulation figures ---\n")

df_hr15 <- sim_df |> filter(HR_true > 1.4 & HR_true < 1.6)

fig_coverage <- plot_sim_summary(df_hr15, "hdi_covers",
  ylab  = "HDI 95% coverage",
  title = "Coverage: HDI of pi(beta) vs. n  (HR = 1.5)") +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey50")

fig_width <- plot_sim_summary(df_hr15, "hdi_width",
  ylab  = "HDI width",
  title = "Interval width: HDI of pi(beta) vs. n  (HR = 1.5)") +
  stat_summary(
    data = df_hr15 |> filter(method == "ci_match"),
    aes(x = factor(n), y = wald_width, group = 1),
    fun = mean, geom = "line", linetype = "dashed",
    colour = "grey30", linewidth = 0.7, inherit.aes = FALSE
  )

fig_mass <- plot_sim_summary(sim_df, "mass_in_wald_ci",
  ylab  = "Mass of pi(beta) within Wald CI",
  title = "Calibration accuracy: mass within Wald CI") +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey50")

fig_lambda <- plot_sim_summary(df_hr15, "lambda",
  ylab  = "lambda",
  title = "Concentration coefficient lambda vs. n  (HR = 1.5)")

sim_panel <- (fig_coverage | fig_width) / (fig_mass | fig_lambda) +
  plot_annotation(title = "Simulation study summary",
                  theme = theme(plot.title = element_text(face = "bold")))

ggsave(out_path("fig_simulation_panel.pdf"), sim_panel,
       width = 14, height = 10, device = "pdf")
cat("Simulation panel saved to", out_path("fig_simulation_panel.pdf"), "\n")


# =============================================================================
# 11.  Validity checks & summary diagnostics
#      Each check prints PASS / WARN / FAIL.
#      A final tally is printed at the end.
# =============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("  SECTION 11: Validity checks\n")
cat(strrep("=", 70), "\n\n", sep = "")

results_log <- list()   # collects (label, status, detail) for each check

chk <- function(label, condition, detail = "", warn_only = FALSE) {
  status <- if (isTRUE(condition)) "PASS"
             else if (warn_only)   "WARN"
             else                  "FAIL"
  results_log[[length(results_log) + 1]] <<- list(
    label = label, status = status, detail = detail)
  symbol <- c(PASS = "[ OK ]", WARN = "[WARN]", FAIL = "[FAIL]")[status]
  cat(sprintf("  %s  %s\n         %s\n", symbol, label, detail))
}


# ── A. Illustration-dataset checks ───────────────────────────────────────────
cat("── A. Illustration dataset ──────────────────────────────────────────\n")

# A1. p̄ is in [0,1] everywhere
chk("A1: p̄ ∈ [0,1]",
    all(out_demo$pbar >= 0 & out_demo$pbar <= 1),
    sprintf("range = [%.4f, %.4f]",
            min(out_demo$pbar), max(out_demo$pbar)))

# A2. p̄ achieves maximum near β̂
max_at <- out_demo$grid[which.max(out_demo$pbar)]
chk("A2: p̄ maximised near beta_hat",
    abs(max_at - out_demo$bhat) < 3 * out_demo$sehat,
    sprintf("argmax(p̄) = %.4f,  beta_hat = %.4f,  |diff| = %.4f",
            max_at, out_demo$bhat, abs(max_at - out_demo$bhat)))

# A3. Valid-count fraction: at least 80 % of bootstrap reps converged
min_frac <- min(out_demo$valid_counts) / CONFIG$demo_B
chk("A3: Bootstrap convergence rate >= 80 %",
    min_frac >= 0.80,
    sprintf("min valid fraction = %.1f %%", 100 * min_frac))

# A4. π sums to 1 for each method
for (nm in names(out_demo$results)) {
  s <- sum(out_demo$results[[nm]]$pi)
  chk(sprintf("A4: pi sums to 1 (%s)", nm),
      abs(s - 1) < 1e-8,
      sprintf("sum = %.10f", s))
}

# A5. π is non-negative everywhere
for (nm in names(out_demo$results)) {
  neg_count <- sum(out_demo$results[[nm]]$pi < 0)
  chk(sprintf("A5: pi >= 0 (%s)", nm),
      neg_count == 0,
      sprintf("negative values = %d", neg_count))
}

# A6. π achieves mode near β̂
for (nm in names(out_demo$results)) {
  mode_at <- out_demo$grid[which.max(out_demo$results[[nm]]$pi)]
  chk(sprintf("A6: mode of pi near beta_hat (%s)", nm),
      abs(mode_at - out_demo$bhat) < 3 * out_demo$sehat,
      sprintf("mode = %.4f,  beta_hat = %.4f", mode_at, out_demo$bhat))
}

# A7. λ is positive and finite
for (nm in names(out_demo$results)) {
  lam <- out_demo$results[[nm]]$lambda
  chk(sprintf("A7: lambda > 0 and finite (%s)", nm),
      is.finite(lam) && lam > 0,
      sprintf("lambda = %.4f", lam))
}

# A8. HDI lower < HDI upper
for (nm in names(out_demo$results)) {
  hdi <- out_demo$results[[nm]]$hdi
  chk(sprintf("A8: HDI lo < hi (%s)", nm),
      hdi["lo"] < hdi["hi"],
      sprintf("lo = %.4f, hi = %.4f", hdi["lo"], hdi["hi"]))
}

# A9. HDI and QI contain β̂ (the MLE should be inside the interval)
for (nm in names(out_demo$results)) {
  r <- out_demo$results[[nm]]
  in_hdi <- out_demo$bhat >= r$hdi["lo"] & out_demo$bhat <= r$hdi["hi"]
  in_qi  <- out_demo$bhat >= r$qi["lo"]  & out_demo$bhat <= r$qi["hi"]
  chk(sprintf("A9a: beta_hat inside HDI (%s)", nm), in_hdi,
      sprintf("beta_hat = %.4f, HDI = [%.4f, %.4f]",
              out_demo$bhat, r$hdi["lo"], r$hdi["hi"]))
  chk(sprintf("A9b: beta_hat inside QI (%s)", nm), in_qi,
      sprintf("beta_hat = %.4f, QI = [%.4f, %.4f]",
              out_demo$bhat, r$qi["lo"], r$qi["hi"]))
}

# A10. CI-matching: mass within Wald CI close to target
if ("ci_match" %in% names(out_demo$results)) {
  mass <- out_demo$results[["ci_match"]]$mass_in_CI
  tol  <- 0.005
  chk("A10: CI-match mass within tolerance (|mass - 0.95| < 0.005)",
      abs(mass - CONFIG$target_mass) < tol,
      sprintf("mass = %.5f  (target = %.2f, tol = %.3f)", mass,
              CONFIG$target_mass, tol))
}

# A11. True beta covered by HDI on illustration dataset (sanity, not guarantee)
for (nm in names(out_demo$results)) {
  r <- out_demo$results[[nm]]
  beta_true_1 <- CONFIG$demo_beta[1]
  covered <- beta_true_1 >= r$hdi["lo"] & beta_true_1 <= r$hdi["hi"]
  chk(sprintf("A11: true beta covered by HDI (%s) [single sample]", nm),
      covered, warn_only = TRUE,
      sprintf("true beta = %.4f, HDI = [%.4f, %.4f]",
              beta_true_1, r$hdi["lo"], r$hdi["hi"]))
}


# ── B. Simulation-study checks ────────────────────────────────────────────────
cat("\n── B. Simulation study ──────────────────────────────────────────────\n")

# B1. Completeness: number of rows equals expected (2 methods × conditions × reps)
n_conditions <- length(CONFIG$sim_n_vec) *
                length(CONFIG$sim_beta_vec) *
                length(CONFIG$sim_cens_rates)
expected_rows <- 2 * n_conditions * CONFIG$sim_R_reps
chk("B1: Expected number of simulation rows",
    nrow(sim_df) == expected_rows,
    sprintf("got %d, expected %d", nrow(sim_df), expected_rows))

# B2. No NA in key columns
key_cols <- c("lambda", "hdi_covers", "hdi_width", "mass_in_wald_ci",
              "bias_mode", "bias_mean")
na_counts <- sapply(sim_df[key_cols], function(x) sum(is.na(x)))
chk("B2: No NA in key result columns",
    all(na_counts == 0),
    paste(names(na_counts[na_counts > 0]),
          na_counts[na_counts > 0], collapse = ", ") |>
      (\(s) if (nchar(s) == 0) "all clean" else s)())

# B3. All lambda values positive and finite
lam_ok <- with(sim_df, sum(!is.finite(lambda) | lambda <= 0))
chk("B3: All lambda values positive and finite",
    lam_ok == 0,
    sprintf("problematic rows = %d (%.1f %%)",
            lam_ok, 100 * lam_ok / nrow(sim_df)))

# B4. All π-mass values in [0,1]
mass_ok <- with(sim_df, sum(mass_in_wald_ci < 0 | mass_in_wald_ci > 1 + 1e-6))
chk("B4: Mass-in-Wald-CI in [0,1]",
    mass_ok == 0,
    sprintf("out-of-range rows = %d", mass_ok))

# B5. HDI coverage for HR=1.5 within ±10 pp of 95 %
cov_summary <- sim_df |>
  filter(HR_true > 1.4 & HR_true < 1.6) |>
  group_by(method) |>
  summarise(cov = mean(hdi_covers, na.rm = TRUE), .groups = "drop")

for (i in seq_len(nrow(cov_summary))) {
  cov <- cov_summary$cov[i]
  chk(sprintf("B5: HDI coverage in [0.85, 1.00] for HR=1.5 (%s)",
              cov_summary$method[i]),
      cov >= 0.85 & cov <= 1.00, warn_only = TRUE,
      sprintf("coverage = %.3f", cov))
}

# B6. CI-matching mass close to target (mean over all conditions)
ci_mass_mean <- sim_df |>
  filter(method == "ci_match") |>
  summarise(m = mean(mass_in_wald_ci, na.rm = TRUE)) |>
  pull(m)
chk("B6: CI-matching average mass within ±0.02 of target",
    abs(ci_mass_mean - CONFIG$target_mass) < 0.02,
    sprintf("mean mass = %.4f (target = %.2f)",
            ci_mass_mean, CONFIG$target_mass))

# B7. Bias of mode shrinks with n (for HR=1.5, ci_match)
bias_by_n <- sim_df |>
  filter(method == "ci_match", HR_true > 1.4 & HR_true < 1.6) |>
  group_by(n) |>
  summarise(abs_bias = mean(abs(bias_mode), na.rm = TRUE), .groups = "drop") |>
  arrange(n)
chk("B7: |Bias(mode)| weakly decreasing with n (HR=1.5, CI match)",
    all(diff(bias_by_n$abs_bias) <= 0.005),   # allow tiny increases
    paste(sprintf("n=%d: %.4f", bias_by_n$n, bias_by_n$abs_bias),
          collapse = ";  "))

# B8. HDI width shrinks with n (for HR=1.5)
width_by_n <- sim_df |>
  filter(HR_true > 1.4 & HR_true < 1.6) |>
  group_by(method, n) |>
  summarise(w = mean(hdi_width, na.rm = TRUE), .groups = "drop") |>
  arrange(method, n)
for (mth in unique(width_by_n$method)) {
  w <- width_by_n |> filter(method == mth)
  chk(sprintf("B8: HDI width decreasing with n (%s)", mth),
      all(diff(w$w) < 0),
      paste(sprintf("n=%d:%.3f", w$n, w$w), collapse = " "))
}

# B9. Wald CI coverage close to 95 % (sanity check on data generator)
wald_cov <- sim_df |>
  filter(method == "ci_match") |>
  group_by(n) |>
  summarise(wc = mean(wald_covers, na.rm = TRUE), .groups = "drop")
for (i in seq_len(nrow(wald_cov))) {
  chk(sprintf("B9: Wald CI coverage in [0.88, 0.99] (n=%d)", wald_cov$n[i]),
      wald_cov$wc[i] >= 0.88 & wald_cov$wc[i] <= 0.99, warn_only = TRUE,
      sprintf("Wald coverage = %.3f", wald_cov$wc[i]))
}

# B10. lambda increases (weakly) with n for CI-matching
lam_by_n <- sim_df |>
  filter(method == "ci_match", HR_true > 1.4 & HR_true < 1.6) |>
  group_by(n) |>
  summarise(lam = mean(lambda, na.rm = TRUE), .groups = "drop") |>
  arrange(n)
chk("B10: Mean lambda weakly increasing with n (CI match, HR=1.5)",
    all(diff(lam_by_n$lam) >= -5),   # allow small Monte Carlo noise
    paste(sprintf("n=%d:%.1f", lam_by_n$n, lam_by_n$lam), collapse = " "))


# ── C. Final tally ────────────────────────────────────────────────────────────
cat("\n", strrep("─", 70), "\n", sep = "")
statuses <- sapply(results_log, `[[`, "status")
n_pass <- sum(statuses == "PASS")
n_warn <- sum(statuses == "WARN")
n_fail <- sum(statuses == "FAIL")
total  <- length(statuses)

cat(sprintf("  Validity check tally:  %d PASS  |  %d WARN  |  %d FAIL  (of %d)\n",
            n_pass, n_warn, n_fail, total))

if (n_fail > 0) {
  cat("\n  Failed checks:\n")
  for (r in results_log[statuses == "FAIL"])
    cat(sprintf("    [FAIL] %s\n           %s\n", r$label, r$detail))
}
if (n_warn > 0) {
  cat("\n  Warnings (expected on small samples or single runs):\n")
  for (r in results_log[statuses == "WARN"])
    cat(sprintf("    [WARN] %s\n           %s\n", r$label, r$detail))
}
cat(strrep("─", 70), "\n\n", sep = "")

# ── D. Printed summary tables ─────────────────────────────────────────────────
cat("── D. R summary output ──────────────────────────────────────────────\n\n")

cat("HDI coverage and width by method, n, and censoring level:\n")
print(as.data.frame(tbl_hdi), digits = 3, row.names = FALSE)

cat("\nBias of mode and mean of pi by condition:\n")
print(as.data.frame(tbl_bias), digits = 4, row.names = FALSE)

cat("\nConcentration coefficient lambda by method, n, and censoring level:\n")
print(as.data.frame(tbl_lambda), digits = 3, row.names = FALSE)

cat("\nWald CI coverage (sanity check on data generator):\n")
print(as.data.frame(wald_cov), digits = 3, row.names = FALSE)

cat("\n", strrep("=", 70), "\n", sep = "")
cat("  All done.  Outputs written to: ./", OUT_DIR, "/\n", sep = "")
cat(strrep("=", 70), "\n\n", sep = "")
