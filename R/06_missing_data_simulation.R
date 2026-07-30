# 06_missing_data_simulation.R
#
# Reproduces the missing-data simulation study in Section 2.7 of the paper
# ("Recovery Under Missing Data"). Simulated binary matrices (p = 14, k = 5,
# density 0.15, matching the real data example) are masked under MCAR and
# MAR mechanisms at 0%, 10%, and 30% missingness, fitted with the proposed
# PDLB algorithm, and compared to the known simulation ground truth.
#
#
# Takes a few minutes to run (50 model fits: 5 conditions x 10 replicates).

suppressMessages({
  library(BiplotML)
  library(dplyr)
})

p <- 14L; k <- 5L; n <- 2000L; D <- 0.15
n_reps <- 10L
rates <- c(0.00, 0.10, 0.30)
mechanisms <- c("MCAR", "MAR")

rmse <- function(a, b) sqrt(mean((a - b)^2))

make_mask_mcar <- function(n, p, rate) {
  matrix(rbinom(n * p, 1, rate), n, p) == 1
}

# MAR: probability that a cell is missing depends on the row's total number
# of characteristics (rowSums of the TRUE, fully observed X) -- an auxiliary
# that is itself always observed, so this is a valid MAR mechanism
# (missingness depends on other observed data, not on the value of the cell
# itself). Probability increases linearly with rank so the realized overall
# rate matches `rate` exactly in expectation.
make_mask_mar <- function(X, rate) {
  n <- nrow(X); p <- ncol(X)
  aux <- rowSums(X)
  rnk <- rank(aux, ties.method = "average")
  p_row <- pmin(pmax(2 * rate * (rnk - 0.5) / n, 0), 1)
  mask <- matrix(FALSE, n, p)
  for (i in seq_len(n)) mask[i, ] <- rbinom(p, 1, p_row[i]) == 1
  mask
}

good_fit_vec <- function(x_true, x_hat, mask) {
  xt <- x_true[mask]; xh <- x_hat[mask]
  ones <- sum(xt == 1); zeros <- sum(xt == 0)
  sens <- if (ones > 0) mean(xh[xt == 1] == 1) * 100 else NA
  spec <- if (zeros > 0) mean(xh[xt == 0] == 0) * 100 else NA
  acc  <- mean(xh == xt) * 100
  c(sensitivity = sens, specificity = spec, accuracy = acc)
}

run_one <- function(rep_id, mechanism, rate) {
  set.seed(1000 * rep_id + round(rate * 100) + (mechanism == "MAR") * 7)
  sim <- simBin(n = n, p = p, k = k, D = D)

  mask <- if (rate == 0) {
    matrix(FALSE, n, p)
  } else if (mechanism == "MCAR") {
    make_mask_mcar(n, p, rate)
  } else {
    make_mask_mar(sim$X, rate)
  }

  X_obs <- sim$X
  X_obs[mask] <- NA

  t0 <- Sys.time()
  fit <- tryCatch(
    LogBip(X_obs, k = k, method = "PDLB", maxit = 1000, plot = FALSE),
    error = function(e) NULL
  )
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")
  if (is.null(fit)) return(NULL)

  Theta_hat <- fitted_LB(fit, type = "link")
  P_hat     <- fitted_LB(fit, type = "response")
  mu_hat    <- fit$Bhat$bb0

  imp <- if (rate == 0) c(sensitivity = NA, specificity = NA, accuracy = NA) else
    good_fit_vec(sim$X, fit$impute_x, mask)

  data.frame(
    rep = rep_id, mechanism = mechanism, rate = rate,
    realized_missing_pct = 100 * mean(mask),
    imp_sensitivity = imp["sensitivity"], imp_specificity = imp["specificity"],
    imp_accuracy = imp["accuracy"],
    rmse_theta = rmse(Theta_hat, sim$Theta),
    rmse_P = rmse(P_hat, sim$P),
    rmse_mu = rmse(mu_hat, sim$mu),
    iterations = fit$iterations, elapsed_sec = elapsed
  )
}

conditions <- data.frame(rate = 0, mechanism = "MCAR")  # 0% baseline (mechanism irrelevant)
for (r in rates[rates > 0]) for (m in mechanisms) conditions <- rbind(conditions, data.frame(rate = r, mechanism = m))

results <- list(); idx <- 1L
for (cc in seq_len(nrow(conditions))) {
  for (rep_id in seq_len(n_reps)) {
    cat(sprintf("rate=%.0f%% mechanism=%s rep=%d/%d\n",
                100 * conditions$rate[cc], conditions$mechanism[cc], rep_id, n_reps))
    out <- run_one(rep_id, conditions$mechanism[cc], conditions$rate[cc])
    if (!is.null(out)) { results[[idx]] <- out; idx <- idx + 1L }
  }
}

full <- bind_rows(results)
dir.create("output", showWarnings = FALSE)
write.csv(full, "output/missing_data_simulation_raw.csv", row.names = FALSE)

summary_tbl <- full %>%
  group_by(rate, mechanism) %>%
  summarise(
    n_reps = n(),
    realized_missing_pct = mean(realized_missing_pct),
    imp_sensitivity = mean(imp_sensitivity, na.rm = TRUE),
    imp_specificity = mean(imp_specificity, na.rm = TRUE),
    imp_accuracy = mean(imp_accuracy, na.rm = TRUE),
    rmse_theta_mean = mean(rmse_theta), rmse_theta_sd = sd(rmse_theta),
    rmse_P_mean = mean(rmse_P), rmse_P_sd = sd(rmse_P),
    rmse_mu_mean = mean(rmse_mu), rmse_mu_sd = sd(rmse_mu),
    iterations_mean = mean(iterations),
    elapsed_mean = mean(elapsed_sec),
    .groups = "drop"
  ) %>%
  arrange(rate, mechanism)

print(as.data.frame(summary_tbl))
write.csv(summary_tbl, "output/missing_data_simulation_summary.csv", row.names = FALSE)
cat("\nSaved output/missing_data_simulation_raw.csv and output/missing_data_simulation_summary.csv\n")
