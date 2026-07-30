# 05_accuracy_comparison.R
#
# Statistical accuracy of PDLB versus MM (Section 2.6 of the paper).
#
# Section 2.5 shows PDLB is faster because it estimates p*k+p parameters
# instead of (n+p)*k+p. This script checks whether that parameter reduction
# costs statistical accuracy, using BiplotML::simBin(), which returns the
# TRUE probability matrix P used to generate the data. This makes it
# possible to measure recovery of the underlying signal (RMSE against the
# true P), not just fit to the observed, noisy sample (in-sample deviance).
#
# Both algorithms are run at n = 200, 1,000, 5,000, with 3 replicates each
# (fresh simulated data every time), matching p = 14, k = 5, and the
# marginal probability of ones (0.15) used in the benchmark of Section 2.5.
#
# Takes about 10-15 minutes (MM needs up to ~3,000 iterations at these n).

suppressMessages({
  library(BiplotML)
  library(dplyr)
})

auc <- function(prob, label) {
  # Mann-Whitney U based AUC, pooled over all cells
  r <- rank(prob)
  n1 <- sum(label == 1); n0 <- sum(label == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[label == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

deviance_bernoulli <- function(P, X) {
  eps <- 1e-10
  P <- pmin(pmax(P, eps), 1 - eps)
  -2 * sum(X * log(P) + (1 - X) * log(1 - P))
}

p <- 14L
k <- 5L
D <- 0.15
n_grid <- c(200, 1000, 5000)
n_reps <- 3L

run_pair <- function(n, rep_id) {
  set.seed(2026 * 100 + rep_id * 17 + n)
  sim <- simBin(n = n, p = p, k = k, D = D)
  X <- sim$X
  Ptrue <- sim$P

  fit_mm <- LogBip(X, k = k, method = "MM",   maxit = 3000, plot = FALSE)
  fit_pd <- LogBip(X, k = k, method = "PDLB", maxit = 1000, plot = FALSE)

  P_mm <- fitted_LB(fit_mm, type = "response")
  P_pd <- fitted_LB(fit_pd, type = "response")

  data.frame(
    n = n, rep = rep_id,
    dev_mm = deviance_bernoulli(P_mm, X),
    dev_pd = deviance_bernoulli(P_pd, X),
    rmse_true_mm = sqrt(mean((P_mm - Ptrue)^2)),
    rmse_true_pd = sqrt(mean((P_pd - Ptrue)^2)),
    auc_mm = auc(as.vector(P_mm), as.vector(X)),
    auc_pd = auc(as.vector(P_pd), as.vector(X)),
    cor_mm_pd = cor(as.vector(P_mm), as.vector(P_pd)),
    rmse_mm_pd = sqrt(mean((P_mm - P_pd)^2)),
    iterations_mm = fit_mm$iterations,
    iterations_pd = fit_pd$iterations,
    converged_mm = fit_mm$iterations < 3000,
    converged_pd = fit_pd$iterations < 1000,
    max_abs_B_mm = max(abs(as.matrix(fit_mm$Bhat))),
    max_abs_B_pd = max(abs(as.matrix(fit_pd$Bhat)))
  )
}

results <- list()
i <- 1L
for (n in n_grid) {
  for (r in seq_len(n_reps)) {
    cat("n =", n, "rep =", r, "...\n")
    results[[i]] <- run_pair(n, r)
    i <- i + 1L
  }
}

raw <- bind_rows(results)
dir.create("output", showWarnings = FALSE)
write.csv(raw, "output/accuracy_comparison_raw.csv", row.names = FALSE)

summary_tbl <- raw %>%
  group_by(n) %>%
  summarise(
    dev_mm_mean = mean(dev_mm), dev_pd_mean = mean(dev_pd),
    rmse_true_mm_mean = mean(rmse_true_mm), rmse_true_mm_sd = sd(rmse_true_mm),
    rmse_true_pd_mean = mean(rmse_true_pd), rmse_true_pd_sd = sd(rmse_true_pd),
    auc_mm_mean = mean(auc_mm), auc_pd_mean = mean(auc_pd),
    cor_mm_pd_mean = mean(cor_mm_pd),
    rmse_mm_pd_mean = mean(rmse_mm_pd),
    pct_converged_mm = mean(converged_mm) * 100,
    pct_converged_pd = mean(converged_pd) * 100,
    max_abs_B_mm_mean = mean(max_abs_B_mm),
    max_abs_B_pd_mean = mean(max_abs_B_pd),
    .groups = "drop"
  )

print(as.data.frame(raw))
cat("\n\n=== SUMMARY (Table in Section 2.6 of the paper) ===\n")
print(as.data.frame(summary_tbl))
write.csv(summary_tbl, "output/accuracy_comparison_summary.csv", row.names = FALSE)

cat("\nDone. See output/accuracy_comparison_summary.csv (Table 2.6) and\n")
cat("output/accuracy_comparison_raw.csv (all", n_reps, "replicates).\n")
