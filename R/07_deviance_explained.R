# 07_deviance_explained.R
#
# Justifies the choice of k = 5 latent dimensions used to fit the real
# armed-conflict data (Section 3 of the paper; reproduces Figure 2). Fits
# the model on the full data (no held-out cells, no imputation) for
# k = 0..13 (the maximum possible given p = 14) and computes the
# percentage of the null-model deviance explained by each additional
# dimension -- a scree-plot-style analysis analogous to proportion of
# variance explained in classical PCA.
#
# Takes several minutes (13 full-data PDLB fits at n = 7,165, p = 14).

suppressMessages({
  library(BiplotML)
  library(dplyr)
})

conflict <- read.csv("data/conflict_arm_data.csv", stringsAsFactors = FALSE)
x <- as.matrix(conflict %>% select(-id, -armed_group))
n <- nrow(x); p <- ncol(x)

neg_log_lik <- function(x, theta) {
  P <- plogis(theta)
  eps <- 1e-10
  P <- pmin(pmax(P, eps), 1 - eps)
  -sum(x * log(P) + (1 - x) * log(1 - P))
}

# k = 0: null (independence) model, theta_ij = mu_j = logit(colMeans(x))
mu0 <- colMeans(x)
theta0 <- matrix(qlogis(pmin(pmax(mu0, 1e-6), 1 - 1e-6)), n, p, byrow = TRUE)
dev0 <- 2 * neg_log_lik(x, theta0)

results <- data.frame(k = 0, deviance = dev0, pct_explained = 0)

for (k in 1:13) {
  fit <- LogBip(x, k = k, method = "PDLB", maxit = 1000, plot = FALSE)
  theta_hat <- fitted_LB(fit, type = "link")
  dev_k <- 2 * neg_log_lik(x, theta_hat)
  pct <- 100 * (dev0 - dev_k) / dev0
  cat(sprintf("k=%2d | deviance=%.1f | deviance explained=%.2f%%\n", k, dev_k, pct))
  results <- rbind(results, data.frame(k = k, deviance = dev_k, pct_explained = pct))
}

results$marginal_gain <- c(NA, diff(results$pct_explained))
print(results)

dir.create("output", showWarnings = FALSE)
write.csv(results, "output/deviance_explained.csv", row.names = FALSE)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  p_plot <- ggplot(results, aes(x = k, y = pct_explained)) +
    geom_line(color = "#0E185F") +
    geom_point(size = 2.5, color = "#0E185F") +
    geom_vline(xintercept = 5, linetype = 2, color = "grey40") +
    scale_x_continuous(breaks = 0:13) +
    labs(x = "Number of dimensions (k)", y = "Deviance explained (%)") +
    theme_bw(base_size = 13)
  ggsave("output/deviance_explained.pdf", p_plot, width = 6, height = 4.5)
  ggsave("output/deviance_explained.eps", p_plot, width = 6, height = 4.5, device = cairo_ps)
  cat("Saved output/deviance_explained.{csv,pdf,eps}\n")
}
