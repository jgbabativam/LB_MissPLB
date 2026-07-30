# 04_benchmark.R
#
# Computational benchmark: wall-clock fitting time and number of estimated
# parameters, comparing the algorithm proposed in this paper ("PDLB",
# block coordinate descent with data projection) against the previous
# coordinate descent MM algorithm ("MM", Babativa-Marquez & Vicente-
# Villardon 2021), as a function of the sample size n.
#
# Data are simulated with BiplotML::simBin().
#
# NOTE: this script takes roughly 20-30 minutes to run in full (MM needs
# hundreds to thousands of iterations to converge, especially at large n).

library(BiplotML)
library(dplyr)
library(ggplot2)

p <- 14L   # number of variables, matching the real conflict data example
k <- 5L    # number of latent dimensions, matching the real data example
D <- 0.15  # marginal probability of a one, approximately matching the real data (~16.5%)

n_grid <- c(200, 500, 1000, 2000, 5000, 10000, 20000, 50000)
n_reps <- 3L

run_one <- function(n, method, maxit, rep_id) {
  set.seed(2026 * 100 + rep_id)
  sim <- simBin(n = n, p = p, k = k, D = D)
  t <- system.time(
    fit <- LogBip(sim$X, k = k, method = method, maxit = maxit, plot = FALSE)
  )
  data.frame(
    n = n, method = method, rep = rep_id,
    elapsed_sec = as.numeric(t["elapsed"]),
    iterations  = fit$iterations,
    converged   = fit$iterations < maxit
  )
}

results <- list()
i <- 1L

for (n in n_grid) {
  for (r in seq_len(n_reps)) {
    cat("n =", n, "rep =", r, "- fitting MM...\n")
    results[[i]] <- run_one(n, "MM", maxit = 3000, rep_id = r); i <- i + 1L
    cat("n =", n, "rep =", r, "- fitting PDLB...\n")
    results[[i]] <- run_one(n, "PDLB", maxit = 1000, rep_id = r); i <- i + 1L
  }
}

raw <- bind_rows(results) %>%
  mutate(
    n_parameters = ifelse(method == "MM", (n + p) * k + p, p * k + p)
  )

dir.create("output", showWarnings = FALSE)
write.csv(raw, "output/benchmark_results_raw.csv", row.names = FALSE)

bench <- raw %>%
  group_by(n, method, n_parameters) %>%
  summarise(
    elapsed_mean = mean(elapsed_sec), elapsed_sd = sd(elapsed_sec),
    iterations_mean = mean(iterations), iterations_sd = sd(iterations),
    n_reps = n(),
    .groups = "drop"
  ) %>%
  arrange(n, method)

print(as.data.frame(bench))
write.csv(bench, "output/benchmark_results.csv", row.names = FALSE)

p_time <- ggplot(bench, aes(x = n, y = elapsed_mean, color = method, shape = method)) +
  geom_ribbon(aes(ymin = pmax(elapsed_mean - elapsed_sd, 1e-3), ymax = elapsed_mean + elapsed_sd,
                   fill = method), alpha = 0.15, color = NA) +
  geom_point(size = 2.5) +
  geom_line() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Sample size (n, log scale)", y = "Fitting time in seconds (log scale)",
       color = "Method", shape = "Method", fill = "Method") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top")

ggsave("output/benchmark_time.pdf", p_time, width = 6, height = 5)

## --- Figure: number of parameters vs n ---------------------------------
param_curve <- expand.grid(n = seq(min(n_grid), max(n_grid), length.out = 200),
                            method = c("MM", "PDLB")) %>%
  mutate(n_parameters = ifelse(method == "MM", (n + p) * k + p, p * k + p))

p_param <- ggplot(param_curve, aes(x = n, y = n_parameters, color = method)) +
  geom_line(linewidth = 1) +
  labs(x = "Sample size (n)", y = "Number of estimated parameters",
       color = "Method") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top")

ggsave("output/benchmark_parameters.pdf", p_param, width = 6, height = 5)

cat("\nDone. See output/benchmark_results.csv (summary) and\n")
cat("output/benchmark_results_raw.csv (all", n_reps, "replicates).\n")
