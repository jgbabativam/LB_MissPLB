# 01_fit_model.R
#
# Fits the logistic biplot model to the armed-conflict data using the
# PDLB (Projection-based Data Logistic Biplot) algorithm 
#
# This reproduces the model fitted in Section 3 ("Application with Real
# Data") of the paper. 
#
# Requires: install.packages("BiplotML")

library(BiplotML)
library(dplyr)

conflict <- read.csv("data/conflict_arm_data.csv", stringsAsFactors = FALSE)

x <- conflict %>% select(-id, -armed_group)

set.seed(12345)
fit <- LogBip(x, k = 5, method = "PDLB", maxit = 1000, plot = FALSE)

dir.create("output", showWarnings = FALSE)
saveRDS(fit, "output/fitted_model.rds")

cat("Model fitted with method =", fit$method, "in", fit$iterations, "iterations.\n")
cat("Saved to output/fitted_model.rds\n")
