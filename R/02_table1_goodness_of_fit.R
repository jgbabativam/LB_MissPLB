# 02_table1_goodness_of_fit.R
#
# Reproduces Table 1 of the paper: sensitivity, specificity, and overall
# correct-classification rate ("Global") for each of the 14 violence
# types.
#
# Run 01_fit_model.R first.

library(BiplotML)
library(dplyr)

conflict <- read.csv("data/conflict_arm_data.csv", stringsAsFactors = FALSE)
x <- conflict %>% select(-id, -armed_group)

fit <- readRDS("output/fitted_model.rds")

Pi  <- fitted_LB(fit, type = "response")
thr <- BiplotML:::thresholds(x = x, Pi)

good_fit <- function(x, xhat) {
  pcc   <- ifelse((x == 1 & xhat == 1) | (x == 0 & xhat == 0), 1, 0)
  ones  <- apply(x, 2, sum)
  zeros <- nrow(x) - ones
  data.frame(
    Sensitivity = round(100 * apply((xhat == 1) & (x == 1), 2, sum) / ones, 1),
    Specificity = round(100 * apply((xhat == 0) & (x == 0), 2, sum) / zeros, 1),
    Global      = round(100 * colSums(pcc) / nrow(pcc), 1)
  )
}

table1 <- good_fit(x, thr$pred)
table1 <- table1[order(-table1$Global), ]

print(table1)

dir.create("output", showWarnings = FALSE)
write.csv(table1, "output/table1_goodness_of_fit.csv")
