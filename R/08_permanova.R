# 08_permanova.R
#
# Tests whether the responsible-armed-group label explains variation in the
# fitted row markers (Section 3 of the paper), using a permutation-based
# multivariate analysis of variance (PERMANOVA; Anderson 2001) on the
# k = 5 dimensional row markers from the model fitted in 01_fit_model.R.
#
# Requires the "vegan" package: install.packages("vegan")
#
# Takes a few minutes (100 permutations on n = 7,165 rows).

suppressMessages({
  library(BiplotML)
  library(dplyr)
  library(vegan)
})

conflict <- read.csv("data/conflict_arm_data.csv", stringsAsFactors = FALSE)
x <- conflict %>% select(-id, -armed_group)

fit <- readRDS("output/fitted_model.rds")
A <- as.matrix(fit$Ahat)
group <- factor(conflict$armed_group)

cat("Group sizes:\n")
print(table(group))

set.seed(2026)
res <- adonis2(A ~ group, method = "euclidean", permutations = 100)
print(res)

dir.create("output", showWarnings = FALSE)
sink("output/permanova_result.txt")
print(table(group))
cat("\n")
print(res)
sink()
cat("\nSaved output/permanova_result.txt\n")
