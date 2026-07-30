# 03_figure1_biplot.R
#
# Reproduces Figure 1 of the paper: the logistic biplot for the
# armed-conflict data, with row markers coloured by the responsible
# armed group (a variable external to the model, used only for display).
#
# Run 01_fit_model.R first.

library(BiplotML)
library(ggplot2)

conflict <- read.csv("data/conflict_arm_data.csv", stringsAsFactors = FALSE)
fit <- readRDS("output/fitted_model.rds")

p <- plotBLB(fit, dim = c(1, 2), xylim = c(-90, 90), escala = 70,
             col.ind = conflict$armed_group) +
  theme(legend.position = "top", legend.title = element_blank()) +
  labs(title = "", subtitle = "")

dir.create("output", showWarnings = FALSE)
ggsave("output/figure1_biplot.pdf", p, width = 8, height = 8)
ggsave("output/figure1_biplot.eps", p, width = 8, height = 8, device = cairo_ps)

print(p)
