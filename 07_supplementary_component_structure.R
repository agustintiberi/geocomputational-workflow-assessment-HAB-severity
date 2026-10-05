# =====================================================================
# 07_supplementary_component_structure.R
# ---------------------------------------------------------------------
# Dimensional structure of the Bloom Severity (BS) components:
# pairwise Spearman correlations among BB, RBA and RBF (Table S3).
# Script 7 of the workflow.
#
# Associated publication:
#   Tiberi, A. E., A. A. Drozd, H. R. Zerda & P. de Tezanos Pinto, 2026.
#   An open geocomputational workflow for enabling long-term and
#   large-scale assessment of harmful algal bloom severity.
#   Hydrobiologia. https://doi.org/10.1007/s10750-026-06415-5
#
# Copyright (C) 2026 Agustín Eduardo Tiberi
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version. See the LICENSE file for details.
# ---------------------------------------------------------------------
# Input : MODIS_monthly_BB_RBA_RBF.csv (from 02_MODIS_daily_to_monthly.R)
#
# Tests whether the three components of the temporal BS index are three
# independent dimensions or fewer. Structure reported in the paper:
#   RBA and RBF co-vary  -> a single spatio-temporal ACTIVITY dimension
#   BB is independent    -> a separate INTENSITY dimension
# Under this structure, equal (1:1:1) weighting represents activity twice
# and intensity once, implicitly weighting activity over intensity by
# about 2:1.
#
# RBF is derived from RBA by construction (it counts the dates on which
# RBA exceeds a threshold), so their association is both conceptual and
# constructional; this reinforces rather than weakens the single-
# dimension interpretation.
#
# Outputs:
#   TableS3_component_correlations.csv   Table S3
#   component_correlations.png           scatterplots of the three pairs
# =====================================================================

# ------------------------------ CONFIG -------------------------------
in_csv  <- "MODIS_monthly_BB_RBA_RBF.csv"
out_dir <- "."

# ---------------------------- LOAD -----------------------------------
d <- read.csv(in_csv, stringsAsFactors = FALSE)
d$ym <- as.Date(d$ym)
d <- d[is.finite(d$BB) & is.finite(d$RBA) & is.finite(d$RBF), ]

cat(sprintf("Months with all three components defined: %d (%s to %s)\n",
            nrow(d), min(d$ym), max(d$ym)))

# ------------------------ CORRELATIONS -------------------------------
pairs <- list(c("RBA", "RBF"), c("BB", "RBA"), c("BB", "RBF"))

res <- do.call(rbind, lapply(pairs, function(p) {
  ct <- cor.test(d[[p[1]]], d[[p[2]]], method = "spearman", exact = FALSE)
  data.frame(pair    = paste(p[1], p[2], sep = "-"),
             rho     = round(unname(ct$estimate), 3),
             p_value = signif(ct$p.value, 3),
             n       = nrow(d))
}))

cat("\n=== Spearman correlations among BS components (Table S3) ===\n")
print(res, row.names = FALSE)

cat("\nReading: a high RBA-RBF value with low BB-RBA and BB-RBF values\n")
cat("indicates two dimensions (activity and intensity), not three.\n")

# --------------------------- FIGURE ----------------------------------
png(file.path(out_dir, "component_correlations.png"),
    width = 9, height = 3.2, units = "in", res = 300, type = "cairo")
op <- par(mfrow = c(1, 3), mar = c(4, 4, 2.4, 1), mgp = c(2.2, 0.7, 0))

sp <- function(x, y, xl, yl) {
  ct <- cor.test(x, y, method = "spearman", exact = FALSE)
  plot(x, y, pch = 19, col = "#2c7fb8aa", cex = 0.8,
       xlab = xl, ylab = yl,
       main = sprintf("rho = %.2f", unname(ct$estimate)))
}
with(d, sp(RBA, RBF, "RBA (area)",     "RBF (frequency)"))
with(d, sp(BB,  RBA, "BB (intensity)", "RBA (area)"))
with(d, sp(BB,  RBF, "BB (intensity)", "RBF (frequency)"))

par(op); dev.off()

write.csv(res, file.path(out_dir, "TableS3_component_correlations.csv"),
          row.names = FALSE)
cat("\n>>> Written: TableS3_component_correlations.csv, component_correlations.png\n")
