# =====================================================================
# 02_MODIS_daily_to_monthly.R
# ---------------------------------------------------------------------
# Daily MODIS-FAI metrics -> filtered daily series and monthly series.
# Script 2 of the workflow.
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
# Input : MODIS_FAI_daily_metrics.csv  (exported by 01_MODIS_FAI_daily_metrics.js)
#         columns: date, n_valid, n_bloom, RBA, BB_bloom, BB_allwater
#
# Steps :
#   1. Load and check the daily series
#   2. Apply the minimum-coverage filter (kept here, not in GEE, so that
#      the criterion is explicit and auditable)
#   3. Diagnostics for the two candidate BB definitions
#   4. Aggregate to monthly BB and RBA, and compute monthly RBF
#   5. Write outputs and a diagnostic figure
#
# Outputs:
#   MODIS_daily_filtered.csv       daily series after the coverage filter
#   MODIS_monthly_BB_RBA_RBF.csv   monthly BB, RBA and RBF
#   MODIS_monthly_series.png       diagnostic plot of the monthly series
#
# Base R only (no additional packages); PNG written through Cairo.
# =====================================================================

# ------------------------------ CONFIG -------------------------------
in_csv   <- "MODIS_FAI_daily_metrics.csv"
out_dir  <- "."

# Minimum-coverage filter. A date is usable when the number of valid
# water cells reaches a given fraction of the reference extent.
min_cov_frac  <- 0.65   # fraction of the reference extent required (65%)

# Reference extent: number of 250 m cells in the minimum (permanent)
# reservoir extent polygon (AOI), which is also the maximum number of
# valid cells observable on any date. Set to NA to estimate it from the
# data as the 99th percentile of n_valid (gives the same value here).
n_total_cells <- 2154

# Monthly RBF: a date counts as a significant-bloom date when RBA
# exceeds this value.
rbf_rba_thr <- 0.25

# Minimum number of usable images for a month to be retained.
min_img_per_month <- 3

# BB definition carried forward. "bloom" (median FAI over bloom cells
# only) is the definition used in the paper; "allwater" is kept for
# diagnostic purposes only.
bb_definition <- "bloom"

# --------------------------- 1. LOAD DATA ----------------------------
d <- read.csv(in_csv, stringsAsFactors = FALSE)
d$date <- as.Date(d$date)
for (v in c("n_valid", "n_bloom", "RBA", "BB_bloom", "BB_allwater")) {
  d[[v]] <- suppressWarnings(as.numeric(d[[v]]))
}
d <- d[!is.na(d$date) & is.finite(d$n_valid), ]
d <- d[order(d$date), ]

cat(sprintf("Daily records loaded: %d (%s to %s)\n",
            nrow(d), min(d$date), max(d$date)))

# ----------------------- 2. COVERAGE FILTER --------------------------
if (is.na(n_total_cells)) {
  n_total_cells <- as.numeric(quantile(d$n_valid, 0.99, na.rm = TRUE))
  cat(sprintf("Reference extent estimated from data: %.0f cells (99th pct of n_valid)\n",
              n_total_cells))
} else {
  cat(sprintf("Reference extent supplied: %.0f cells\n", n_total_cells))
}

min_cells <- min_cov_frac * n_total_cells
d$coverage <- d$n_valid / n_total_cells
keep <- d$n_valid >= min_cells

cat(sprintf("Coverage threshold: %.0f cells (%.0f%% of reference extent)\n",
            min_cells, 100 * min_cov_frac))
cat(sprintf("Dates retained: %d of %d (%.1f%%)\n",
            sum(keep), nrow(d), 100 * mean(keep)))

d <- d[keep, ]

# ------------------- 3. BB DEFINITION DIAGNOSTICS --------------------
# BB_bloom   : median FAI over bloom cells only (undefined on bloom-free dates)
# BB_allwater: median FAI over all valid water cells
# Diagnostic only: reproduces the comparison reported in the paper
# (Methods, "Algal bloom temporal variables computation"); it does not
# alter the definition used downstream.
n_nobloom <- sum(!is.finite(d$BB_bloom))
cat(sprintf("\nDates with no bloom cells (BB_bloom undefined): %d (%.1f%%)\n",
            n_nobloom, 100 * n_nobloom / nrow(d)))

both <- is.finite(d$BB_bloom) & is.finite(d$BB_allwater)
if (sum(both) > 10) {
  rho_defs <- cor(d$BB_bloom[both], d$BB_allwater[both], method = "spearman")
  cat(sprintf("Spearman rho, BB_bloom vs BB_allwater: %.3f (n = %d)\n",
              rho_defs, sum(both)))
}
rho_all_rba <- cor(d$BB_allwater, d$RBA, method = "spearman",
                   use = "complete.obs")
cat(sprintf("Spearman rho, BB_allwater vs RBA: %.3f\n", rho_all_rba))

d$BB <- if (bb_definition == "bloom") d$BB_bloom else d$BB_allwater
cat(sprintf("\nBB definition carried forward: %s\n", bb_definition))

# --------------------- 4. MONTHLY AGGREGATION ------------------------
# BB  : median of daily BB over dates on which BB is defined
#       (dates with no bloom cell are treated as missing, not as zero)
# RBA : mean of daily RBA
# RBF : number of dates with RBA > rbf_rba_thr, divided by the number of
#       usable images in that month (normalises for cloud-limited months)
d$ym <- as.Date(cut(d$date, "month"))

months <- sort(unique(d$ym))
monthly <- do.call(rbind, lapply(months, function(m) {
  s <- d[d$ym == m, ]
  n_img <- nrow(s)
  data.frame(
    ym    = m,
    n_img = n_img,
    BB    = if (any(is.finite(s$BB))) median(s$BB[is.finite(s$BB)]) else NA_real_,
    RBA   = mean(s$RBA, na.rm = TRUE),
    RBF   = sum(s$RBA > rbf_rba_thr, na.rm = TRUE) / max(n_img, 1)
  )
}))

n_before <- nrow(monthly)
monthly  <- monthly[monthly$n_img >= min_img_per_month, ]
cat(sprintf("\nMonths retained: %d of %d (>= %d usable images)\n",
            nrow(monthly), n_before, min_img_per_month))
cat(sprintf("Monthly coverage: %s to %s\n", min(monthly$ym), max(monthly$ym)))
cat(sprintf("Months with undefined BB (no bloom cell all month): %d\n",
            sum(!is.finite(monthly$BB))))

# ----------------------- 5. OUTPUTS + FIGURE -------------------------
write.csv(d[, c("date", "n_valid", "n_bloom", "coverage", "RBA",
                "BB_bloom", "BB_allwater")],
          file.path(out_dir, "MODIS_daily_filtered.csv"), row.names = FALSE)
write.csv(monthly, file.path(out_dir, "MODIS_monthly_BB_RBA_RBF.csv"),
          row.names = FALSE)

png(file.path(out_dir, "MODIS_monthly_series.png"),
    width = 9, height = 7, units = "in", res = 300, type = "cairo")
op <- par(mfrow = c(3, 1), mar = c(3.2, 4.2, 1.6, 1), mgp = c(2.3, 0.7, 0))
plot(monthly$ym, monthly$BB,  type = "l", col = "#d95f02", lwd = 1.4,
     xlab = "", ylab = "BB (median FAI)")
plot(monthly$ym, monthly$RBA, type = "l", col = "#1b9e77", lwd = 1.4,
     xlab = "", ylab = "RBA (bloom fraction)")
plot(monthly$ym, monthly$RBF, type = "l", col = "#7570b3", lwd = 1.4,
     xlab = "", ylab = "RBF (significant-bloom date fraction)")
par(op); dev.off()

cat("\n>>> Written: MODIS_daily_filtered.csv, MODIS_monthly_BB_RBA_RBF.csv,",
    "MODIS_monthly_series.png\n")
