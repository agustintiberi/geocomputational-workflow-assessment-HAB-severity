# =====================================================================
# 08_fig_temporal_severity.R
# ---------------------------------------------------------------------
# Monthly temporal Bloom Severity (BS): normalisation, trend and change
# point analysis, and Figure 7 of the paper. Script 8 of the workflow.
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
# BS = BB_n + RBA_n + RBF_n  (equal 1:1:1 weighting, range 0-3; Eq. 11a)
# Each component is winsorised at the 1st/99th percentiles and rescaled
# to 0-1 before summing (Eqs. 10a and 10b), so that a few extreme months
# do not dominate the composite.
#
# Months missing from the input (fewer than three usable images; one
# month, April 2017, in the published record) are handled as follows:
# short gaps are linearly interpolated ONLY to complete the regular
# monthly sequence over which the normalisation range is set; the
# interpolated months are then excluded from every statistic, summary
# and output value. All tests are therefore computed on the 309 months
# with observed data.
#
# Statistics reported in the paper (Results, "Long-term severity"):
#   - modified Mann-Kendall (Hamed & Rao, 1998) on BS and on each
#     normalised component
#   - Pettitt change point on BS
#   - mean and median BS, and mean of each component, before and after
#     the change point
#
# A 12-month centred moving mean is overlaid in the figure because the
# raw monthly series is dominated by the seasonal cycle.
#
# Outputs:
#   fig_temporal_severity.png      Figure 7
#   MODIS_monthly_severity.csv     monthly components, BS and moving means
#
# Requires: trend (Pettitt), modifiedmk (Hamed-Rao modified Mann-Kendall)
# =====================================================================

# ------------------------------ CONFIG -------------------------------
in_csv    <- "MODIS_monthly_BB_RBA_RBF.csv"
out_dir   <- "."
win_lo    <- 0.01     # winsorising percentiles for each component
win_hi    <- 0.99
max_gap   <- 2        # months; gaps up to this length are interpolated
                      # to set the normalisation range only (see header)
roll_win  <- 12       # moving-mean window, months

library(trend)
library(modifiedmk)

# ---------------------------- 1. LOAD --------------------------------
d <- read.csv(in_csv, stringsAsFactors = FALSE)
d$ym <- as.Date(d$ym)
d <- d[order(d$ym), ]
cat(sprintf("Months in file: %d (%s to %s)\n", nrow(d), min(d$ym), max(d$ym)))

# Complete the monthly sequence and flag the months with no data
full <- data.frame(ym = seq(min(d$ym), max(d$ym), by = "month"))
d    <- merge(full, d, by = "ym", all.x = TRUE)
d$missing <- is.na(d$BB) | is.na(d$RBA) | is.na(d$RBF)
cat(sprintf("Months with no data after completing the sequence: %d of %d (%s)\n",
            sum(d$missing), nrow(d),
            paste(format(d$ym[d$missing], "%b %Y"), collapse = ", ")))

# Interpolate short gaps; longer gaps stay NA
fill_short <- function(x, max_gap) {
  na <- is.na(x)
  if (!any(na)) return(x)
  r <- rle(na)
  y <- approx(seq_along(x)[!na], x[!na], xout = seq_along(x), rule = 2)$y
  keep <- rep(r$values & r$lengths > max_gap, r$lengths)
  y[keep] <- NA
  y
}
for (v in c("BB", "RBA", "RBF")) d[[v]] <- fill_short(d[[v]], max_gap)

# ------------------- 2. NORMALISE AND BUILD BS -----------------------
# The normalisation range is set over the completed monthly sequence.
normalise <- function(x, name) {
  qs <- quantile(x, c(win_lo, win_hi), na.rm = TRUE)
  cat(sprintf("%-4s winsorised at [%.4f, %.4f]\n", name, qs[1], qs[2]))
  y <- pmin(pmax(x, qs[1]), qs[2])
  (y - qs[1]) / (qs[2] - qs[1])
}
cat("\n=== Normalisation ===\n")
d$BB_n  <- normalise(d$BB,  "BB")
d$RBA_n <- normalise(d$RBA, "RBA")
d$RBF_n <- normalise(d$RBF, "RBF")
d$BS    <- d$BB_n + d$RBA_n + d$RBF_n

# Interpolated months are excluded from all statistics and outputs.
for (v in c("BB", "RBA", "RBF", "BB_n", "RBA_n", "RBF_n", "BS")) {
  d[[v]][d$missing] <- NA
}

# --------------------- 3. TREND AND CHANGE POINT ---------------------
ok <- !is.na(d$BS)
x  <- d$BS[ok]; t <- d$ym[ok]
cat(sprintf("\nMonths used in the tests: %d\n", length(x)))

cat("\n=== Modified Mann-Kendall (Hamed-Rao) on BS ===\n")
mk <- mmkh(as.numeric(x))
print(round(mk, 6))

cat("\n=== Pettitt change point on BS ===\n")
pt <- pettitt.test(x)
cp_idx  <- as.integer(pt$estimate[1])
cp_date <- t[cp_idx]
cat(sprintf("Change point: %s (p = %.3g)\n", format(cp_date, "%b %Y"), pt$p.value))

before <- ok & d$ym <  cp_date
after  <- ok & d$ym >= cp_date
pct    <- function(a, b) 100 * (b - a) / a

pre  <- mean(d$BS[before]); post <- mean(d$BS[after])
cat(sprintf("Mean BS before: %.3f | after: %.3f | change: %+.1f%%\n",
            pre, post, pct(pre, post)))
cat(sprintf("Median BS before: %.3f | after: %.3f\n",
            median(d$BS[before]), median(d$BS[after])))

cat("\n=== Normalised components: change between periods and trend ===\n")
for (v in c("RBA_n", "RBF_n", "BB_n")) {
  a  <- mean(d[[v]][before]); b <- mean(d[[v]][after])
  mv <- mmkh(as.numeric(d[[v]][ok]))
  cat(sprintf("%-6s mean %.3f -> %.3f (%+.0f%%) | Hamed-Rao tau = %+.3f, p = %.3g\n",
              v, a, b, pct(a, b), mv["Tau"], mv["new P-value"]))
}

# Moving mean, centred
roll_mean <- function(x, k) {
  n <- length(x); y <- rep(NA_real_, n); h <- floor(k / 2)
  for (i in seq_len(n)) {
    j <- max(1, i - h):min(n, i + h)
    if (sum(!is.na(x[j])) >= k / 2) y[i] <- mean(x[j], na.rm = TRUE)
  }
  y
}
d$BS_roll <- roll_mean(d$BS, roll_win)

# --------------------------- 4. FIGURE -------------------------------
col_bb  <- "#d95f02"; col_rba <- "#1b9e77"; col_rbf <- "#7570b3"

png(file.path(out_dir, "fig_temporal_severity.png"),
    width = 9, height = 6.5, units = "in", res = 300, type = "cairo")
op <- par(mfrow = c(2, 1), mar = c(3.2, 4.2, 1.8, 1), mgp = c(2.4, 0.7, 0))

# (a) BS with change point
plot(d$ym, d$BS, type = "l", col = "grey65", lwd = 0.9,
     xlab = "", ylab = "Bloom severity (BS)")
lines(d$ym, d$BS_roll, col = "black", lwd = 2)
abline(v = cp_date, lty = 2, lwd = 1.4)
segments(min(d$ym), pre,  cp_date,   pre,  col = "#c51b7d", lwd = 2)
segments(cp_date,   post, max(d$ym), post, col = "#c51b7d", lwd = 2)
text(cp_date, par("usr")[4], paste0("  ", format(cp_date, "%b %Y")),
     adj = c(0, 1.3), cex = 0.85)
legend("topright", c("monthly", paste0(roll_win, "-month mean"), "period mean"),
       col = c("grey65", "black", "#c51b7d"), lwd = c(0.9, 2, 2),
       bty = "n", cex = 0.8)
mtext("(a)", side = 3, line = 0.4, adj = 0, font = 2)

# (b) normalised components, smoothed with the same 12-month window as (a).
# The raw monthly series are dominated by the seasonal cycle and overlap
# into an unreadable band; smoothing exposes the level shifts, the
# co-movement of RBA and RBF, and the divergence of BB. RBF is dashed so
# that where it sits on top of RBA the overlap remains visible.
d$BB_r  <- roll_mean(d$BB_n,  roll_win)
d$RBA_r <- roll_mean(d$RBA_n, roll_win)
d$RBF_r <- roll_mean(d$RBF_n, roll_win)

plot(d$ym, d$RBA_r, type = "l", col = col_rba, lwd = 2,
     ylim = range(c(d$BB_r, d$RBA_r, d$RBF_r), na.rm = TRUE),
     xlab = "", ylab = paste0("Normalised component (", roll_win, "-month mean)"))
lines(d$ym, d$RBF_r, col = col_rbf, lwd = 2, lty = 2)
lines(d$ym, d$BB_r,  col = col_bb,  lwd = 2)
abline(v = cp_date, lty = 2, lwd = 1.2)
legend("topright", c("BB (intensity)", "RBA (area)", "RBF (frequency)"),
       col = c(col_bb, col_rba, col_rbf), lwd = 2, lty = c(1, 1, 2),
       bty = "n", cex = 0.8, horiz = TRUE)
mtext("(b)", side = 3, line = 0.4, adj = 0, font = 2)

par(op); dev.off()

write.csv(d[, c("ym", "missing", "n_img", "BB", "RBA", "RBF",
                "BB_n", "RBA_n", "RBF_n", "BB_r", "RBA_r", "RBF_r",
                "BS", "BS_roll")],
          file.path(out_dir, "MODIS_monthly_severity.csv"), row.names = FALSE)
cat("\n>>> Written: fig_temporal_severity.png, MODIS_monthly_severity.csv\n")
