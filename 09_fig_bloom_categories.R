# =====================================================================
# 09_fig_bloom_categories.R
# ---------------------------------------------------------------------
# Composition of bloom categories by calendar month and by year:
# Figure 6 and Table S4 of the paper. Script 9 of the workflow.
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
# Input : MODIS_daily_filtered.csv (from 02_MODIS_daily_to_monthly.R)
#         columns used: date, RBA
#
# Each usable date is classified by the fraction of the reservoir in bloom:
#   Significant : RBA >  0.25
#   Slight      : 0.05 < RBA <= 0.25
#   Null        : RBA <= 0.05
#
# The "Significant" class is the criterion behind the temporal Relative
# Bloom Frequency: RBF is the fraction of usable dates in a month
# classified as Significant (Eq. 7). Only the 0.25 boundary therefore
# propagates into the severity analysis; the 0.05 boundary affects the
# descriptive composition of categories only.
#
# Outputs:
#   fig_categories_monthly.png          Figure 6 (all years pooled)
#   TableS4a_categories_annual.csv      Table S4a (proportions by year)
#   TableS4b_categories_monthly.csv     Table S4b (proportions by month)
#   fig_categories_annual.png           year-by-year version of Figure 6
#                                       (not shown in the paper; the same
#                                       information is given in Table S4a)
#
# Base R only; PNG written through Cairo.
# =====================================================================

# ------------------------------ CONFIG -------------------------------
in_csv   <- "MODIS_daily_filtered.csv"
out_dir  <- "."
thr_sig  <- 0.25    # RBA above this -> Significant
thr_null <- 0.05    # RBA at or below this -> Null

col_sig  <- "#E03A32"   # Significant
col_sli  <- "#0E8B3D"   # Slight
col_nul  <- "#4A9BD4"   # Null

month_lab <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# ---------------------------- 1. LOAD --------------------------------
d <- read.csv(in_csv, stringsAsFactors = FALSE)
d$date <- as.Date(d$date)
d$RBA  <- suppressWarnings(as.numeric(d$RBA))
d <- d[is.finite(d$RBA) & !is.na(d$date), ]

cat(sprintf("Usable dates: %d (%s to %s)\n", nrow(d), min(d$date), max(d$date)))

# ------------------------ 2. CLASSIFY DATES --------------------------
d$category <- cut(d$RBA,
                  breaks = c(-Inf, thr_null, thr_sig, Inf),
                  labels = c("Null", "Slight", "Significant"),
                  right  = TRUE)

cat("\n=== Overall composition (%) ===\n")
print(round(100 * prop.table(table(d$category)), 1))

d$year  <- as.integer(format(d$date, "%Y"))
d$month <- as.integer(format(d$date, "%m"))

# Stacking order: Significant at the bottom, Null on top
ord  <- c("Significant", "Slight", "Null")
cols <- c(col_sig, col_sli, col_nul)

# Helper: proportion matrix (categories x groups) plus the group sample size
compose <- function(group) {
  tb <- table(factor(d$category, levels = ord), group)
  n  <- colSums(tb)
  list(p = sweep(tb, 2, n, "/"), n = n)
}

# ------------------- 3. BY CALENDAR MONTH (FIGURE 6) -----------------
mo <- compose(d$month)
colnames(mo$p) <- month_lab[as.integer(colnames(mo$p))]

cat("\n=== Usable dates per calendar month (all years pooled) ===\n")
print(mo$n)

png(file.path(out_dir, "fig_categories_monthly.png"),
    width = 8, height = 4.2, units = "in", res = 300, type = "cairo")
par(mar = c(3.4, 4.2, 1.2, 8.5), mgp = c(2.4, 0.7, 0), xpd = NA)
bp <- barplot(mo$p, col = cols, border = NA, space = 0.35,
              ylab = "Proportion of usable dates", xlab = "Month",
              ylim = c(0, 1), las = 1, cex.names = 0.9)
legend(x = max(bp) + 1.2, y = 1,
       legend = ord, fill = cols, border = NA, bty = "n",
       title = "Bloom category", cex = 0.9)
dev.off()

# --------------------------- 4. BY YEAR ------------------------------
yr <- compose(d$year)

cat("\n=== Usable dates per year ===\n")
print(yr$n)

png(file.path(out_dir, "fig_categories_annual.png"),
    width = 9, height = 4.2, units = "in", res = 300, type = "cairo")
par(mar = c(3.8, 4.2, 1.2, 8.5), mgp = c(2.4, 0.7, 0), xpd = NA)
bp <- barplot(yr$p, col = cols, border = NA, space = 0.3,
              ylab = "Proportion of usable dates", xlab = "",
              ylim = c(0, 1), las = 2, cex.names = 0.75)
mtext("Year", side = 1, line = 2.6)
legend(x = max(bp) + 1.5, y = 1,
       legend = ord, fill = cols, border = NA, bty = "n",
       title = "Bloom category", cex = 0.9)
dev.off()

# ---------------------------- 5. TABLES ------------------------------
# Proportions (0-1); Table S4 reports them as percentages.
write.csv(round(t(yr$p), 4), file.path(out_dir, "TableS4a_categories_annual.csv"))
write.csv(round(t(mo$p), 4), file.path(out_dir, "TableS4b_categories_monthly.csv"))

cat("\n>>> Written: fig_categories_monthly.png, fig_categories_annual.png,",
    "TableS4a_categories_annual.csv, TableS4b_categories_monthly.csv\n")
