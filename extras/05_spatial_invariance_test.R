# =====================================================================
# 05_spatial_invariance_test.R
# ---------------------------------------------------------------------
# Seasonal invariance of the spatial bloom pattern across the twelve
# monthly composites. Script 5 of the workflow (optional: exploratory
# support for the Results section "Seasonality of the spatial pattern";
# no value reported in the paper is computed here).
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
# Question: do the twelve monthly composites show the SAME spatial
# pattern at different intensities, or does the pattern itself
# reorganise through the year?
#
#   - High pairwise correlations across months indicate a seasonally
#     invariant pattern: seasonality changes how much, not where.
#   - Months with low correlation to the common pattern indicate a
#     seasonal spatial reorganisation.
#
# Method: cells valid in ALL twelve months are compared across months
# (Spearman rank correlation), for each variable separately. The monthly
# mean is reported alongside, to separate intensity ("how much") from
# spatial structure ("where").
#
# Input : the twelve monthly composites from script 03 (full record,
#         CALENDAR_MONTH = 1 ... 12). Preferably use the versions with
#         edge cells and the islet excluded (script 12), since the
#         seasonally exposed islet would otherwise introduce a spurious
#         seasonal signal.
#
# Requires: terra
# =====================================================================

# ------------------------------ CONFIG -------------------------------
in_dir   <- "."
file_tpl <- "MODIS_spatial_RBF_BB_20000101_20260101_m%d.tif"  # %d = month
out_dir  <- "."
bands    <- c("RBF", "BB")                        # variables to test
band_order <- c("n_obs", "n_bloom", "RBF", "BB")  # band order in the GeoTIFF

# Correlation above which a month is considered to track the common
# pattern. Descriptive cut-off, not a statistical test.
inv_thr  <- 0.90

month_lab <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

library(terra)

# ---------------------------- 1. LOAD --------------------------------
paths <- file.path(in_dir, sprintf(file_tpl, 1:12))
missing <- !file.exists(paths)
if (any(missing)) stop("Missing files: ", paste(basename(paths[missing]), collapse = ", "))

rasters <- lapply(paths, function(p) {
  r <- rast(p)
  if (!all(band_order %in% names(r))) {
    if (nlyr(r) == length(band_order)) names(r) <- band_order
    else stop("Unexpected band count in ", basename(p))
  }
  r
})

# All months must share the same grid for a cell-wise comparison
geom_ref <- c(nrow(rasters[[1]]), ncol(rasters[[1]]))
for (i in 2:12) {
  if (!all(c(nrow(rasters[[i]]), ncol(rasters[[i]])) == geom_ref))
    stop("Raster geometry differs in month ", i)
}
cat(sprintf("Loaded 12 monthly composites (%d x %d cells)\n", geom_ref[1], geom_ref[2]))

# ------------------- 2. TEST, ONE VARIABLE AT A TIME -----------------
results <- list()

for (b in bands) {
  cat(sprintf("\n=============== %s ===============\n", b))

  # cells x months matrix
  M <- sapply(rasters, function(r) values(r[[b]], mat = FALSE))
  colnames(M) <- month_lab

  # keep only cells valid in every month, so all pairs use the same cells
  complete <- apply(M, 1, function(x) all(is.finite(x)))
  M <- M[complete, , drop = FALSE]
  cat(sprintf("Cells valid in all 12 months: %d\n", nrow(M)))
  if (nrow(M) < 50) { cat("Too few common cells; skipping.\n"); next }

  # monthly intensity: how much, independent of where
  intensity <- colMeans(M)
  cat("\nMonthly mean (intensity):\n")
  print(round(intensity, 4))

  # pairwise spatial correlation: where, independent of how much
  C <- cor(M, method = "spearman")
  off <- C[upper.tri(C)]
  cat(sprintf("\nPairwise spatial correlation (Spearman, %d month pairs):\n", length(off)))
  cat(sprintf("  min = %.3f | median = %.3f | max = %.3f\n",
              min(off), median(off), max(off)))
  cat(sprintf("  pairs above %.2f: %d of %d (%.0f%%)\n",
              inv_thr, sum(off > inv_thr), length(off),
              100 * mean(off > inv_thr)))

  # correlation of each month with the mean pattern across months
  mean_map <- rowMeans(M)
  vs_mean  <- apply(M, 2, function(x) cor(x, mean_map, method = "spearman"))
  cat("\nCorrelation of each month with the average spatial pattern:\n")
  print(round(vs_mean, 3))

  departing <- names(vs_mean)[vs_mean < inv_thr]
  if (length(departing) == 0) {
    cat(sprintf("\n=> All months track the common pattern (all >= %.2f):\n", inv_thr))
    cat("   the spatial structure is seasonally invariant for this variable.\n")
  } else {
    cat(sprintf("\n=> Months departing from the common pattern: %s\n",
                paste(departing, collapse = ", ")))
    cat("   The spatial structure reorganises seasonally for this variable.\n")
  }

  results[[b]] <- list(C = C, intensity = intensity, vs_mean = vs_mean)
  write.csv(round(C, 4), file.path(out_dir, sprintf("invariance_cormat_%s.csv", b)))
}

# -------------------------- 3. FIGURE --------------------------------
png(file.path(out_dir, "spatial_invariance.png"),
    width = 5 * length(results), height = 7.5, units = "in", res = 300,
    type = "cairo")
op <- par(mfrow = c(2, length(results)), mar = c(4, 4, 3, 4))

heat <- colorRampPalette(c("#f7fbff", "#9ecae1", "#2171b5", "#08306b"))(100)

# top row: correlation matrices
for (b in names(results)) {
  C <- results[[b]]$C
  image(1:12, 1:12, C[, 12:1], col = heat, zlim = c(min(C), 1),
        axes = FALSE, xlab = "", ylab = "",
        main = sprintf("%s - pairwise spatial correlation", b))
  axis(1, at = 1:12, labels = month_lab, las = 2, cex.axis = 0.7)
  axis(2, at = 1:12, labels = rev(month_lab), las = 2, cex.axis = 0.7)
  box()
}

# bottom row: monthly intensity
for (b in names(results)) {
  int <- results[[b]]$intensity
  plot(1:12, int, type = "b", pch = 19, col = "#d95f02", lwd = 1.6,
       xaxt = "n", xlab = "", ylab = sprintf("Mean %s", b),
       main = sprintf("%s - seasonal intensity", b))
  axis(1, at = 1:12, labels = month_lab, las = 2, cex.axis = 0.8)
}

par(op); dev.off()

cat("\n>>> Written: spatial_invariance.png and invariance_cormat_*.csv\n")
