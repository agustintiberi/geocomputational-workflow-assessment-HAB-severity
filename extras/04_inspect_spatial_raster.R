# =====================================================================
# 04_inspect_spatial_raster.R
# ---------------------------------------------------------------------
# Quality-control inspection of a spatial product exported by script 03.
# Script 4 of the workflow (optional: it checks the products but does
# not generate any result reported in the paper).
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
# Input : multi-band GeoTIFF exported by 03_MODIS_spatial_RBF_BB.js
#         bands (in order): n_obs, n_bloom, RBF, BB
#
# Purpose:
#   1. Report raster geometry and per-band statistics
#   2. Run consistency checks (RBF in [0, 1]; n_bloom <= n_obs;
#      RBF equals n_bloom / n_obs)
#   3. Test whether observation frequency still drives the RBF pattern,
#      which is the bias the per-cell normalisation is meant to remove
#   4. Write a four-panel map and an RBF histogram
#
# Run it on any product (full record, annual or monthly) before using it
# in scripts 06, 10 or 12.
#
# Requires: terra
# =====================================================================

# ------------------------------ CONFIG -------------------------------
# Example: September composite pooled over the full record.
tif_path <- "MODIS_spatial_RBF_BB_20000101_20260101_m9.tif"  # adjust
out_dir  <- "."
label    <- "September composite (2000-2025)"                # adjust

library(terra)

# ---------------------------- 1. LOAD --------------------------------
r <- rast(tif_path)

# GEE sometimes drops band names on export; assign them by position.
expected <- c("n_obs", "n_bloom", "RBF", "BB")
if (!all(expected %in% names(r))) {
  if (nlyr(r) == length(expected)) {
    names(r) <- expected
    cat("Band names were missing and have been assigned by position.\n")
  } else {
    stop(sprintf("Expected %d bands, found %d.", length(expected), nlyr(r)))
  }
}

cat(sprintf("\nRaster: %d x %d cells, %d bands\n", nrow(r), ncol(r), nlyr(r)))
cat(sprintf("Resolution: %.6f x %.6f (CRS units)\n", res(r)[1], res(r)[2]))
cat(sprintf("CRS: %s\n", crs(r, describe = TRUE)$name))
print(ext(r))

# ------------------------ 2. BAND STATISTICS -------------------------
cat("\n=== Per-band statistics (unmasked cells only) ===\n")
stats <- do.call(rbind, lapply(expected, function(b) {
  v <- values(r[[b]], mat = FALSE)
  v <- v[is.finite(v)]
  data.frame(band = b, n_cells = length(v),
             min = round(min(v), 4), median = round(median(v), 4),
             mean = round(mean(v), 4), max = round(max(v), 4))
}))
print(stats, row.names = FALSE)

n_masked <- sum(!is.finite(values(r[["n_obs"]], mat = FALSE)))
cat(sprintf("\nMasked cells (outside the AOI or below MIN_OBS): %d\n", n_masked))

# ------------------------ 3. CONSISTENCY CHECKS ----------------------
cat("\n=== Consistency checks ===\n")
v_obs <- values(r[["n_obs"]],   mat = FALSE)
v_blo <- values(r[["n_bloom"]], mat = FALSE)
v_rbf <- values(r[["RBF"]],     mat = FALSE)
ok    <- is.finite(v_obs) & is.finite(v_blo) & is.finite(v_rbf)

bad_range  <- sum(v_rbf[ok] < 0 | v_rbf[ok] > 1)
bad_count  <- sum(v_blo[ok] > v_obs[ok])
recomputed <- v_blo[ok] / v_obs[ok]
max_dev    <- max(abs(recomputed - v_rbf[ok]))

cat(sprintf("RBF outside [0, 1]         : %d cells %s\n", bad_range,
            ifelse(bad_range == 0, "(OK)", "<-- CHECK")))
cat(sprintf("n_bloom greater than n_obs : %d cells %s\n", bad_count,
            ifelse(bad_count == 0, "(OK)", "<-- CHECK")))
cat(sprintf("Max |RBF - n_bloom/n_obs|  : %.2e %s\n", max_dev,
            ifelse(max_dev < 1e-5, "(OK)", "<-- CHECK")))

# --------------- 4. DOES OBSERVATION COUNT DRIVE RBF? ----------------
# A strong association would indicate that the spatial pattern still
# reflects how often each cell is observed rather than how often it blooms.
rho <- cor(v_obs[ok], v_rbf[ok], method = "spearman")
cat(sprintf("\nSpearman rho, n_obs vs RBF: %.3f\n", rho))
cat("Interpretation: values near zero indicate that the normalisation\n")
cat("removed the observation-frequency bias; a strong positive or negative\n")
cat("value would indicate residual sampling structure worth examining.\n")

cat(sprintf("\nObservation range across cells: %.0f to %.0f (ratio %.1fx)\n",
            min(v_obs[ok]), max(v_obs[ok]), max(v_obs[ok]) / min(v_obs[ok])))

# -------------------------- 5. FIGURES -------------------------------
pal <- colorRampPalette(c("#171cb1", "#9ef9fd", "#41ce13",
                          "#eefd0e", "#f33214", "#b41f07"))(100)

png(file.path(out_dir, "spatial_product_panels.png"),
    width = 9, height = 8, units = "in", res = 300, type = "cairo")
op <- par(mfrow = c(2, 2), mar = c(2.5, 2.5, 2.5, 4))
plot(r[["n_obs"]],   col = pal, main = paste("n_obs -", label))
plot(r[["n_bloom"]], col = pal, main = "n_bloom")
plot(r[["RBF"]],     col = pal, main = "RBF (bloom frequency)", range = c(0, 1))
plot(r[["BB"]],      col = pal, main = "BB (median FAI over bloom observations)")
par(op); dev.off()

png(file.path(out_dir, "spatial_RBF_histogram.png"),
    width = 6, height = 4, units = "in", res = 300, type = "cairo")
hist(v_rbf[ok], breaks = 40, col = "#2c7fb8", border = "white",
     main = paste("RBF distribution -", label),
     xlab = "Bloom frequency (0-1)")
dev.off()

cat("\n>>> Written: spatial_product_panels.png, spatial_RBF_histogram.png\n")
