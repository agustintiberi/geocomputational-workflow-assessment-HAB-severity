# =====================================================================
# 06_fig_spatial_products.R
# ---------------------------------------------------------------------
# Long-term spatial products (BB, RBF, BS) over the full record:
# normalisation, spatial bloom severity and Figure 9 of the paper.
# Script 6 of the workflow.
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
# Input : full-record multi-band GeoTIFF from 03_MODIS_spatial_RBF_BB.js
#         (START = 2000-01-01, END = 2026-01-01, CALENDAR_MONTH = null),
#         AFTER exclusion of edge cells and of the seasonally exposed
#         islet by script 12 (written to the 'masked/' folder).
#         >>> Run script 12 before this one. <<<
#         bands: n_obs, n_bloom, RBF, BB
#
# Steps:
#   BB_n   spatial bloom biomass, winsorised (1st/99th pct) and rescaled
#          to 0-1 (Eqs. 10a and 10b)
#   RBF_n  spatial bloom frequency, winsorised and rescaled to 0-1
#   BS     spatial bloom severity = BB_n + RBF_n, range 0-2 (Eq. 11b)
#
# Spatial BS uses two components only: RBA is a reservoir-wide quantity
# (bloom fraction of the water body on a given date) with no per-cell
# counterpart, so it enters the temporal severity index (Eq. 11a) but
# not the spatial one.
#
# Outputs:
#   fig_spatial_products.png   Figure 9 (three panels)
#   MODIS_spatial_BS.tif       BB_n, RBF_n and BS rasters
#
# Requires: terra
# =====================================================================

# ------------------------------ CONFIG -------------------------------
tif_path  <- "masked/MODIS_spatial_RBF_BB_20000101_20260101.tif"  # from script 12
aoi_path  <- "data/aoi.geojson"   # reservoir outline drawn on the maps; NA to skip
out_dir   <- "."
out_png   <- "fig_spatial_products.png"
out_tif   <- "MODIS_spatial_BS.tif"

win_lo    <- 0.01        # winsorising percentiles
win_hi    <- 0.99
panel_lab <- c("(a)", "(b)", "(c)")   # set to c("", "", "") to omit
sbar_km   <- 6           # scale bar length, km

library(terra)

# ---------------------------- 1. LOAD --------------------------------
r <- rast(tif_path)
expected <- c("n_obs", "n_bloom", "RBF", "BB")
if (!all(expected %in% names(r))) {
  if (nlyr(r) == length(expected)) names(r) <- expected
  else stop(sprintf("Expected %d bands, found %d.", length(expected), nlyr(r)))
}
cat(sprintf("Raster loaded: %d x %d cells, %d bands\n", nrow(r), ncol(r), nlyr(r)))

outline <- if (!is.na(aoi_path)) vect(aoi_path) else NULL

# ------------------ 2. WINSORISE AND RESCALE TO 0-1 ------------------
# Clipping at the 1st/99th percentiles before rescaling prevents a few
# extreme cells (mostly along the shoreline) from compressing the range
# of the whole map (Eqs. 10a and 10b).
normalise <- function(x, lo = win_lo, hi = win_hi, name = "") {
  v  <- values(x, mat = FALSE)
  v  <- v[is.finite(v)]
  qs <- quantile(v, c(lo, hi), na.rm = TRUE)
  cat(sprintf("%-4s winsorised at [%.4f, %.4f]; raw range [%.4f, %.4f]\n",
              name, qs[1], qs[2], min(v), max(v)))
  y <- clamp(x, qs[1], qs[2], values = TRUE)
  (y - qs[1]) / (qs[2] - qs[1])
}

cat("\n=== Normalisation ===\n")
BB_n  <- normalise(r[["BB"]],  name = "BB")
RBF_n <- normalise(r[["RBF"]], name = "RBF")

# ------------------------- 3. SEVERITY -------------------------------
# Equally weighted sum of the two normalised components (range 0-2; Eq. 11b).
BS <- BB_n + RBF_n
names(BB_n) <- "BB_n"; names(RBF_n) <- "RBF_n"; names(BS) <- "BS"

cat("\n=== Summary of plotted layers ===\n")
for (lyr in list(BB_n, RBF_n, BS)) {
  v <- values(lyr, mat = FALSE); v <- v[is.finite(v)]
  cat(sprintf("%-6s n = %d | min %.3f | median %.3f | mean %.3f | max %.3f\n",
              names(lyr), length(v), min(v), median(v), mean(v), max(v)))
}

writeRaster(c(BB_n, RBF_n, BS), file.path(out_dir, out_tif), overwrite = TRUE)

# --------------------------- 4. FIGURE -------------------------------
# Distinct palettes so the composite is not mistaken for a component.
pal_bb  <- colorRampPalette(c("#1a9850", "#a6d96a", "#ffffbf", "#fdae61", "#d73027"))(100)
pal_rbf <- colorRampPalette(c("#fff5f0", "#fcbba1", "#fb6a4a", "#cb181d", "#67000d"))(100)
pal_bs  <- colorRampPalette(c("#3288bd", "#99d594", "#e6f598", "#fee08b", "#fc8d59", "#d53e4f"))(100)

png(file.path(out_dir, out_png),
    width = 12, height = 4.4, units = "in", res = 300, type = "cairo")
op <- par(mfrow = c(1, 3), mar = c(2, 2, 2.5, 4.5))

# Scale bar drawn explicitly. terra::sbar() takes the length in MAP UNITS,
# which for a lon/lat raster means degrees, not metres; drawing it here
# keeps the length correct and the placement predictable.
scale_bar <- function(x, km, frac_x = 0.50, frac_y = 0.07) {
  e    <- as.vector(ext(x))                       # xmin, xmax, ymin, ymax
  lat  <- mean(c(e[3], e[4]))
  deg  <- km / (111.32 * cos(lat * pi / 180))     # degrees of longitude per km
  x0   <- e[1] + (e[2] - e[1]) * frac_x
  y0   <- e[3] + (e[4] - e[3]) * frac_y
  tick <- (e[4] - e[3]) * 0.018
  segments(x0, y0, x0 + deg, y0, lwd = 5, lend = 1)
  segments(c(x0, x0 + deg), y0 - tick, c(x0, x0 + deg), y0 + tick,
           lwd = 5, lend = 1)
  text(x0 + deg / 2, y0, sprintf("%g km", km), pos = 3, cex = 1.1,
       font = 2, offset = 0.45)
}

# North arrow. Geographic north is up in EPSG:4326, so a vertical arrow is
# correct without further rotation.
north_arrow <- function(x, frac_x = 0.90, frac_y = 0.80, size = 0.11) {
  e  <- as.vector(ext(x))
  w  <- e[2] - e[1]; h <- e[4] - e[3]
  cx <- e[1] + w * frac_x
  y0 <- e[3] + h * frac_y
  y1 <- y0 + h * size
  arrows(cx, y0, cx, y1, length = 0.10, lwd = 3, lend = 1)
  text(cx, y1, "N", pos = 3, font = 2, cex = 1.2, offset = 0.25)
}

draw <- function(x, pal, main, lab, rng, add_sbar = FALSE, add_north = FALSE) {
  plot(x, col = pal, range = rng, main = main, axes = FALSE, mar = NA,
       plg = list(title = "", cex = 0.8))
  if (!is.null(outline)) plot(outline, add = TRUE, border = "grey40", lwd = 0.7)
  if (nchar(lab)) mtext(lab, side = 3, line = 0.4, adj = 0, font = 2, cex = 1)
  if (add_north) north_arrow(x)
  if (add_sbar)  scale_bar(x, sbar_km)
}

draw(BB_n,  pal_bb,  "Bloom biomass (BB)",    panel_lab[1], c(0, 1),
     add_north = TRUE)
draw(RBF_n, pal_rbf, "Bloom frequency (RBF)", panel_lab[2], c(0, 1))
draw(BS,    pal_bs,  "Bloom severity (BS)",   panel_lab[3], c(0, 2),
     add_sbar = TRUE)

par(op); dev.off()

cat(sprintf("\n>>> Written: %s and %s\n", out_png, out_tif))
