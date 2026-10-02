# Open geocomputational workflow for long-term assessment of harmful algal bloom severity

Google Earth Engine (JavaScript) and R code supporting the paper:

> Tiberi, A.E., Drozd, A.A., Zerda, H.R., de Tezanos, P.P. An open geocomputational workflow for enabling long-term and large-scale
> assessment of harmful algal bloom severity. Hydrobiologia (2026). https://doi.org/10.1007/s10750-026-06415-5

The workflow reconstructs long-term spatial and temporal patterns of harmful algal
blooms from open satellite archives, using the Floating Algae Index (FAI) on MODIS
imagery with Sentinel-2 and Landsat-8 as higher-resolution references. It is
demonstrated on a 26-year daily record (2000–2025) of the Río Hondo reservoir,
Argentina, a system with no prior chlorophyll-a monitoring.

## Availability 

⚠️ **Scripts are currently curated and will be available soon** ⚠️
Zenodo repository will be also available soon. 

## Contents

| Script | Purpose |
|---|---|
| `01_MODIS_FAI_daily_metrics.js` | Daily MODIS metrics: relative bloom area, bloom biomass, valid-cell counts |
| `02_MODIS_daily_to_monthly.R` | Coverage filtering, monthly aggregation, relative bloom frequency |
| `03_MODIS_spatial_RBF_BB.js` | Per-cell spatial composites (observations, bloom observations, RBF, BB) |
| `04_inspect_spatial_raster.R` | Raster inspection and consistency checks |
| `05_spatial_invariance_test.R` | Spatial invariance test across months |
| `06_fig_spatial_products.R` | Figure: spatial composites |
| `07_supplementary_component_structure.R` | Correlation structure among severity components |
| `08_fig_temporal_severity.R` | Figure: temporal severity and change point |
| `09_fig_bloom_categories.R` | Figure: monthly bloom categories |
| `10_fig_raster_grids.R` | Figures: annual and monthly raster grids |
| `11_mask_edges_and_islet.R` | Post-hoc exclusion of edge cells and the seasonally exposed islet |
| `12_GEE_daily_FAI_event.js` | Daily FAI for a single bloom event |
| `13_fig_daily_FAI_event.R` | Figure: daily FAI event sequence |

Validation and calibration scripts will be provided separately under `validation/`.

## Requirements

- A Google Earth Engine account for the `.js` scripts, run in the Code Editor.
- R (≥ 4.0) for the `.R` scripts.

Data sources are the MOD09GQ and MOD09GA Collection 6.1 products, Sentinel-2
Level-2A and Landsat-8 Level-2, all accessed through Google Earth Engine. No
proprietary data are required.

## Citation

If you use this code, please cite the paper above.

## License

This code is licensed under the
[GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html).

You are free to use, modify and redistribute it, including for commercial
purposes, provided that any derivative work is distributed under the same
license and its source code is made available. See the `LICENSE` file for
the full terms. 

