# Open geocomputational workflow for long-term assessment of harmful algal bloom severity

Google Earth Engine (JavaScript) and R code supporting the paper:

> [2026] Tiberi, A. E; Drozd, A. A; Zerda, H. R; de Tezanos Pinto, P. An open geocomputational workflow 
> for enabling long-term and large-scale assessment of harmful algal bloom severity. *Hydrobiologia*.

The workflow reconstructs long-term spatial and temporal patterns of harmful algal
blooms from open satellite archives, using the Floating Algae Index (FAI) on MODIS
imagery with Sentinel-2 and Landsat-8 as higher-resolution references. It is
demonstrated on a 26-year daily record (2000–2025) of the Río Hondo reservoir,
Argentina, a system with no prior chlorophyll-a monitoring.

## Availability 

⚠️ **These scripts will be made publicly available upon publication of the paper.** ⚠️
Until then the repository is provided for editorial and review purposes only.

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
| `12_mask_edges_and_islet.R` | Post-hoc exclusion of edge cells and the seasonally exposed islet |
| `13_GEE_daily_FAI_event.js` | Daily FAI for a single bloom event |
| `14_fig_daily_FAI_event.R` | Figure: daily FAI event sequence |

Validation and calibration scripts are provided separately under `validation/`.

## Requirements

- A Google Earth Engine account for the `.js` scripts, run in the Code Editor.
- R (≥ 4.0) for the `.R` scripts.

Data sources are the MOD09GQ and MOD09GA Collection 6.1 products, Sentinel-2
Level-2A and Landsat-8 Level-2, all accessed through Google Earth Engine. No
proprietary data are required.

## Citation

If you use this code, please cite the paper above.

## License

This work is licensed under CC BY-NC-SA 4.0
https://creativecommons.org/licenses/by-nc-sa/4.0/

