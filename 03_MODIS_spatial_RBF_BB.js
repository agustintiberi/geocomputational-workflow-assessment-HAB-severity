/* =====================================================================
   03_MODIS_spatial_RBF_BB.js
   ---------------------------------------------------------------------
   Per-cell spatial bloom products for a selected period.
   Script 3 of the workflow.

   Associated publication:
     Tiberi, A. E., A. A. Drozd, H. R. Zerda & P. de Tezanos Pinto, 2026.
     An open geocomputational workflow for enabling long-term and
     large-scale assessment of harmful algal bloom severity.
     Hydrobiologia. https://doi.org/10.1007/s10750-026-06415-5

   Copyright (C) 2026 Agustín Eduardo Tiberi
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or (at
   your option) any later version. See the LICENSE file for details.
   ---------------------------------------------------------------------
   Produces a single multi-band raster for the selected period:

     n_obs    number of valid (unmasked) observations per cell
     n_bloom  number of observations with FAI >= FAI_THR per cell
     RBF      n_bloom / n_obs -> spatial relative bloom frequency (0-1; Eq. 9)
     BB       median FAI over bloom observations only -> spatial bloom
              biomass (Eq. 8)

   Why RBF is normalised per cell:
   cells differ in how often they are observed (clouds, masking, reservoir
   level). A raw SUM of bloom occurrences therefore mixes "blooms often"
   with "observed often". Dividing by each cell's own n_obs removes that
   bias and yields a comparable frequency. n_obs and n_bloom are exported
   as well so that the ratio can be audited.

   Products used in the paper:
     - full-record composite: START = '2000-01-01', END = '2026-01-01'
     - annual composites:     one run per year (START = 'YYYY-01-01',
                              END = 'YYYY+1-01-01')
     - monthly composites:    full record with CALENDAR_MONTH = 1 ... 12

   Exclusion of edge cells and of the seasonally exposed islet is applied
   afterwards, in R (script 12). Normalisation to 0-1 and the spatial
   severity index are computed in R (script 06).

   Requires an imported geometry asset named 'AOI': the minimum
   (permanent) reservoir extent polygon.
   ===================================================================== */

//============================ PARAMETERS ==============================
var START = '2000-01-01';   // period start (inclusive)
var END   = '2026-01-01';   // period end   (exclusive)

// Calendar-month filter: null = all months; 9 = September only, etc.
// Used to build monthly composites pooled over the whole record.
var CALENDAR_MONTH = null;

var FAI_THR       = -0.001;  // bloom / non-bloom threshold
var CLOUD_B01_MAX = 1500;    // cloud mask (b01, red 620-670 nm; raw value)
var LAND_B07_MAX  = 1800;    // land mask (b07, SWIR3 2105-2155 nm; raw value)

// Per-image coverage filter (same criterion as the temporal series,
// script 02). Set to false if the export runs out of memory.
// >>> CHECK: must match the setting used for the published products. <<<
var USE_COVERAGE_FILTER = true;
var MIN_COV_FRAC        = 0.65;  // fraction of the reference extent per image

// Cells observed fewer than MIN_OBS times are masked out of RBF and BB,
// because a frequency computed from very few observations is unstable.
var MIN_OBS = 10;

// Reducer for the BB raster: 'median' (the definition used in the paper,
// consistent with the temporal BB) or 'mean' (lighter on memory over
// long periods, but NOT the published definition).
var BB_REDUCER = 'median';

var SCALE = 250;
var CRS   = 'EPSG:4326';
var aoi   = ee.FeatureCollection(AOI).geometry();

Map.centerObject(aoi, 12);

//========================= INPUT COLLECTIONS ==========================
// MOD09GQ 250 m: b01 = red (645 nm), b02 = NIR (859 nm)
var MODGQ = ee.ImageCollection('MODIS/061/MOD09GQ')
  .filterBounds(aoi).filterDate(START, END)
  .select(['sur_refl_b01', 'sur_refl_b02'], ['b01', 'b02']);

// MOD09GA 500 m: b05 = SWIR1 (1240 nm), b06 = SWIR2 (1640 nm),
// b07 = SWIR3 (2130 nm). Resampled to 250 m to align with MOD09GQ.
var MODGA = ee.ImageCollection('MODIS/061/MOD09GA')
  .filterBounds(aoi).filterDate(START, END)
  .select(['sur_refl_b05', 'sur_refl_b06', 'sur_refl_b07'],
          ['b05', 'b06', 'b07'])
  .map(function (img) {
    return img.resample('bilinear').reproject({crs: CRS, scale: SCALE});
  });

if (CALENDAR_MONTH !== null) {
  MODGQ = MODGQ.filter(ee.Filter.calendarRange(CALENDAR_MONTH, CALENDAR_MONTH, 'month'));
  MODGA = MODGA.filter(ee.Filter.calendarRange(CALENDAR_MONTH, CALENDAR_MONTH, 'month'));
}

// combine() pairs images by system:index (the acquisition date in both
// products). Check that both collections have the same size.
print('Number of MOD09GQ images:', MODGQ.size());
print('Number of MOD09GA images:', MODGA.size());

var IC = MODGQ.combine(MODGA);

//==================== ATMOSPHERIC CORRECTION (SWIR) ===================
// Rrsc(i) = Rrc(i) - min(Rrc_SWIR1, Rrc_SWIR2, Rrc_SWIR3)   (Eq. 1)
function atmosphericCorrection(image) {
  var minSWIR = image.select('b05')
    .min(image.select('b06'))
    .min(image.select('b07'))
    .rename('minSWIR');
  var x = minSWIR.multiply(0.0001);
  var B01 = image.select('b01').multiply(0.0001).subtract(x).rename('B01'); // 645 nm
  var B02 = image.select('b02').multiply(0.0001).subtract(x).rename('B02'); // 859 nm
  var B05 = image.select('b05').multiply(0.0001).subtract(x).rename('B05'); // 1240 nm
  return image.addBands([minSWIR, B01, B02, B05]);
}

//=========================== MASKS + FAI ==============================
// Cloud and land masks are applied to the uncorrected reflectance (Rrc).
function computeFAI(image) {
  var img = atmosphericCorrection(image).clip(aoi);
  img = img.updateMask(img.select('b01').lte(CLOUD_B01_MAX))
           .updateMask(img.select('b07').lte(LAND_B07_MAX));
  // FAI (Hu, 2009; Eqs. 2A and 2B)
  var fai = img.expression(
    'nir - (red + (swir - red) * ((859 - 645) / (1240 - 645)))',
    { nir:  img.select('B02'),
      red:  img.select('B01'),
      swir: img.select('B05') }
  ).rename('FAI');
  return ee.Image(fai.copyProperties(image, ['system:time_start']));
}

var IC_FAI = ee.ImageCollection(IC.map(computeFAI));

//===================== PER-IMAGE COVERAGE FILTER ======================
// Tags each image with its number of valid water cells and keeps images
// reaching MIN_COV_FRAC of the reference extent. The reference extent is
// the 99th percentile of the per-image counts within the period, which
// equals the number of cells in the AOI (2154) whenever the period
// contains fully observed dates.
if (USE_COVERAGE_FILTER) {
  IC_FAI = IC_FAI.map(function (img) {
    var n = img.reduceRegion({
      reducer: ee.Reducer.count(), geometry: aoi,
      scale: SCALE, crs: CRS, maxPixels: 1e10
    }).get('FAI');
    return img.set('n_valid', n);
  });

  var fullExtent = ee.Number(IC_FAI.reduceColumns({
    reducer: ee.Reducer.percentile([99]), selectors: ['n_valid']
  }).get('p99'));

  print('Reference extent (99th pct of valid cells):', fullExtent);
  print('Images before coverage filter:', IC_FAI.size());

  IC_FAI = IC_FAI.filter(ee.Filter.gte('n_valid', fullExtent.multiply(MIN_COV_FRAC)));
  print('Images after coverage filter:', IC_FAI.size());
} else {
  print('Images used (no coverage filter):', IC_FAI.size());
}

//========================== SPATIAL PRODUCTS ==========================
// Valid observations per cell
var nObs = IC_FAI.select('FAI').count().rename('n_obs');

// Bloom occurrences per cell (gte preserves the mask, so masked
// observations are not counted as zeros)
var nBloom = IC_FAI.map(function (img) {
  return ee.Image(img).gte(FAI_THR).rename('bloom');
}).sum().rename('n_bloom');

// Spatial relative bloom frequency (Eq. 9)
var RBF = nBloom.divide(nObs).rename('RBF');

// Spatial bloom biomass: FAI over bloom observations only (Eq. 8)
var IC_bloom = IC_FAI.map(function (img) {
  img = ee.Image(img);
  return img.updateMask(img.gte(FAI_THR));
});
var BB = (BB_REDUCER === 'median' ? IC_bloom.median() : IC_bloom.mean()).rename('BB');

// Mask cells with too few observations to support a stable statistic
var stable = nObs.gte(MIN_OBS);
var product = nObs.addBands(nBloom).addBands(RBF).addBands(BB)
                  .updateMask(stable)
                  .toFloat();

//============================ VISUAL CHECK ============================
/*
Map.addLayer(nObs.updateMask(stable), {min: 0, max: 500}, 'n_obs', false);
Map.addLayer(RBF.updateMask(stable),
  {min: 0, max: 1,
   palette: ['171cb1', '9ef9fd', '41ce13', 'eefd0e', 'f33214', 'b41f07']}, 'RBF');
Map.addLayer(BB.updateMask(stable),
  {min: FAI_THR, max: 0.2,
   palette: ['171cb1', '9ef9fd', '41ce13', 'eefd0e', 'f33214', 'b41f07']}, 'BB', false);
*/

//=============================== EXPORT ===============================
var tag = START.replace(/-/g, '') + '_' + END.replace(/-/g, '') +
          (CALENDAR_MONTH !== null ? '_m' + CALENDAR_MONTH : '');

Export.image.toDrive({
  image:          product,
  description:    'MODIS_spatial_RBF_BB_' + tag,
  fileNamePrefix: 'MODIS_spatial_RBF_BB_' + tag,
  region:         aoi.bounds(),
  scale:          SCALE,
  crs:            CRS,
  maxPixels:      1e10
});
