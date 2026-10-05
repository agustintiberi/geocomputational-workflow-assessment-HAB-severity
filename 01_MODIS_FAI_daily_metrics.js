/* =====================================================================
   01_MODIS_FAI_daily_metrics.js
   ---------------------------------------------------------------------
   Daily MODIS Floating Algae Index (FAI) metrics over the area of
   interest (AOI). Script 1 of the workflow.

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
   Output: one row per date, with
     date         image date (YYYY-MM-dd)
     n_valid      number of valid water cells (after masking)
     n_bloom      number of cells with FAI >= FAI_THR
     RBA          n_bloom / n_valid (Relative Bloom Area)
     BB_bloom     median FAI over bloom cells only (Bloom Biomass, BB)
     BB_allwater  median FAI over all valid water cells (diagnostic only)
   ---------------------------------------------------------------------
   NOTE 1: BB_bloom is the Bloom Biomass (BB) variable used in the paper.
   BB_allwater is exported for diagnostic purposes only: it is strongly
   associated with RBA and is not used as a bloom variable (see Methods).

   NOTE 2: no minimum-coverage filter is applied here. n_valid is
   exported and the filter is applied in R (script 02), where it is
   explicit and easy to audit.

   Requires an imported geometry asset named 'AOI': the minimum
   (permanent) reservoir extent polygon.
   ===================================================================== */

//============================ PARAMETERS ==============================
var START = '2000-01-01';
var END   = '2026-01-01';   // exclusive

var FAI_THR       = -0.001; // bloom / non-bloom threshold
var CLOUD_B01_MAX = 1500;   // cloud mask (b01, red 620-670 nm; raw value)
var LAND_B07_MAX  = 1800;   // land mask (b07, SWIR3 2105-2155 nm; raw value)

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

// combine() pairs images by system:index (the acquisition date in both
// products). Check that both collections have the same size before
// relying on the pairing.
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
  // cloud: high red reflectance; land: high SWIR3 reflectance
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

// Visual check
Map.addLayer(IC_FAI.select('FAI'),
  {min: FAI_THR, max: 0.2,
   palette: ['171cb1', '9ef9fd', '41ce13', 'eefd0e', 'f33214', 'fd2b0a', 'b41f07']},
  'FAI', false);

//========================= PER-IMAGE METRICS ==========================
var countAndMedian = ee.Reducer.count()
  .combine({reducer2: ee.Reducer.median(), sharedInputs: true});

function dailyMetrics(fai) {
  fai = ee.Image(fai);

  // over all valid water cells
  var statsAll = fai.reduceRegion({
    reducer: countAndMedian, geometry: aoi,
    scale: SCALE, crs: CRS, maxPixels: 1e10
  });

  // over bloom cells only
  var faiBloom = fai.updateMask(fai.gte(FAI_THR));
  var statsBloom = faiBloom.reduceRegion({
    reducer: countAndMedian, geometry: aoi,
    scale: SCALE, crs: CRS, maxPixels: 1e10
  });

  var nValid = ee.Number(ee.Algorithms.If(statsAll.get('FAI_count'),
                                          statsAll.get('FAI_count'), 0));
  var nBloom = ee.Number(ee.Algorithms.If(statsBloom.get('FAI_count'),
                                          statsBloom.get('FAI_count'), 0));

  // BB_bloom is null on dates with no bloom cell; it is treated as
  // missing (not zero) downstream.
  return ee.Feature(null, {
    date:        ee.Date(fai.get('system:time_start')).format('YYYY-MM-dd'),
    n_valid:     nValid,
    n_bloom:     nBloom,
    RBA:         nBloom.divide(nValid.max(1)),
    BB_bloom:    statsBloom.get('FAI_median'),
    BB_allwater: statsAll.get('FAI_median')
  });
}

var series = ee.FeatureCollection(IC_FAI.map(dailyMetrics))
  .filter(ee.Filter.gt('n_valid', 0));

print('Dates with at least one valid cell:', series.size());
print('Example:', series.first());

//=============================== EXPORT ===============================
Export.table.toDrive({
  collection:     series.select(['date', 'n_valid', 'n_bloom', 'RBA',
                                 'BB_bloom', 'BB_allwater']),
  description:    'MODIS_FAI_daily_metrics',
  fileNamePrefix: 'MODIS_FAI_daily_metrics',
  fileFormat:     'CSV'
});
