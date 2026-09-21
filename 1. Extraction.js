// ============================================================
// NAGPUR DISTRICT - MULTI-INDEX ANALYSIS
// Includes: LST, NDBI, NDVI, NDWI, LULC, Albedo, Terrain
// Summer: March 1 – September 1
// ============================================================

// ------------------------------------------------------------
// 1. LOAD NAGPUR DISTRICT
// ------------------------------------------------------------
var nagpur = ee.FeatureCollection('projects/vidarbha-lst-project/assets/Nagpur_District');
var region = nagpur.geometry();
Map.centerObject(nagpur, 9);
Map.addLayer(nagpur, {color: 'red'}, 'Nagpur District');

// ------------------------------------------------------------
// 2. COMMON SETTINGS
// ------------------------------------------------------------
var START_YEAR = 2017;
var END_YEAR = 2025;
var START_MONTH = 3;
var START_DAY = 1;
var END_MONTH = 9;
var END_DAY = 1;               // September 1 is exclusive
var MAX_CLOUD_PROBABILITY = 40;

// ============================================================
// 3. LANDSAT LST (LANDSAT 8 & 9)
// ============================================================

// Cloud mask for Landsat Level-2
function maskLandsat(image) {
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 1).eq(0)
               .and(qa.bitwiseAnd(1 << 2).eq(0))
               .and(qa.bitwiseAnd(1 << 3).eq(0))
               .and(qa.bitwiseAnd(1 << 4).eq(0))
               .and(qa.bitwiseAnd(1 << 5).eq(0));
  var saturation = image.select('QA_RADSAT').eq(0);
  var lst = image.select('ST_B10')
                .multiply(0.00341802)
                .add(149.0)
                .subtract(273.15)
                .rename('LST');
  var validLST = lst.gt(-20).and(lst.lt(70));
  return lst
    .updateMask(mask)
    .updateMask(saturation)
    .updateMask(validLST)
    .copyProperties(image, ['system:time_start']);
}

// Get Landsat collection for a given year (summer)
function getLandsat(year) {
  var startDate = ee.Date.fromYMD(year, START_MONTH, START_DAY);
  var endDate = ee.Date.fromYMD(year, END_MONTH, END_DAY);
  var landsat8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2');
  var landsat9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2');
  var collection = landsat8
    .merge(landsat9)
    .filterBounds(region)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.lt('CLOUD_COVER', 70))
    .map(maskLandsat);
  return collection;
}

// Create median LST composite for a year
function createLST(year) {
  var collection = getLandsat(year);
  return collection.median().clip(region).rename('LST');
}

// Print image counts and export LST composites
print('==============================');
print('LST IMAGE COUNTS');
print('==============================');
for (var y = START_YEAR; y <= END_YEAR; y++) {
  var coll = getLandsat(y);
  print('Summer ' + y + ' image count:', coll.size());
}

for (var year = START_YEAR; year <= END_YEAR; year++) {
  var lst = createLST(year);
  // Display first year only
  if (year === START_YEAR) {
    Map.addLayer(lst, {
      min: 20,
      max: 55,
      palette: ['040274','2c7bb6','00a6ca','00ccbc','90eb9d','ffff8c','f9d057','f29e2e','e76818','d7191c']
    }, 'LST Summer ' + year);
  }
  Export.image.toDrive({
    image: lst,
    description: 'LST_' + year + '_Summer_Nagpur',
    folder: 'Nagpur_LST_2017_2025',
    fileNamePrefix: 'LST_' + year + '_Summer_Nagpur',
    region: region,
    scale: 30,
    crs: 'EPSG:4326',
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF',
    formatOptions: {cloudOptimized: true}
  });
}

// ============================================================
// 4. SENTINEL-2 INDICES (NDBI, NDVI, NDWI)
// ============================================================

// Mask edges (B8A and B9 bands)
function maskEdges(image) {
  var edgeMask1 = image.select('B8A').mask();
  var edgeMask2 = image.select('B9').mask();
  return image.updateMask(edgeMask1).updateMask(edgeMask2);
}

// Cloud probability mask
function maskClouds(image) {
  var cloudProbability = ee.Image(image.get('cloud_mask')).select('probability');
  var cloudMask = cloudProbability.lt(MAX_CLOUD_PROBABILITY);
  return image.updateMask(cloudMask);
}

// Get Sentinel-2 SR collection with cloud probability masking
function getSentinel2(year) {
  var startDate = ee.Date.fromYMD(year, START_MONTH, START_DAY);
  var endDate = ee.Date.fromYMD(year, END_MONTH, END_DAY);
  var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(region)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 80))
    .map(maskEdges);
  var clouds = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
    .filterBounds(region)
    .filterDate(startDate, endDate);
  var joined = ee.Join.saveFirst('cloud_mask').apply({
    primary: s2,
    secondary: clouds,
    condition: ee.Filter.equals({leftField: 'system:index', rightField: 'system:index'})
  });
  var masked = ee.ImageCollection(joined).map(maskClouds);
  return masked;
}

// Create NDBI (B11 – B8) / (B11 + B8)
function createNDBI(year) {
  var collection = getSentinel2(year);
  var ndbiCollection = collection.map(function(image) {
    var ndbi = image.normalizedDifference(['B11', 'B8']).rename('NDBI');
    return ndbi.copyProperties(image, ['system:time_start']);
  });
  return ndbiCollection.median().clip(region).rename('NDBI').toFloat();
}

// Create NDVI (B8 – B4) / (B8 + B4)
function createNDVI(year) {
  var collection = getSentinel2(year);
  var ndviCollection = collection.map(function(image) {
    var ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI');
    return ndvi.copyProperties(image, ['system:time_start']);
  });
  return ndviCollection.median().clip(region).rename('NDVI').toFloat();
}

// Create NDWI (B3 – B8) / (B3 + B8)
function createNDWI(year) {
  var collection = getSentinel2(year);
  var ndwiCollection = collection.map(function(image) {
    var ndwi = image.normalizedDifference(['B3', 'B8']).rename('NDWI');
    return ndwi.copyProperties(image, ['system:time_start']);
  });
  return ndwiCollection.median().clip(region).rename('NDWI').toFloat();
}

// Print Sentinel‑2 image counts
print('==============================');
print('SENTINEL-2 IMAGE COUNTS');
print('==============================');
for (var year = START_YEAR; year <= END_YEAR; year++) {
  var coll = getSentinel2(year);
  print('Summer ' + year + ' images:', coll.size());
}

// Display and export NDBI, NDVI, NDWI for all years
var ndbi2017 = createNDBI(START_YEAR);
Map.addLayer(ndbi2017, {
  min: -0.5,
  max: 0.5,
  palette: ['006400','7FFF00','FFFF00','FFA500','FF0000','800000']
}, 'NDBI Summer ' + START_YEAR);

var ndvi2017 = createNDVI(START_YEAR);
Map.addLayer(ndvi2017, {
  min: -0.2,
  max: 0.8,
  palette: ['8c510a','d8b365','f6e8c3','c7eae5','5ab4ac','01665e']
}, 'NDVI Summer ' + START_YEAR);

var ndwi2017 = createNDWI(START_YEAR);
Map.addLayer(ndwi2017, {
  min: -0.5,
  max: 0.5,
  palette: ['8c510a','d8b365','f6e8c3','c7eae5','5ab4ac','01665e']
}, 'NDWI Summer ' + START_YEAR);

// Export NDBI, NDVI, NDWI for each year
for (var year = START_YEAR; year <= END_YEAR; year++) {
  var ndbi = createNDBI(year);
  Export.image.toDrive({
    image: ndbi,
    description: 'NDBI_' + year + '_Summer_Nagpur',
    folder: 'Nagpur_NDBI_2017_2025',
    fileNamePrefix: 'NDBI_' + year + '_Summer_Nagpur',
    region: region,
    scale: 30,
    crs: 'EPSG:4326',
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF',
    formatOptions: {cloudOptimized: true}
  });

  var ndvi = createNDVI(year);
  Export.image.toDrive({
    image: ndvi,
    description: 'NDVI_' + year + '_Summer_Nagpur',
    folder: 'Nagpur_NDVI_2017_2025',
    fileNamePrefix: 'NDVI_' + year + '_Summer_Nagpur',
    region: region,
    scale: 30,
    crs: 'EPSG:4326',
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF',
    formatOptions: {cloudOptimized: true}
  });

  var ndwi = createNDWI(year);
  Export.image.toDrive({
    image: ndwi,
    description: 'NDWI_' + year + '_Summer_Nagpur',
    folder: 'Nagpur_NDWI_2017_2025',
    fileNamePrefix: 'NDWI_' + year + '_Summer_Nagpur',
    region: region,
    scale: 30,
    crs: 'EPSG:4326',
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF',
    formatOptions: {cloudOptimized: true}
  });
}

// ============================================================
// 5. LULC (DYNAMIC WORLD)
// ============================================================

var CLASS_NAMES = ['water','trees','grass','flooded_vegetation','crops','shrub_and_scrub','built','bare','snow_and_ice'];
var LULC_PALETTE = ['419BDF','397D49','88B053','7A87C6','E49635','DFC35A','C4281B','A59B8F','B39FE1'];

function createAnnualLULC(year) {
  var startDate = ee.Date.fromYMD(year, 1, 1);
  var endDate = ee.Date.fromYMD(year + 1, 1, 1);
  var dw = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1')
    .filterBounds(region)
    .filterDate(startDate, endDate);
  print('Dynamic World image count - ' + year, dw.size());
  var probabilities = dw.select(CLASS_NAMES).median();
  var lulc = probabilities.toArray()
    .arrayArgmax()
    .arrayGet([0])
    .rename('LULC')
    .clip(region)
    .toByte();
  return lulc;
}

// Display 2017 for checking
var lulc2017 = createAnnualLULC(START_YEAR);
Map.addLayer(lulc2017, {min: 0, max: 8, palette: LULC_PALETTE}, 'LULC 2017');
print('LULC 2017', lulc2017);

// Export LULC for each year
for (var year = START_YEAR; year <= END_YEAR; year++) {
  var lulc = createAnnualLULC(year);
  Export.image.toDrive({
    image: lulc,
    description: 'LULC_' + year + '_Nagpur',
    folder: 'Nagpur_LULC_Annual_2017_2025',
    fileNamePrefix: 'LULC_' + year + '_Nagpur',
    region: region,
    scale: 30,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF',
    formatOptions: {cloudOptimized: true}
  });
}

// ============================================================
// 6. ALBEDO (MODIS MCD43A3 V6.1)
// ============================================================

function maskAlbedoQuality(image) {
  var qa = image.select('BRDF_Albedo_Band_Mandatory_Quality_shortwave');
  var good = qa.bitwiseAnd(1).eq(0);
  return image
    .updateMask(good)
    .select('Albedo_WSAshortwave')
    .multiply(0.001)
    .copyProperties(image, ['system:time_start']);
}

function createSummerAlbedo(year) {
  var startDate = ee.Date.fromYMD(year, START_MONTH, START_DAY);
  var endDate = ee.Date.fromYMD(year, END_MONTH, END_DAY);
  var albedoColl = ee.ImageCollection('MODIS/061/MCD43A3')
    .filterBounds(region)
    .filterDate(startDate, endDate)
    .map(maskAlbedoQuality);
  print('Summer Albedo images - ' + year, albedoColl.size());
  return albedoColl.median().rename('Albedo').clip(region);
}

// Display 2017
var albedo2017 = createSummerAlbedo(START_YEAR);
Map.addLayer(albedo2017, {
  min: 0,
  max: 0.6,
  palette: ['0000FF','00FFFF','00FF00','FFFF00','FF0000']
}, 'Albedo Summer 2017');

// Export albedo for each year
for (var year = START_YEAR; year <= END_YEAR; year++) {
  var albedo = createSummerAlbedo(year);
  Export.image.toDrive({
    image: albedo,
    description: 'Albedo_Summer_' + year + '_Nagpur',
    folder: 'Nagpur_Albedo_Summer_2017_2025',
    fileNamePrefix: 'Albedo_Summer_' + year + '_Nagpur',
    region: region,
    scale: 30,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF',
    formatOptions: {cloudOptimized: true}
  });
}

// ============================================================
// 7. TERRAIN (SRTM DEM)
// ============================================================

var dem = ee.Image('USGS/SRTMGL1_003').select('elevation').clip(region);
var elevation = dem.rename('Elevation');
var slope = ee.Terrain.slope(dem).rename('Slope');
var aspect = ee.Terrain.aspect(dem).rename('Aspect');

Map.addLayer(elevation, {min: 200, max: 500, palette: ['006400','7FFF00','FFFF00','FFA500','FF0000']}, 'Elevation');
Map.addLayer(slope, {min: 0, max: 30, palette: ['FFFFFF','FFFF00','FFA500','FF0000']}, 'Slope');
Map.addLayer(aspect, {min: 0, max: 360, palette: ['0000FF','00FFFF','00FF00','FFFF00','FF0000','FF00FF','0000FF']}, 'Aspect');

print('DEM:', dem);
print('Elevation:', elevation);
print('Slope:', slope);
print('Aspect:', aspect);

// Export elevation
Export.image.toDrive({
  image: elevation,
  description: 'Elevation_Nagpur',
  folder: 'Nagpur_Terrain',
  fileNamePrefix: 'Elevation_Nagpur',
  region: region,
  scale: 30,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF',
  formatOptions: {cloudOptimized: true}
});

// Export slope
Export.image.toDrive({
  image: slope,
  description: 'Slope_Nagpur',
  folder: 'Nagpur_Terrain',
  fileNamePrefix: 'Slope_Nagpur',
  region: region,
  scale: 30,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF',
  formatOptions: {cloudOptimized: true}
});

// Export aspect
Export.image.toDrive({
  image: aspect,
  description: 'Aspect_Nagpur',
  folder: 'Nagpur_Terrain',
  fileNamePrefix: 'Aspect_Nagpur',
  region: region,
  scale: 30,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF',
  formatOptions: {cloudOptimized: true}
});
