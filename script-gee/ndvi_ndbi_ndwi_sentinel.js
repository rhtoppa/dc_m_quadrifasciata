/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var table = ee.FeatureCollection("projects/ee-toppa/assets/limite_8500");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
// 1. Define your area of ​​interest (this can be imported as a shapefile from Assets).
var aoi = table;

// 2. Function to add indexes to the image
function addIndices(image) {
  var ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var ndwi = image.normalizedDifference(['B8', 'B11']).rename('NDWI');
  var ndbi = image.normalizedDifference(['B11', 'B8']).rename('NDBI');
  return image.addBands(ndvi).addBands(ndwi).addBands(ndbi);
}

// 3. Function to filter, calculate indexes, and generate composite data by period.
function processPeriod(startDate, endDate, periodName) {
  var collection = ee.ImageCollection('COPERNICUS/S2_SR')
    .filterBounds(aoi)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
    .map(addIndices);

  var composite = collection.median().clip(aoi);

  // Export NDVI
  Export.image.toDrive({
    image: composite.select('NDVI'),
    description: 'NDVI_' + periodName,
    folder: 'GEE_exports',
    fileNamePrefix: 'NDVI_' + periodName,
    region: aoi.geometry(),
    scale: 10,
    crs: 'EPSG:31983', // SIRGAS 2000 / UTM zone 23S
    maxPixels: 1e13
  });

  // Export NDWI
  Export.image.toDrive({
    image: composite.select('NDWI'),
    description: 'NDWI_' + periodName,
    folder: 'GEE_exports',
    fileNamePrefix: 'NDWI_' + periodName,
    region: aoi.geometry(),
    scale: 10,
    crs: 'EPSG:31983',
    maxPixels: 1e13
  });

  // Export NDBI
  Export.image.toDrive({
    image: composite.select('NDBI'),
    description: 'NDBI_' + periodName,
    folder: 'GEE_exports',
    fileNamePrefix: 'NDBI_' + periodName,
    region: aoi.geometry(),
    scale: 10,
    crs: 'EPSG:31983',
    maxPixels: 1e13
  });
}

// 4. Run for each of the four release periods.
processPeriod('2023-10-09', '2024-01-24', 'Period1');
processPeriod('2024-03-04', '2024-05-03', 'Period2');
processPeriod('2024-05-27', '2024-09-28', 'Period3');
processPeriod('2024-10-07', '2025-01-12', 'Period4');
