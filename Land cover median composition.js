Map.setCenter(125.72239942629965,38.990859676683726);
Map.setOptions("SATELLITE");

var boundary= ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017");
var nkBorder = boundary.filter(ee.Filter.eq('country_co', 'KN'));

var watershed = ee.FeatureCollection(Watershed);

var styling = { color: 'black', fillColor: '#ffffff00' };
Map.addLayer(watershed.style(styling), {}, 'Taedong Watershed');

var start = '2014-01-01'; 
var end = '2019-12-31';
var calstart = 5; 
var calend = 9;
var cloudcover = 1;

function maskL8sr(image) {

  var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11111', 2)).eq(0);
  var saturationMask = image.select('QA_RADSAT').eq(0);
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);

  return image.addBands(opticalBands, null, true)
      .addBands(thermalBands, null, true)
      .updateMask(qaMask)
      .updateMask(saturationMask);
}

var L8SR = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2");

// 1) 117/33
var L8pr1 = ee.ImageCollection(L8SR
  .filter(ee.Filter.eq('WRS_PATH', 117))
  .filter(ee.Filter.eq('WRS_ROW', 33))
  .filterDate(start, end)
   .filter(ee.Filter.calendarRange(calstart, calend, 'month'))
  .filterMetadata('CLOUD_COVER', 'less_than', cloudcover)
  .map(maskL8sr)
  .sort('system:time_start')
); 
// 2) 116/33
var L8pr2 = ee.ImageCollection(L8SR
  .filter(ee.Filter.eq('WRS_PATH', 116))
  .filter(ee.Filter.eq('WRS_ROW', 33))
  .filterDate(start, end)
   .filter(ee.Filter.calendarRange(calstart, calend, 'month'))
  .filterMetadata('CLOUD_COVER', 'less_than', cloudcover)
  .map(maskL8sr)
  .sort('system:time_start')
);
// 3) 116/32
var L8pr3 = ee.ImageCollection(L8SR
  .filter(ee.Filter.eq('WRS_PATH', 116))
  .filter(ee.Filter.eq('WRS_ROW', 32))
  .filterDate(start, end)
   .filter(ee.Filter.calendarRange(calstart, calend, 'month'))
  .filterMetadata('CLOUD_COVER', 'less_than', cloudcover)
  .map(maskL8sr)
  .sort('system:time_start')
);
// 4) 117/32
var L8pr4 = ee.ImageCollection(L8SR
  .filter(ee.Filter.eq('WRS_PATH', 117))
  .filter(ee.Filter.eq('WRS_ROW', 32))
  .filterDate(start, end)
   .filter(ee.Filter.calendarRange(calstart, calend, 'month'))
  .filterMetadata('CLOUD_COVER', 'less_than', cloudcover)
  .map(maskL8sr)
  .sort('system:time_start')
);

// Landsat Collection
var landsatCollection = L8pr1.merge(L8pr2).merge(L8pr3).merge(L8pr4).sort('system:time_start');


// path/row collection's Median
var median_pr1 = L8pr1.median();
Map.addLayer(median_pr1, l8visualization, 'Median_pr1');

var median_pr2 = L8pr2.median();
Map.addLayer(median_pr2, l8visualization, 'Median_pr2'); 

var median_pr3 = L8pr3.median();
Map.addLayer(median_pr3, l8visualization, 'Median_pr3'); 

var median_pr4 = L8pr4.median();
Map.addLayer(median_pr4, l8visualization, 'Median_pr4'); 

// Median Image Collection
var median_col = ee.ImageCollection([median_pr1, median_pr2, median_pr3, median_pr4]);
var median_clip =median_col.median().clip(watershed);

    
// Export each bands    
Export.image.toDrive({ 
      image: median_clip.select('SR_B2'),
      description: 'SR_B2',
      folder:'your_folder',
      scale: 30,
      fileFormat: 'GeoTIFF',
      maxPixels: 57614177730,
      region: watershed
    });


Export.image.toDrive({ 
      image: median_clip.select('SR_B3'),
      description: 'SR_B3',
      folder:'your_folder',
      scale: 30,
      fileFormat: 'GeoTIFF',
      maxPixels: 57614177730,
      region: watershed
    });

Export.image.toDrive({ 
      image: median_clip.select('SR_B4'),
      description: 'SR_B4',
      folder:'your_folder',
      scale: 30,
      fileFormat: 'GeoTIFF',
      maxPixels: 57614177730,
      region: watershed
    });


Export.image.toDrive({ 
      image: median_clip.select('SR_B5'),
      description: 'SR_B5',
      folder:'your_folder',
      scale: 30,
      fileFormat: 'GeoTIFF',
      maxPixels: 57614177730,
      region: watershed
    });

Export.image.toDrive({ 
      image: median_clip.select('SR_B6'),
      description: 'SR_B6',
      folder:'your_folder',
      scale: 30,
      fileFormat: 'GeoTIFF',
      maxPixels: 57614177730,
      region: watershed
    });


Export.image.toDrive({ 
      image: median_clip.select('SR_B7'),
      description: 'SR_B7',
      folder:'your_folder',
      scale: 30,
      fileFormat: 'GeoTIFF',
      maxPixels: 57614177730,
      region: watershed
    });




// Compute NDWI and Tasseled Cap Trasformation
var L8_col = L8pr1.merge(L8pr2).merge(L8pr3).merge(L8pr4);


// NDWI
// The Normalized Difference Water Index (NDWI) is derived from the Near-Infrared (NIR) and Green (G) channels.
var addNDWI = function(image) {
  var ndwi = image.normalizedDifference(['SR_B3', 'SR_B5']).rename('NDWI');
  return image.addBands(ndwi);
};

var NDWI = L8_col.map(addNDWI).select('NDWI');

Map.addLayer(NDWI.mean().clip(watershed), {min: -1, max: 1, palette: ['white', 'skyblue', 'blue']}, 'Mean NDWIpr2', false); // mean


var NDWI_mean = NDWI.mean().clip(watershed);
Export.image.toDrive({ 
      image: NDWI_mean,
      description: 'NDWI_mean',
      folder:'your_folder',
      scale: 30,
      fileFormat: 'GeoTIFF',
      maxPixels: 57614177730,
      region: watershed
    });

    

// Tasseled cap brightness (DeVries et al,, 2016)
var addTCB = function(image) {
  var TCB = image.expression(
     '0.2043*Blue + 0.4158*Green + 0.5524*Red + 0.5741*NIR + 0.3124*SWIR1 + 0.2303*SWIR2', {
     'Blue' : image.select(['SR_B2']), 
     'Green' : image.select(['SR_B3']),
     'Red' : image.select(['SR_B4']), 
     'NIR' : image.select(['SR_B5']), 
     'SWIR1' : image.select(['SR_B6']),
     'SWIR2' : image.select(['SR_B7'])
     } 
    ).rename('TCB');
  return image.addBands(TCB);
};

var TCB = L8_col.map(addTCB).select('TCB');
var TCB_mean = TCB.mean().clip(watershed);
Export.image.toDrive({ 
      image: TCB_mean,
      description: 'TCB_mean',
      folder:'your_folder',
      scale: 30,
      fileFormat: 'GeoTIFF',
      maxPixels: 57614177730,
      region: watershed
    });
    
    
// Tasseled cap greenness (DeVries et al,, 2016)
var addTCG = function(image) {
  var TCG = image.expression(
     '(-0.1603)*Blue + 0.2819*Green + (-0.4934)*Red + 0.7940*NIR + (-0.0002)*SWIR1 + (-0.1446)*SWIR2', {
     'Blue' : image.select(['SR_B2']), 
     'Green' : image.select(['SR_B3']),
     'Red' : image.select(['SR_B4']), 
     'NIR' : image.select(['SR_B5']), 
     'SWIR1' : image.select(['SR_B6']),
     'SWIR2' : image.select(['SR_B7'])
     } 
    ).rename('TCG');
  return image.addBands(TCG);
};

var TCG = L8_col.map(addTCG).select('TCG');
var TCG_mean = TCG.mean().clip(watershed);
Export.image.toDrive({ 
      image: TCG_mean,
      description: 'TCG_mean',
      folder:'your_folder',
      scale: 30,
      fileFormat: 'GeoTIFF',
      maxPixels: 57614177730,
      region: watershed
    });

// Tasseled cap wetness(DeVries et al,, 2016)
var addTCW = function(image) {
  var TCW = image.expression(
     '0.0315*Blue + 0.2021*Green + 0.3102*Red + 0.1594*NIR + (-0.6806)*SWIR1 + (-0.6109)*SWIR2', {
     'Blue' : image.select(['SR_B2']), 
     'Green' : image.select(['SR_B3']),
     'Red' : image.select(['SR_B4']), 
     'NIR' : image.select(['SR_B5']), 
     'SWIR1' : image.select(['SR_B6']),
     'SWIR2' : image.select(['SR_B7'])
     } 
    ).rename('TCW');
  return image.addBands(TCW); 
};

var TCW = L8_col.map(addTCW).select('TCW');
var TCW_mean = TCW.mean().clip(watershed);
Export.image.toDrive({ 
      image: TCW_mean,
      description: 'TCW_mean',
      folder:'your_folder',
      scale: 30,
      fileFormat: 'GeoTIFF',
      maxPixels: 57614177730,
      region: watershed
    });


// Display Tributaries
Map.addLayer(Piryu.style({
  fillColor: '#ffffff00',
  color: 'red',
  width: 2.0
}), {}, 'Piryu', false); 

Map.addLayer(Hwangju.style({
  fillColor: '#ffffff00',
  color: 'red',
  width: 2.0
}), {}, 'Hwangju', false); 

Map.addLayer(Chaeryoung.style({
  fillColor: '#ffffff00',
  color: 'red',
  width: 2.0
}), {}, 'Chaeryoung', false); 

Map.addLayer(Sunhwa.style({
  fillColor: '#ffffff00',
  color: 'red',
  width: 2.0
}), {}, 'Sunhwa', false); 

Map.addLayer(Nam.style({
  fillColor: '#ffffff00',
  color: 'red',
  width: 2.0
}), {}, 'Nam', false); 

Map.addLayer(Hapjang.style({
  fillColor: '#ffffff00',
  color: 'red',
  width: 2.0
}), {}, 'Hapjang', false); 

Map.addLayer(Potong.style({
  fillColor: '#ffffff00',
  color: 'red',
  width: 2.0
}), {}, 'Potong', false); 


Map.addLayer(Konyang.style({
  fillColor: '#ffffff00',
  color: 'red',
  width: 2.0
}), {}, 'Konyang', false);
