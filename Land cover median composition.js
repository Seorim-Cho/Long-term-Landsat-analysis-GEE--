/*
Author: Seorim Cho (srcho@khu.ac.kr, chosrm15@gmail.com )

This code is free and open. 
By using this code and any data derived with it, 
you agree to cite the following reference 
in any publications derived from them: ""

This function applies a Tasseled cap transformation developed by DeVires et al. (2016) to Landsat series imagery.
in order to obtain a Tasseled cap brightness value at each pixel
: DeVries, B., Pratihast, A. K., Verbesselt, J., Kooistra, L., & Herold, M. (2016). Characterizing forest change using community-based monitoring data and Landsat time series. PloS one, 11(3), e0147121.

Also, you can access to Google Earth Engine (GEE) code
- GEE Git: https://earthengine.googlesource.com/users/SEORIM/Taedong_River
- GEE repository: https://code.earthengine.google.com/?accept_repo=users/SEORIM/Taedong_River
*/

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


///////////////////////////////////////////////////////////////////////////////////
// Compute  Tasseled Cap Trasformation
var L8_col = L8pr1.merge(L8pr2).merge(L8pr3).merge(L8pr4);

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
