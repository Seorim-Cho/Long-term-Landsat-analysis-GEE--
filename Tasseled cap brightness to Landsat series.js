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
Map.setOptions("TERRAIN");

var Piryu_point    = ee.Geometry.Point(126.04023893474744,39.24102786430433);
var Nam_point      = ee.Geometry.Point(125.89172802656319,39.033722152341696);
var Hapjang_point  = ee.Geometry.Point(125.79830917494823,39.05523089830885);
var Potong_point   = ee.Geometry.Point(125.70177180068879,39.00556592964109);
var Sunhwa_point   = ee.Geometry.Point(125.64195538937057,38.99606860054983);
var Konyang_point  = ee.Geometry.Point(125.6002535008059,38.79597867277304);
var Hwangju_point  = ee.Geometry.Point(125.63616979507148,38.71280508163958);
var Chaeryou_point = ee.Geometry.Point(125.64176351769854,38.635680384467925);


var p1 = ee.Geometry.Point(125.9833695502994, 39.16637098866668);
var p2 = ee.Geometry.Point(125.86521687265774,39.038861558863445);
var p3 = ee.Geometry.Point(125.7545391680441,39.00567607744487);
var p4 = ee.Geometry.Point(125.63792063366898,38.970542388267845);
var p5 = ee.Geometry.Point(125.5814210271079,38.886198950185985);
var p6 = ee.Geometry.Point(125.5837825920421,38.75687976045603);
var p7 = ee.Geometry.Point(125.51168481372179,38.714968299848024);
var p8 = ee.Geometry.Point(125.22791731743744,38.687874430041376);


var path = 117;
var row =  33;
var calstart = 3;
var calend = 11;
function maskL457sr(image) {

  var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11111', 2)).eq(0);
  var saturationMask = image.select('QA_RADSAT').eq(0);
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  var thermalBand = image.select('ST_B6').multiply(0.00341802).add(149.0);
 
  return image.addBands(opticalBands, null, true)
      .addBands(thermalBand, null, true)
      .updateMask(qaMask)
      .updateMask(saturationMask);
}
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
//Select Landsat 5 images
var L5SR = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2");
var L5 = ee.ImageCollection(L5SR
  .filter(ee.Filter.eq('WRS_PATH', path))
  .filter(ee.Filter.eq('WRS_ROW', row))  
  .filterDate('1984-01-01', '1998-12-31')
  .filter(ee.Filter.calendarRange(calstart, calend, 'month'))
  .map(maskL457sr)
);
function renameBandsETM(image) {
    var bands =['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'QA_PIXEL'];
    var new_bands = ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2', 'QA_PIXEL'];
    return image.select(bands).rename(new_bands);
}
var L5reName = L5.map(renameBandsETM);
var L7SR = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2");
var L7 = ee.ImageCollection(L7SR
  .filter(ee.Filter.eq('WRS_PATH', path))
  .filter(ee.Filter.eq('WRS_ROW', row))
  .filterDate('1999-01-01', '2013-12-31') 
  .filter(ee.Filter.calendarRange(calstart, calend, 'month'))
  .map(maskL457sr)
);
function renameBandsETM(image) {
    var bands =['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'QA_PIXEL'];
    var new_bands = ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2', 'QA_PIXEL'];
    return image.select(bands).rename(new_bands);
}
var L7reName = L7.map(renameBandsETM);

// Select Landsat 8 images
var L8SR = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2");
var L8 = ee.ImageCollection(L8SR
  .filter(ee.Filter.eq('WRS_PATH', path))
  .filter(ee.Filter.eq('WRS_ROW', row))
  .filterDate('2014-01-01', '2021-12-31')
  .filter(ee.Filter.calendarRange(calstart, calend, 'month'))
  .map(maskL8sr)
);  
function renameBandsOLI(image) {
    var bands = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'QA_PIXEL'];
    var new_bands = ['UltraBlue', 'Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2', 'QA_PIXEL'];
    return image.select(bands).rename(new_bands);
}
var L8reName = L8.map(renameBandsOLI);

// Merge Landat 5, 7, 8 
var ls = L5reName.merge(L7reName).merge(L8reName);

// TCB-Brightness(All landsat sensor, DeVries et al., 2016)
var addTCB = function(image) {
  var TCB = image.expression(
     '0.2043*Blue + 0.4158*Green + 0.5524*Red + 0.5741*NIR + 0.3124*SWIR1 + 0.2303*SWIR2', {
     'Blue' : image.select(['Blue']), 
     'Green' : image.select(['Green']),
     'Red' : image.select(['Red']), 
     'NIR' : image.select(['NIR']), 
     'SWIR1' : image.select(['SWIR1']),
     'SWIR2' : image.select(['SWIR2'])
     } 
    ).rename('TCB');
  return image.addBands(TCB); 
};
var ls_addTCB = ee.ImageCollection(ls.map(addTCB));
var landsat_TC = ls_addTCB.select(['Blue', 'Green', 'Red', 'NIR', 'SWIR1' , 'SWIR2', 'TCB']);
/////////////////////////////////////////////////////////////////////////////////////////////////
// 1. Piryu River
var Piryu_extr = function(image) {
  return image.sampleRegions({collection: Piryu_point, scale: 30, geometries: true});
};
var Piryu_TC = landsat_TC.map(Piryu_extr).flatten();
print(Piryu_TC, 'Piryu_TC');
Export.table.toDrive({
      collection: Piryu_TC,
      description: 'Piryu_TC',
      folder: 'your_folder',
      fileFormat: 'CSV',
    });
// 2. Sunhwa River
var Sunhwa_extr = function(image) {
  return image.sampleRegions({collection: Sunhwa_point, scale: 30, geometries: true});
};
var Sunhwa_TC = landsat_TC.map(Sunhwa_extr).flatten();
print(Sunhwa_TC, 'Sunhwa_TC');
Export.table.toDrive({
      collection: Sunhwa_TC,
      description: 'Sunhwa_TC',
      folder: 'your_folder',
      fileFormat: 'CSV',
    });
// 3. Potong River
var Potong_extr = function(image) {
  return image.sampleRegions({collection: Potong_point, scale: 30, geometries: true});
};
var Potong_TC = landsat_TC.map(Potong_extr).flatten();
print(Potong_TC, 'Potong_TC');
Export.table.toDrive({
      collection: Potong_TC,
      description: 'Potong_TC',
      folder: 'your_folder',
      fileFormat: 'CSV',
    });
// 4. Hapjang River
var Hapjang_extr = function(image) {
  return image.sampleRegions({collection: Hapjang_point, scale: 30, geometries: true});
};
var Hapjang_TC = landsat_TC.map(Hapjang_extr).flatten();
print(Hapjang_TC, 'Hapjang_TC');
Export.table.toDrive({
      collection: Hapjang_TC,
      description: 'Hapjang_TC',
      folder: 'your_folder',
      fileFormat: 'CSV',
    });
// 5. Chaeryoug River
var Chaeryoug_extr = function(image) {
  return image.sampleRegions({collection: Chaeryoung_point, scale: 30, geometries: true});
};
var Chaeryoug_TC = landsat_TC.map(Chaeryoug_extr).flatten();
print(Chaeryoug_TC, 'Chaeryoug_TC');
Export.table.toDrive({
      collection: Chaeryoug_TC,
      description: 'Chaeryoug_TC',
      folder: 'your_folder',
      fileFormat: 'CSV',
    });
// 6. Hwangju River
var Hwangju_extr = function(image) {
  return image.sampleRegions({collection: Hwangju_point, scale: 30, geometries: true});
};
var Hwangju_TC = landsat_TC.map(Hwangju_extr).flatten();
print(Hwangju_TC, 'Hwangju_TC');
Export.table.toDrive({
      collection: Hwangju_TC,
      description: 'Hwangju_TC',
      folder: 'your_folder',
      fileFormat: 'CSV',
    });
// 7. Konyang River
var Konyang_extr = function(image) {
  return image.sampleRegions({collection: Konyang_point, scale: 30, geometries: true});
};
var Konyang_TC = landsat_TC.map(Konyang_extr).flatten();
print(Konyang_TC, 'Konyang_TC');
Export.table.toDrive({
      collection: Konyang_TC,
      description: 'Konyang_TC',
      folder: 'your_folder',
      fileFormat: 'CSV',
    });
    
    
// 8. Nam River
var Nam_extr = function(image) {
  return image.sampleRegions({collection: Nam_point, scale: 30, geometries: true});
};
var Nam_TC = landsat_TC.map(Nam_extr).flatten();
print(Nam_TC, 'Nam_TC');
Export.table.toDrive({
      collection: Nam_TC,
      description: 'Nam_TC',
      folder: 'your_folder',
      fileFormat: 'CSV',
    });
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
// p1
var p1_extr = function(image) {
  return image.sampleRegions({collection: p1, scale: 30, geometries: true});
};
var p1_TC = landsat_TC.map(p1_extr).flatten();
Export.table.toDrive({
      collection: p1_TC,
      description: 'p1_TC',
      folder: 'your_folder',
      fileFormat: 'CSV',
    });
print(p1_TC, 'p1_TC');
// p2
var p2_extr = function(image) {
  return image.sampleRegions({collection: p2, scale: 30, geometries: true});
};
var p2_TC = landsat_TC.map(p2_extr).flatten();
print(p2_TC, 'p2_TC');
Export.table.toDrive({
      collection: p2_TC,
      description: 'p2_TC',
      folder: 'your_folder',
      fileFormat: 'CSV',
    });
// p3
var p3_extr = function(image) {
  return image.sampleRegions({collection: p3, scale: 30, geometries: true});
};
var p3_TC = landsat_TC.map(p3_extr).flatten();
print(p3_TC, 'p3_TC');
Export.table.toDrive({
      collection: p3_TC,
      description: 'p3_TC',
      folder: 'your_folder',
      fileFormat: 'CSV',
    });
    
// p4
var p4_extr = function(image) {
  return image.sampleRegions({collection: p4, scale: 30, geometries: true});
};
var p4_TC = landsat_TC.map(p4_extr).flatten();
print(p4_TC, 'p4_TC');
Export.table.toDrive({
      collection: p4_TC,
      description: 'p4_TC',
      folder: 'your_folder',
      fileFormat: 'CSV',
    });
// p5
var p5_extr = function(image) {
  return image.sampleRegions({collection: p5, scale: 30, geometries: true});
};
var p5_TC = landsat_TC.map(p5_extr).flatten();
print(p5_TC, 'p5_TC');
Export.table.toDrive({
      collection: p5_TC,
      description: 'p5_TC',
      folder: 'your_folder',
      fileFormat: 'CSV',
    });
    
// p6
var p6_extr = function(image) {
  return image.sampleRegions({collection: p6, scale: 30, geometries: true});
};
var p6_TC = landsat_TC.map(p6_extr).flatten();
print(p6_TC, 'p6_TC');
Export.table.toDrive({
      collection: p6_TC,
      description: 'p6_TC',
      folder: 'your_folder',
      fileFormat: 'CSV',
    });
    
// p7
var p7_extr = function(image) {
  return image.sampleRegions({collection: p7, scale: 30, geometries: true});
};
var p7_TC = landsat_TC.map(p7_extr).flatten();
print(p7_TC, 'p7_TC');
Export.table.toDrive({
      collection: p7_TC,
      description: 'p7_TC',
      folder: 'your_folder',
      fileFormat: 'CSV',
    });
// p8
var p8_extr = function(image) {
  return image.sampleRegions({collection: p8, scale: 30, geometries: true});
};
var p8_TC = landsat_TC.map(p8_extr).flatten();
print(p8_TC, 'p8_TC');
Export.table.toDrive({
      collection: p8_TC,
      description: 'p8_TC',
      folder: 'your_folder',
      fileFormat: 'CSV',
    });
