//load copernicus polygons in the GEE asset

//=======================================================================================================================================
// EC-JRC 2021 - Code to create S1 VV, VH and VH/VV 10-days time series to select features and parametrize the 
// classification models for the EU crop map
// Related to the paper: "From parcel to continental scale -- 
// A first European crop type map based on Sentinel-1 and LUCAS Copernicus in-situ observations" 
// by Raphaël d'Andrimont, Astrid Verhegghen, Guido Lemoine, Pieter Kempeneers, Michele Meroni, Marijn van der Velde. 
//=======================================================================================================================================

///////////////////////////
//   INPUTS
///////////////////////////////////////////

// 1. Date
var start_date='2018-01-01';
var end_date='2018-12-31';  

// 2. Time step
var step = 10//12//6//7 // in days (time window for meaning)

// 3. Pixel spacing 10 or 20 m
var pix_export=10

// 4. Buffer polygones lucas in meters - not in the final version
//var buffer=-10

// 5. Polygons
var parcel = lucas_training_all ;
print('parcel',parcel.limit(10))

// 6. Select classtype of interest - Run for each crop type

var classtype='urban'
var parcel = parcel.filter(ee.Filter.inList('LC1_L1',["A"])) ;

//var classtype='crop'
//var parcel = parcel.filter(ee.Filter.inList('LC1_L1',["B"])) ;

//var classtype='woodland'
//var parcel = parcel.filter(ee.Filter.inList('LC1_L1',["C"])) ;

//var classtype='shrubland'
//var parcel = parcel.filter(ee.Filter.inList('LC1_L1',["D"])) ;

//var classtype='grassland'
//var parcel = parcel.filter(ee.Filter.inList('LC1_L1',["E"])) ;

//var classtype='bareland'
//var parcel = parcel.filter(ee.Filter.inList('LC1_L1',["F"])) ;

// var classtype='water'
// var parcel = parcel.filter(ee.Filter.inList('LC1_L1',["G"])) ;

//var classtype='wetlands'
//var parcel = parcel.filter(ee.Filter.inList('LC1_L1',["H"])) ;

print('number of parcels of '+classtype, parcel.size())
Map.addLayer(parcel.draw('red'),{},'lucas_training_raster')


// 7. Regions for the export

var EU_NW1 = ee.Geometry.Rectangle(-13.69,48,0,70.1);
var EU_NW2a = ee.Geometry.Rectangle(0,48,13,50);
var EU_NW2b = ee.Geometry.Rectangle(0,50,13,70.1);
var EU_NE1a = ee.Geometry.Rectangle(13,48,23.5,51);
var EU_NE1b = ee.Geometry.Rectangle(13,51,23.5,56);
var EU_NE1c = ee.Geometry.Rectangle(13,56,23.5,60);
var EU_NE1d = ee.Geometry.Rectangle(13,60,23.5,70.1);
var EU_NE2 = ee.Geometry.Rectangle(23.5,48,34.7,70.1);
var EU_SW1 = ee.Geometry.Rectangle(-13.69,32.63,0,48);
var EU_SW2 = ee.Geometry.Rectangle(0,35.5,13,48);
var EU_SE1 = ee.Geometry.Rectangle(13,32.63,23.5,48);
var EU_SE2 = ee.Geometry.Rectangle(23.5,32.63,34.7,48);


Map.addLayer(EU_NW1,{},'EU_NW1')
Map.addLayer(EU_NW2a,{},'EU_NW2a')
Map.addLayer(EU_NW2b,{},'EU_NW2b')
Map.addLayer(EU_NE1a,{},'EU_NE1a')
Map.addLayer(EU_NE1b,{},'EU_NE1b')
Map.addLayer(EU_NE1c,{},'EU_NE1c')
Map.addLayer(EU_NE1d,{},'EU_NE1d')
Map.addLayer(EU_NE2,{},'EU_NE2')
Map.addLayer(EU_SW1,{},'EU_SW1')
Map.addLayer(EU_SW2,{},'EU_SW2')
Map.addLayer(EU_SE1,{},'EU_SE1')
Map.addLayer(EU_SE2,{},'EU_SE2')


///////////////////////////////////////////
// A / Data preparation
///////////////////////////////////////////

// 1. Parcels

//with the buffer we loose many polygons
/*
var parcel = parcel.map(function(f) { return f.buffer(buffer)})
var parcel = parcel.map(function(f) { return f.set({'bufferedarea': f.area()}) });*/
//Map.addLayer(parcel.draw('blue'),{},'lucas_training_raster_buffered')

// 2. Regions

// additional filter if needed - biome
// var parcel_biom1 = parcel.filterMetadata('BIOME_N','equals',1) ;
// var parcel_biom2 = parcel.filterMetadata('BIOME_N','equals',2) ;
// var parcel_biom3 = parcel.filterMetadata('BIOME_N','equals',3) ;
// var parcel_biom4 = parcel.filterMetadata('BIOME_N','equals',4) ;

var parcel_EU_NW1 = parcel.filterBounds(EU_NW1) ;
var parcel_EU_NW2a = parcel.filterBounds(EU_NW2a) ;
var parcel_EU_NW2b = parcel.filterBounds(EU_NW2b) ;
var parcel_EU_NE1a = parcel.filterBounds(EU_NE1a) ;
var parcel_EU_NE1b = parcel.filterBounds(EU_NE1b) ;
var parcel_EU_NE1c = parcel.filterBounds(EU_NE1c) ;
var parcel_EU_NE1d = parcel.filterBounds(EU_NE1d) ;
var parcel_EU_NE2 = parcel.filterBounds(EU_NE2) ;
var parcel_EU_SW1 = parcel.filterBounds(EU_SW1) ;
var parcel_EU_SW2 = parcel.filterBounds(EU_SW2) ;
var parcel_EU_SE1 = parcel.filterBounds(EU_SE1) ;
var parcel_EU_SE2 = parcel.filterBounds(EU_SE2) ;

Map.addLayer(parcel_EU_SW1.draw('white'),{},'lucas_training_raster_buffered')

print('number of parcels of '+classtype+ ' EU_NW1',parcel_EU_NW1.size())
print('number of parcels of '+classtype+ ' EU_NW2a',parcel_EU_NW2a.size())
print('number of parcels of '+classtype+ ' EU_NW2b',parcel_EU_NW2b.size())
print('number of parcels of '+classtype+ ' EU_NE1a',parcel_EU_NE1a.size())
print('number of parcels of '+classtype+ ' EU_NE1b',parcel_EU_NE1b.size())
print('number of parcels of '+classtype+ ' EU_NE1c',parcel_EU_NE1c.size())
print('number of parcels of '+classtype+ ' EU_NE1d',parcel_EU_NE1d.size())
print('number of parcels of '+classtype+ ' EU_NE2',parcel_EU_NE2.size())
print('number of parcels of '+classtype+ ' EU_SW1',parcel_EU_SW1.size())
print('number of parcels of '+classtype+ ' EU_SW2',parcel_EU_SW2.size())
print('number of parcels of '+classtype+ ' EU_SE1',parcel_EU_SE1.size())
print('number of parcels of '+classtype+ ' EU_SE2',parcel_EU_SE2.size())


///////////////////////////////////////////
// B / FUNCTIONS
///////////////////////////////////////////

// Functions to convert from/to dB
function toNatural(img) {
  return ee.Image(10.0).pow(img.select('..').divide(10.0)).copyProperties(img, ['system:time_start'])
}

function toDBmodif(img) {
  var todB = img.select('VV','VH').log10().multiply(10.0).rename('VV_db','VH_db')
  return img.addBands(todB,['VV_db','VH_db']).copyProperties(img, ['system:time_start']);
}

function toDB(img) {
  return ee.Image(img).log10().multiply(10.0).copyProperties(img, ['system:time_start']);
}

// Remove ugly edges
function maskEdge(img) {
  //connected component is giving a unique label to connected pixels of the same value and we are masking the ones bigger than 100px
  var mask = img.select(0).unitScale(-25, 5).multiply(255).toByte().connectedComponents(ee.Kernel.rectangle(1,1), 100);
 return img.updateMask(mask.select(0).abs());  
}

function S1VHVV (img) { var ratio = img.select(['VH']).divide(img.select(['VV'])).rename('VHVV');
  return img.addBands(ratio,['VHVV'])
} 

// Convert ImageCollection to image stack
function stack(i1, i2)
{
  return ee.Image(i1).addBands(ee.Image(i2))
}


///////////////////////////////////////////
// C / PROCESSING
///////////////////////////////////////////

// get the data from S1 (VV pol.)
var s1 = ee.ImageCollection('COPERNICUS/S1_GRD').filterMetadata('instrumentMode', 'equals', 'IW').
  filter(ee.Filter.eq('transmitterReceiverPolarisation', ['VV', 'VH'])).
  filterBounds(parcel).filterDate(start_date, end_date)
  .sort('system:time');

Map.addLayer(s1,{},'s1')

// Remove ugly edges
s1 = s1.map(maskEdge)


// Ratio and mean are made from natural (non-logarithmic) values 
s1 = s1.map(toNatural)

// Olha Danylo's procedure to create weekly means (adapted)

var days = ee.List.sequence(0, ee.Date(end_date).difference(ee.Date(start_date), 'day'), step).
  map(function(d) { return ee.Date(start_date).advance(d, "day") })

print ('days', days)

var dates = days.slice(0,-1).zip(days.slice(1))

print ('dates', dates)

var s1res = dates.map(function(range) {
  var dstamp = ee.Date(ee.List(range).get(0)).format('YYYYMMdd')
  var temp_collection = s1.filterDate(ee.List(range).get(0),
  ee.List(range).get(1)).mean().select(['VV', 'VH'], [ee.String('VV_').cat(dstamp), ee.String('VH_').cat(dstamp)])
  return temp_collection
})

print('s1res',s1res)

//transform back to DB 
s1res=s1res.map(toDB)

print('s1res',s1res)

//VHVV are computed separately as they should not be converted to dB

s1 = s1.map(S1VHVV)

var s1resRatio = dates.map(function(range) {
  var dstamp = ee.Date(ee.List(range).get(0)).format('YYYYMMdd')
  var temp_collection = s1.filterDate(ee.List(range).get(0),
  ee.List(range).get(1)).mean().select(['VHVV'], [ee.String('VHVV_').cat(dstamp)])
  return temp_collection
})

print('s1resRatio',s1resRatio)

//put the two lists together 
var combine_db_ratio = s1res.zip(s1resRatio).flatten()
print(combine_db_ratio)

// Convert ImageCollection to image stack
var s1stack = combine_db_ratio.slice(1).iterate(stack, combine_db_ratio.get(0))
print(s1stack)

// OPTIONAL : Smooth the image by convolving with the boxcar kernel
// var boxcar = ee.Kernel.square({
//   radius: 3, units: 'pixels', normalize: true
// });
//s1stack = ee.Image(s1stack).convolve(boxcar).clip(aoi)

//transform the image to float to reduce size
s1stack = ee.Image(s1stack).toFloat()

print('s1stack',s1stack)

//Display the stack 
Map.addLayer(s1stack,{},'s1stack_europe')
Map.addLayer(s1stack.clip(parcel),{},'s1stack_europe parcel')

///////////////////////////////////////////
// D / EXPORT
///////////////////////////////////////////

var trainingexportdrive_EU_NW1 = s1stack.sampleRegions({
  collection: parcel_EU_NW1,
  properties: ['POINT_ID','stratum','LC1','LU1'],
  tileScale:16,
  scale: pix_export,
  geometries:false
});

var trainingexportdrive_EU_NW2a = s1stack.sampleRegions({
  collection: parcel_EU_NW2a,
  properties: ['POINT_ID','stratum','LC1','LU1'],
  tileScale:16,
  scale: pix_export,
  geometries:false
});

var trainingexportdrive_EU_NW2b = s1stack.sampleRegions({
  collection: parcel_EU_NW2b,
  properties: ['POINT_ID','stratum','LC1','LU1'],
  tileScale:16,
  scale: pix_export,
  geometries:false
});

var trainingexportdrive_EU_NE1a = s1stack.sampleRegions({
  collection: parcel_EU_NE1a,
  properties: ['POINT_ID','stratum','LC1','LU1'],
  tileScale:16,
  scale: pix_export,
  geometries:false
});
var trainingexportdrive_EU_NE1b = s1stack.sampleRegions({
  collection: parcel_EU_NE1b,
  properties:['POINT_ID','stratum','LC1','LU1'],
  tileScale:16,
  scale: pix_export,
  geometries:false
});
var trainingexportdrive_EU_NE1c = s1stack.sampleRegions({
  collection: parcel_EU_NE1c,
  properties:['POINT_ID','stratum','LC1','LU1'],
  tileScale:16,
  scale: pix_export,
  geometries:false
});
var trainingexportdrive_EU_NE1d = s1stack.sampleRegions({
  collection: parcel_EU_NE1d,
  properties:['POINT_ID','stratum','LC1','LU1'],
  tileScale:16,
  scale: pix_export,
  geometries:false
});
var trainingexportdrive_EU_NE2 = s1stack.sampleRegions({
  collection: parcel_EU_NE2,
  properties: ['POINT_ID','stratum','LC1','LU1'],
  tileScale:16,
  scale: pix_export,
  geometries:false
});
var trainingexportdrive_EU_SW1 = s1stack.sampleRegions({
  collection: parcel_EU_SW1,
  properties: ['POINT_ID','stratum','LC1','LU1'],
  tileScale:16,
  scale: pix_export,
  geometries:false
});


var trainingexportdrive_EU_SW2 = s1stack.sampleRegions({
  collection: parcel_EU_SW2,
  properties:['POINT_ID','stratum','LC1','LU1'],
  tileScale:16,
  scale: pix_export,
  geometries:false
});

var trainingexportdrive_EU_SE1 = s1stack.sampleRegions({
  collection: parcel_EU_SE1,
  properties: ['POINT_ID','stratum','LC1','LU1'],
  tileScale:16,
  scale: pix_export,
  geometries:false
});
var trainingexportdrive_EU_SE2 = s1stack.sampleRegions({
  collection: parcel_EU_SE2,
  properties: ['POINT_ID','stratum','LC1','LU1'],
  tileScale:16,
  scale: pix_export,
  geometries:false
});

print('trainingexportdrive',trainingexportdrive_EU_SW1.limit(10))
print('trainingexportdrive size',trainingexportdrive_EU_SW1.size())

Export.table.toDrive({'collection': trainingexportdrive_EU_SE1, 
'description': 'S1_point_'+classtype+'_'+step+'days_'+pix_export+'m_1Jan-31Dec_EU_SE1', 
'fileNamePrefix':'S1_point_'+classtype+'_'+step+'days_'+pix_export+'m_1Jan-31Dec_EU_SE1_ratio-db',
folder:'GEE'
})

Export.table.toDrive({'collection': trainingexportdrive_EU_SE2, 
'description': 'S1_point_'+classtype+'_'+step+'days_'+pix_export+'m_1Jan-31Dec_EU_SE2', 
'fileNamePrefix':'S1_point_'+classtype+'_'+step+'days_'+pix_export+'m_1Jan-31Dec_EU_SE2_ratio-db',
folder:'GEE'
})

Export.table.toDrive({'collection': trainingexportdrive_EU_SW2, 
'description': 'S1_point_'+classtype+'_'+step+'days_'+pix_export+'m_1Jan-31Dec_EU_SW2', 
'fileNamePrefix':'S1_point_'+classtype+'_'+step+'days_'+pix_export+'m_1Jan-31Dec_EU_SW2_ratio-db', 
folder:'GEE'
})

Export.table.toDrive({'collection': trainingexportdrive_EU_SW1, 
'description': 'S1_point_'+classtype+'_'+step+'days_'+pix_export+'m_1Jan-31Dec_EU_SW1', 
'fileNamePrefix':'S1_point_'+classtype+'_'+step+'days_'+pix_export+'m_1Jan-31Dec_EU_SW1_ratio-db', 
folder:'GEE'
})

Export.table.toDrive({'collection': trainingexportdrive_EU_NE1a, 
'description': 'S1_point_'+classtype+'_'+step+'days_'+pix_export+'m_1Jan-31Dec_EU_NE1a', 
'fileNamePrefix':'S1_point_'+classtype+'_'+step+'days_'+pix_export+'m_1Jan-31Dec_EU_NE1a_ratio-db', 
folder:'GEE'
})

Export.table.toDrive({'collection': trainingexportdrive_EU_NE1b, 
'description': 'S1_point_'+classtype+'_'+step+'days_'+pix_export+'m_1Jan-31Dec_EU_NE1b', 
'fileNamePrefix':'S1_point_'+classtype+'_'+step+'days_'+pix_export+'m_1Jan-31Dec_EU_NE1b_ratio-db', 
folder:'GEE'
})

Export.table.toDrive({'collection': trainingexportdrive_EU_NE1c, 
'description': 'S1_point_'+classtype+'_'+step+'days_'+pix_export+'m_1Jan-31Dec_EU_NE1c', 
'fileNamePrefix':'S1_point_'+classtype+'_'+step+'days_'+pix_export+'m_1Jan-31Dec_EU_NE1c_ratio-db', 
folder:'GEE'
})

Export.table.toDrive({'collection': trainingexportdrive_EU_NE1d, 
'description': 'S1_point_'+classtype+'_'+step+'days_'+pix_export+'m_1Jan-31Dec_EU_NE1d', 
'fileNamePrefix':'S1_point_'+classtype+'_'+step+'days_'+pix_export+'m_1Jan-31Dec_EU_NE1d_ratio-db', 
folder:'GEE'
})

Export.table.toDrive({'collection': trainingexportdrive_EU_NE2, 
'description': 'S1_point_'+classtype+'_'+step+'days_'+pix_export+'m_1Jan-31Dec_EU_NE2', 
'fileNamePrefix':'S1_point_'+classtype+'_'+step+'days_'+pix_export+'m_1Jan-31Dec_EU_NE2_ratio-db', 
folder:'GEE'
})

Export.table.toDrive({'collection': trainingexportdrive_EU_NW2a, 
'description': 'S1_point_'+classtype+'_'+step+'days_'+pix_export+'m_1Jan-31Dec_EU_NW2a', 
'fileNamePrefix':'S1_point_'+classtype+'_'+step+'days_'+pix_export+'m_1Jan-31Dec_EU_NW2a_ratio-db',
folder:'GEE'
})

Export.table.toDrive({'collection': trainingexportdrive_EU_NW2b,
'description': 'S1_point_'+classtype+'_'+step+'days_'+pix_export+'m_1Jan-31Dec_EU_NW2b', 
'fileNamePrefix':'S1_point_'+classtype+'_'+step+'days_'+pix_export+'m_1Jan-31Dec_EU_NW2b_ratio-db',
folder:'GEE'
})

Export.table.toDrive({'collection': trainingexportdrive_EU_NW1, 
'description': 'S1_point_'+classtype+'_'+step+'days_'+pix_export+'m_1Jan-31Dec_EU_NW1', 
'fileNamePrefix':'S1_point_'+classtype+'_'+step+'days_'+pix_export+'m_1Jan-31Dec_EU_NW1_ratio-db',
folder:'GEE'
})

