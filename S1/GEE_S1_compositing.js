//=======================================================================================================================================
// EC-JRC 2021 - Code to create a S1 10-days VV and VH composite to Google Cloud Storage
// Related to the paper: "From parcel to continental scale -- 
// A first European crop type map based on Sentinel-1 and LUCAS Copernicus in-situ observations" 
// by Raphaël d'Andrimont, Astrid Verhegghen, Guido Lemoine, Pieter Kempeneers, Michele Meroni, Marijn van der Velde. 
//=======================================================================================================================================

// A / FUNCTIONS
//=============================================
// Functions to convert from/to dB
function toNatural(img) {
    return ee.Image(10.0).pow(img.select('..').divide(10.0)).copyProperties(img, ['system:time_start'])
  }
  
  function toDB(img) {
    return ee.Image(img).log10().multiply(10.0).copyProperties(img, ['system:time_start']);
  }
  
  // Remove ugly edges
  function maskEdge(img) {
    var mask = img.select(0).unitScale(-25, 5).multiply(255).toByte().connectedComponents(ee.Kernel.rectangle(1,1), 100);
    return img.updateMask(mask.select(0).abs());  
  }
  
  // B/ S1 COMPOSITE
  //=============================================
  
  // 1. Select the area of interest
  //=============================================
  //list of countries
  
  var europe = ee.List(['Austria','France','Belgium','Denmark','Greece','Germany','United Kingdom','Italy','Lithuania','Malta','Netherlands','Poland','Romania','Portugal','Slovenia','Bulgaria','Croatia',
          'Hungary','Cyprus','Czechia','Sweden','Ireland','Finland','Spain','Latvia','Estonia','Luxembourg','Slovakia'])
  
  //var europe = ee.List(['Poland'])
  
  var country=countries.filter(ee.Filter.inList('country_na',europe))
  print('country',country)
  
  Map.addLayer(country,{},'country',true)
  var aoi = country.union();
  print(aoi)
  
  //2. Select the dates and time step
  //=============================================
  // select dates
  var start_date = '2018-01-01'
  var end_date = '2018-12-31'
  
  //Set  the time step
  var step = 10 // in days (time window for averaging)
  
  //Spatial resolution
  var res=10
  
  //3. Processing
  //============================================
  // 3.1 / Data pre-processing
  
  // get the data from S1 (VV pol.)
  var s1 = ee.ImageCollection('COPERNICUS/S1_GRD').filterMetadata('instrumentMode', 'equals', 'IW')
    .filter(ee.Filter.eq('transmitterReceiverPolarisation', ['VV', 'VH']))//this seems to be always the case for IW mode
    .filterBounds(aoi)
    .filterDate(start_date, end_date)
    .sort('system:time');
  
  
  // Remove ugly edges
  s1 = s1.map(maskEdge)
  
  // Average are made from natural (non-logarithmic) values
  s1 = s1.map(toNatural)
  
  // Olha Danylo's procedure to create weekly means (adapted)
  var days = ee.List.sequence(0, ee.Date(end_date).difference(ee.Date(start_date), 'day'), step).
    map(function(d) { return ee.Date(start_date).advance(d, "day") })
  
  print ('days', days)
  
  var dates = days.slice(0,-1).zip(days.slice(1))
  print ('dates', dates)
  
  // 3.2 / Temporal compositing
  var s1res = dates.map(function(range) {
    var dstamp = ee.Date(ee.List(range).get(0)).format('YYYYMMdd')
    var temp_collection = s1.filterDate(ee.List(range).get(0),
    ee.List(range).get(1)).mean().select(['VV', 'VH'], [ee.String('VV_').cat(dstamp), ee.String('VH_').cat(dstamp)])
    return temp_collection
  })
  
  
  //transform back to DB
  s1res=s1res.map(toDB)
  
  // Convert ImageCollection to image stack
  function stack(i1, i2)
  {
    return ee.Image(i1).addBands(ee.Image(i2))
  }
  
  var s1stack = s1res.slice(1).iterate(stack, s1res.get(0))
  
  //toImage
  s1stack = ee.Image(s1stack)
  
  //transform the image to float to reduce size
  s1stack = s1stack.toFloat()
  
  // 3.3 / Clip to the study area
  //clip to the countries of EU-28
  s1stack = s1stack.clip(country)
  
  var viz ={bands:['VV_20180511','VV_20180610','VV_20180710'],min:-25,max:0}
  
  Map.addLayer(s1stack,viz,'s1stack_europe',0)
  
  