
/*****  *****/

//var geometry = countries
var nigeria =  country.filter(ee.Filter.eq('ADM0_NAME', 'Nigeria'));

Map.centerObject(nigeria, 6);
Map.addLayer(nigeria, {min: 0, max: 1, palette: ['red']}, 'nigeria');


//*****   (1)   PERMANET WATER BODIES  *****//


// Yearly water occurrence (2019–2024) from JRC GSW
var gsw = ee.ImageCollection("JRC/GSW1_4/YearlyHistory")
  .filterBounds(nigeria)
  .filter(ee.Filter.calendarRange(2019, 2024, 'year'))
  // YearlyHistory: 2 = water. Make it binary and name it 'occurrence' for consistency.
  .map(function(img){ return img.eq(2).rename('occurrence'); })
  .mean()                   // mean of yearly binary water -> fraction of years with water
  .multiply(100)            // convert to 0–100 %
  .clip(nigeria);

// Map.addLayer(gsw, {min: 0, max: 100, palette: ['black', 'blue']}, 'Water');

// Permanent water: present in >50% of the years
var permanent = gsw.gt(50);
Map.addLayer(permanent.selfMask(), waterimageVisParam2, 'Permanent water');



//*****   (2)   DISTANCE TO ROAD  *****//

var roads = Roads
var distanceToRoads = roads
    .distance({searchRadius: 100000, maxError: 1})
    .clip(nigeria);
    
Map.addLayer(distanceToRoads, roadimageVisParam, 'Distance_to_Roads');



//*****  (3)  DISTANCE TO RIVER  *****//

var rivers = Rivers
var distanceToRivers = rivers
    .distance({searchRadius: 100000, maxError: 1})
    .clip(nigeria);
    

Map.addLayer(distanceToRivers, roadimageVisParam, 'Distance_to_Rivers');



//*****  (4)  DRAINAGE DENSITY  *****//

var gridScale = 1000; // size of each analysis cell in meters

// Convert rivers to raster with line length burned into pixels
var riverRaster = rivers
  .map(function(f) {
    var lengthKm = f.length().divide(1000);  // convert m to km
    return f.set('length_km', lengthKm);
  })
  .reduceToImage({
    properties: ['length_km'],
    reducer: ee.Reducer.sum()
  })
  .unmask(0);

//Create fishnet/grid image (1 km resolution)
var grid = ee.Image.pixelArea()
  .divide(1e6) // m² to km²
  .reproject({crs: 'EPSG:4326', scale: gridScale})
  .clip(nigeria);

// Compute Drainage Density (length of river per km²)
var drainagedensity = riverRaster.divide(grid)
  .rename('drainage_density_km_per_km2');

Map.addLayer(drainagedensity, {n: 0, max: 5, palette: ['white', 'blue', 'darkblue']},'Drainage Density');



//*****  (4)  TEMPERATURE  *****//


var tempCollection = terraclimate
  .filterBounds(nigeria)
  .filterDate('2019-01-01', '2024-12-31')
  .select('tmmx'); // 'tmmx' = monthly maximum temperature (in tenths of °C)

// Convert from tenths of °C to °C
var tempCelsius = tempCollection.map(function(img) {
  return img.divide(10).copyProperties(img, ['system:time_start']);
});

// Calculate mean max temperature over the period
var temp = tempCelsius.mean().clip(nigeria);

Map.addLayer(temp, {min: 20, max: 45, palette: ['blue', 'lightblue', 'yellow', 'orange', 'red']
}, 'Mean Max Temp (°C)');



//*****  (4)  RAINFALL  *****//

var rainfall = chrips.filterDate('2020-01-01', '2024-12-31')
                          .select('precipitation')
                          .mean()
                          .clip(nigeria);
                          
Map.addLayer(rainfall, rainimageVisParam, 'Rainfall');



//*****  (4)  SOIL  *****//

var soil = Soilcollection.select('b0').clip(nigeria);

Map.addLayer(soil,soilimageVisParam, 'Soil pH');
             


//*****  (5)  LAND USE LAND COVER  *****//


var lulc = esa.mean()
              .clip(nigeria);
                   
                   
// Remap values: 10 - Tree cover,20 - Shrubland, 30 - Grassland, 40 - Cropland
//50 - Built-up, 60 - Bare/sparse vegetation, 80 - Permanent water bodies, 95 - Mangroves

var grouped = lulc.remap(
  [10, 20, 30, 40, 50, 60, 80, 95],  // from
  [1, 2, 3, 4, 5, 6, 7, 8]      // to
);

Map.addLayer(grouped, {min: 1, max: 7, palette: ['#006400','#FFBB22', '#FFFF4C', 
'#F096FF', '#FA0000','#B4B4B4', '#0064C8', '#00CF75']}, 'LULC');




//*****  (6)  ELEVATION  *****//

var srtm = ee.Image('USGS/SRTMGL1_003')
.clip(nigeria);

Map.addLayer(srtm, ele_imageVisParam , 'Elevation');



//*****  (7)  TOPOGRAPHIC POSITION INDEX (TPI)  *****//


var TPI = srtm.subtract(srtm.focalMean(5)).rename('TPI');

Map.addLayer(TPI, {min: -50, max: 50,  palette: ['blue', 'white', 'red']}, 'TPI');



//*****  (8)  SLOPE  *****//

// Compute Slope (degrees)

var slope = ee.Terrain.slope(srtm);

Map.addLayer(slope, slope_imageVisParam, 'Slope');


//*****  (9)  ASPECT  *****//

// Compute Aspect (degrees)

var aspect = ee.Terrain.aspect(srtm);

Map.addLayer(aspect, aspect_imageVisParam, 'Aspect');


//*****  (10)  CURVATURE  *****//

//  5. Compute Curvature (using built-in function)
// Define Sobel kernels for 2nd derivatives
var kernelX = ee.Kernel.fixed(3, 3, [
  [1, -2, 1],
  [2, -4, 2],
  [1, -2, 1]
], -1, -1);

var kernelY = ee.Kernel.fixed(3, 3, [
  [1, 2, 1],
  [-2, -4, -2],
  [1, 2, 1]
], -1, -1);

// Compute 2nd derivatives (approximate curvature)
var d2z_dx2 = srtm.convolve(kernelX);
var d2z_dy2 = srtm.convolve(kernelY);

// Curvature = sum of 2nd derivatives
var curvature = d2z_dx2.add(d2z_dy2).rename('curvature');

// Visualize curvature
Map.addLayer(curvature, curv_imageVisParam, 'Curvature of Nigeria');



//*****  (11)  TOPOGRAPHIC WETNESS INDEX (TWI)  *****//

// Using HydroSHEDS 15 arc-second (~500m) flow accumulation

var slopeRad = slope.multiply(Math.PI).divide(180); // Convert to radians
var tanSlope = slopeRad.tan();

// Avoid divide-by-zero by setting very small value where slope = 0
var safeTanSlope = tanSlope.where(tanSlope.eq(0), 0.001);

var pixelArea = ee.Image.pixelArea();

// -------------------- Load Flow Accumulation ---------------------
var flowAcc = ee.Image("WWF/HydroSHEDS/15ACC")
                .resample('bilinear')
                .reproject({crs: srtm.projection(), scale: 500})
                .clip(nigeria)
                .unmask(1);  // Prevent nulls

// Mask invalid or zero accumulation
flowAcc = flowAcc.updateMask(flowAcc.gt(0));

// -------------------- TWI Calculation ---------------------
var specificCatchmentArea = flowAcc.multiply(pixelArea);
var TWI = specificCatchmentArea.divide(safeTanSlope)
          .add(1).log().rename('TWI')
          .updateMask(safeTanSlope.gt(0));

Map.addLayer(TWI, twi_imageVisParam, 'Topographic Wetness Index');



//*****  (12)  STREAM POWER INDEX (SPI)  *****//

// Specific catchment area (As)
var As = flowAcc.multiply(pixelArea);

// SPI calculation
var SPI = As.multiply(slopeRad.tan()).add(1).log().rename('SPI');

// Visualize SPI
Map.addLayer(SPI, {min: 4, max: 12, palette: ['yellow', 'orange', 'red']}, 'Stream Power Index');



//*****  (13) SCORING EACH MAPS *****//

//DISTANCE TO ROAD SCORING

var riverScore = ee.Image(0)
  .where(distanceToRivers.lte(d5), 5)
  .where(distanceToRivers.gt(d5).and(distanceToRivers.lte(d10)), 4)
  .where(distanceToRivers.gt(d10).and(distanceToRivers.lte(d15)), 3)
  .where(distanceToRivers.gt(d15).and(distanceToRivers.lte(d18)), 2)
  .where(distanceToRivers.gt(d18), 1);
  
Map.addLayer(riverScore, {min: 1, max:5}, 'river_Score');

var roadScore = ee.Image(0)
  .where(distanceToRoad.lte(100), 5)
  .where(distanceToRoad.gt(100).and(distanceToRoad.lte(500)), 4)
  .where(distanceToRoad.gt(500).and(distanceToRoad.lte(2000)), 3)
  .where(distanceToRoad.gt(2000).and(distanceToRoad.lte(5000)), 2)
  .where(distanceToRoad.gt(5000), 1);

Map.addLayer(roadScore, {min: 1, max:5}, 'roadScore');


// // Elevation Scoring

var e5  = 326.6; var e10 = 634.2; var e15 = 941.8; var e18 = 1249.4;

// Low elevation = more vulnerable (5), high = less (1)

var elevScore = ee.Image(0)
  .where(elevation.lte(100), 5)
  .where(elevation.gt(100).and(elevation.lte(300)), 4)
  .where(elevation.gt(300).and(elevation.lte(600)), 3)
  .where(elevation.gt(600).and(elevation.lte(1000)), 2)
  .where(elevation.gt(1000), 1);


Map.addLayer(elevScore, {}, 'elev Score');
  
 
// Rainfall Score
var rainScore = ee.Image(0)
  .where(rainfall.gt(6), 5)
  .where(rainfall.gt(5).and(rainfall.lte(6)), 4)
  .where(rainfall.gt(4).and(rainfall.lte(5)), 3)
  .where(rainfall.gt(3).and(rainfall.lte(4)), 2)
  .where(rainfall.lte(3), 1);

Map.addLayer(rainScore, {}, 'rainScore');

// Temperature Score

var tempScore = ee.Image(0)
  .where(temp.lte(28), 5)
  .where(temp.gt(28).and(temp.lte(30)), 4)
  .where(temp.gt(30).and(temp.lte(32)), 3)
  .where(temp.gt(32).and(temp.lte(34)), 2)
  .where(temp.gt(34), 1);

Map.addLayer(tempScore, {}, 'tempScore');

//*****  (11) EXPORTING THE MAPS  *****//


function exportToDrive(img, name) {
  Export.image.toDrive({
    image: img.clip(nigeria),             // keep it tight to your AOI
    description: name,                    // task name & file base name
    folder: 'GEE',                        // Drive folder
    scale: 500,                           // adjust as needed
    region: nigeria.geometry(),           // geometry, not FC
    fileFormat: 'GeoTIFF',
    crs: 'EPSG:4326',
    maxPixels: 1e13
  });
}

// ---- Add your 8 layers here ----
// Tip: use Int for 1–5 scores (smaller files), Float for continuous layers
var layersToExport = [
  { image: distanceToRivers.toFloat(), name: 'distanceToRivers' },
  { image: distanceToRoads.toFloat(),  name: 'distanceToRoads'  },
  { image: rainfall.toFloat(),     name: 'rainfall'    },
  { image: soil.toFloat(),         name: 'soilph'        },
  { image: srtm.toFloat(),    name: 'elevation'   },
  { image: slope.toFloat(),        name: 'slope'       },
  { image: TPI.toFloat(),          name: 'tpi'         },
  { image: curvature.toFloat(),         name: 'curvature' },
  { image: TWI.toFloat(),         name: 'twi'},
  { image: SPI.toFloat(),         name: 'spi' },
  { image: aspect.toFloat(),         name: 'aspect' },
  { image: temp.toFloat(),         name: 'temperature' },
  { image: lulc.toFloat(),         name: 'lulc' },
  { image: drainagedensity.toFloat(),         name: 'drainagedensity' },
  { image: permanent.toFloat(),         name: 'flood' }
];

// ---- Queue all exports ----
layersToExport.forEach(function (item) {
  exportToDrive(item.image, item.name);
});