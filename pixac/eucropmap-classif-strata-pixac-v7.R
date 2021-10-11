.libPaths(c('~/MyRlibs',.libPaths('/eos/jeodpp/home/users/verheas/')))

# Set the plot dir ----
if (Sys.info()["nodename"] =="d01ri1701915") {
  dir.data<-'/data/Dropbox/JRC/LANDSENSE/CASE-STUDY-8-LUCAS-COPERNICUS-CLASSIF/data/'
  dir.cropclassif<-'/data/Dropbox/JRC/LANDSENSE/CASE-STUDY-8-LUCAS-COPERNICUS-CLASSIF/data/S1_GS/v7/cropclassif/'
  dir.results<-'/data/'
}

if (Sys.info()["nodename"] =="jeodpp-terminal-power-04") { 
  dir.data<-'/eos/jeodpp/data/projects/REFOCUS/data/'
  dir.cropclassif<-'/eos/jeodpp/data/projects/REFOCUS/cropclassif/v7/'
  dir.results<-'/eos/jeodpp/data/projects/REFOCUS/data/S1_GS/'}

library(raster)
library(rgdal)
library(rgeos)

rasterOptions(maxmemory=10e+09,chunksize=5e+08,progress="text",timer = TRUE,tmptime=5) 



###########################
# PREPARE FILE LIST STRATUM
###########################

# Load all S1 upper left points tiles
load(paste0(dir.data,"S1_GS/v7/S1_list_extent__all_v2.Rdata"))

#pts<-SpatialPointsDataFrame(coords = df[,c('xmin','ymax')],data=df,proj4string = CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'))
#shapefile(pts, filename=paste0(dir.data,"S1_GS/v7/S1_list_extent__all.shp"),overwrite=TRUE)

# check duplicate
df$file<-as.character(df$file)


###########################
# CREATE POLYGONS FOOTPRINT OF S1 TILES
###########################

Ps.list<-c()
for (i in 1:nrow(df)) 
{
  print(i)
  # create the rectangle
  coords = matrix(c(df$xmin[i], df$ymax[i],
                    df$xmax[i], df$ymax[i],
                    df$xmax[i], df$ymin[i],
                    df$xmin[i], df$ymin[i],
                    df$xmin[i], df$ymax[i]), 
                  ncol = 2, byrow = TRUE)
  
  Ps = Polygons(list(Polygon(coords)),df[i,c('file')])
  Ps.list<-c(Ps.list,Ps)
}

#convert to spatial polygons
SPs<-SpatialPolygons(Ps.list)
proj4string(SPs)=crs('EPSG:4326')

S1_tiles_polygon<-SpatialPolygonsDataFrame(SPs,data.frame(df),match.ID = FALSE)

shapefile(S1_tiles_polygon, filename=paste0(dir.data,"S1_GS/v7/S1_tiles_polygon.shp"),overwrite=TRUE)
save(S1_tiles_polygon, file=paste0(dir.data,'S1_GS/v7/S1_tiles_polygon.RData'))

###########################
# INTERSECT S1 POLYGONS FOOTPRINT WITH EU
##########################

load(paste0(dir.data,'S1_GS/v7/S1_tiles_polygon.RData'))
#plot(S1_tiles_polygon[aoi,])

# open stratum and land EU28 mask
eu_28<-shapefile(paste0(dir.data,"S1_GS/v7/land_mask/EU_28_MASK/NUTS_RG_01M_2016_3035_LEVEL0_EU28.shp"))
eu_28_wgs84<- spTransform(eu_28,CRS("+init=epsg:4326"))

S1_tiles_polygon_eu<-S1_tiles_polygon[eu_28_wgs84,]

###########################
# INTERSECT S1 POLYGONS FROM EU WITH STRATUM TO EXTRACT THE LIST OF TILES TO PROCESS PER STRATUM
##########################
aoi_all<-shapefile(paste0(dir.data,"S1_GS/v7/land_mask/FinalStratum_astrid/checked_EU_Strata_1and2_Astrid_noHole.shp"))

stratum_i<-1
aoi<-aoi_all[aoi_all@data$stratum==stratum_i,]
S1_tiles_polygon_eu_stratum<-S1_tiles_polygon_eu[aoi,]

#aoi_3035<- spTransform(aoi,crs('EPSG:3035'))

for (stratum_i in seq(1,2)){
  aoi<-aoi_all[aoi_all@data$stratum==stratum_i,]
  S1_tiles_polygon_eu_stratum<-S1_tiles_polygon_eu[aoi,]
  print(paste(stratum_i,nrow(S1_tiles_polygon_eu_stratum)))
  write.table(S1_tiles_polygon_eu_stratum@data$file,
              file=paste0(dir.cropclassif,'list_rasters_eu_stratum_all',stratum_i,'.lst'),
              row.names = FALSE,col.names = FALSE,quote=FALSE)
  shapefile(S1_tiles_polygon_eu_stratum, filename=paste0(dir.cropclassif,'S1_tiles_polygon_process_Stratum_',stratum_i,".shp"),overwrite=TRUE)
}




# interestected tiles selection
# aoi1<-aoi_all[aoi_all@data$stratum==1,]
# aoi2<-aoi_all[aoi_all@data$stratum==2,]
# 
# S1_tiles_polygon_eu_stratum_1<-S1_tiles_polygon_eu[aoi1,]
# S1_tiles_polygon_eu_stratum_2<-S1_tiles_polygon_eu[aoi2,]
# 
# S1_tiles_on_both<-S1_tiles_polygon_eu_stratum_1[S1_tiles_polygon_eu_stratum_2,]


###########################
# IDENTIFY MISSING TILES
###########################


tiles<-as.data.frame(Sys.glob(paste0(dir.results,'S1_classif_v7/','S1_classif_',version,'_stratum_','*','/*level2.tif')))

colnames(tiles)<-c('tiles_path')
tiles$tiles_processed_id<-substr(basename(as.character(tiles$tiles_path) ),1,37)

# File to process
# stratum_i<-1
for (stratum_i in seq(1,2)){
  #tiles_toprocess<-as.data.frame(read.table(paste0('/eos/jeodpp/data/projects/REFOCUS/cropclassif/list_rasters_eu_stratum_all',stratum_i,'.lst'),stringsAsFactors=F))
  tiles_toprocess<-as.data.frame(read.table(paste0(dir.cropclassif,'list_rasters_eu_stratum_all',stratum_i,'.lst'),stringsAsFactors=F))
  
  
  colnames(tiles_toprocess)<-c('tiles_path')
  tiles_toprocess$tiles_toprocess_id<-substr(basename(as.character(tiles_toprocess$tiles_path)) ,1,37)
  
  missing<-tiles_toprocess$tiles_path[!tiles_toprocess$tiles_toprocess_id  %in% tiles$tiles_processed_id ]
  print(length(missing))
  #paste0('/eos/jeodpp/data/projects/REFOCUS/cropclassif/list_rasters_eu_stratum_',stratum_i,'_missing.lst')
  
  write.table(missing,
              file=paste0(dir.cropclassif,'list_rasters_eu_stratum_',stratum_i,'_missing.lst'),
              row.names = FALSE,col.names = FALSE,quote=FALSE)
}




###########################
# BUILD VRT PER STRATUM
###########################

version<-'v7'
stratum<-'1'

outRstName<-paste0(dir.results,'pixac_v7/','pixac_',version,'_stratum_',stratum,'.vrt')
inRstName<-paste0(dir.results,'pixac_',version,'_stratum_',stratum,'/*.tif')

system(
  paste('gdalbuildvrt','-srcnodata 0',outRstName,inRstName))


stratum<-'2'


outRstName<-paste0(dir.results,'pixac_v7/','pixac_',version,'_stratum_',stratum,'.vrt')
inRstName<-paste0(dir.results,'pixac_',version,'_stratum_',stratum,'/*.tif')

system(
  paste('gdalbuildvrt','-srcnodata 0',outRstName,inRstName))

###########################
# BUILD OVERVIEW 5% - TO CHECK MISSING TILES
###########################

stratum<-'1'

inRstName<-paste0(dir.results,'S1_classif_v7/','S1_classif_',version,'_stratum_',stratum,'.vrt')
percent=5
outRstName<-paste0(dir.results,'S1_classif_v7/','S1_classif_',version,'_stratum_',stratum,'-overview1.tif')

system(
  paste('gdal_translate',
        '-of','GTiff','-b 1',"-co", "TILED=YES", "-co", "COMPRESS=LZW", "-co", "BIGTIFF=YES",'--config GDAL_CACHEMAX 4096 ',
        '-outsize',paste(percent,'%',sep=''),paste(percent,'%',sep=''),
        inRstName,outRstName))

stratum<-'2'

inRstName<-paste0(dir.results,'S1_classif_v7/','S1_classif_',version,'_stratum_',stratum,'.vrt')
percent=5
outRstName<-paste0(dir.results,'S1_classif_v7/','S1_classif_',version,'_stratum_',stratum,'-overview1.tif')

system(
  paste('gdal_translate',
        '-of','GTiff','-b 1',"-co", "TILED=YES", "-co", "COMPRESS=LZW", "-co", "BIGTIFF=YES",'--config GDAL_CACHEMAX 4096 ',
        '-outsize',paste(percent,'%',sep=''),paste(percent,'%',sep=''),
        inRstName,outRstName))

# system(
#   paste('gdaladdo',outRstName ,'2 4 8 16'))


###########################
# BUILD FULL MOSAIC PER STRATUM
###########################


version<-'v7'
stratum<-'1'

inRstName<-paste0(dir.results,'pixac_v7/','pixac_',version,'_stratum_',stratum,'.vrt')
outRstName<-paste0(dir.results,'pixac_v7/','pixac_',version,'_stratum_',stratum,'-HR.tif')

system(
  paste('gdal_translate',
        '-of','GTiff','-b 1',"-co", "TILED=YES", "-co", "COMPRESS=LZW", "-co", "BIGTIFF=YES",'--config GDAL_CACHEMAX 4096 ',
        inRstName,outRstName))


stratum<-'2'

inRstName<-paste0(dir.results,'S1_classif_v7/','S1_classif_',version,'_stratum_',stratum,'.vrt')
outRstName<-paste0(dir.results,'S1_classif_v7/','S1_classif_',version,'_stratum_',stratum,'-HR.tif')

system(
  paste('gdal_translate',
        '-of','GTiff','-b 1',"-co", "TILED=YES", "-co", "COMPRESS=LZW", "-co", "BIGTIFF=YES",'--config GDAL_CACHEMAX 4096 ',
        inRstName,outRstName))





###########################
# CLIP BASED ON SHAPEFILE WITH GDALWARP
###########################


# # crop to polygon (test with overview)
# 
# version<-'v7'
# stratum<-'2'
# 
# aoi_all<-shapefile(paste0(dir.data,"S1_GS/v7/land_mask/FinalStratum_astrid/checked_EU_Strata_1and2_Astrid_noHole.shp"))
# 
# 
# shp.path<-paste0(dir.data,"S1_GS/v7/land_mask/FinalStratum_astrid/checked_EU_Strata_1and2_Astrid_noHole.shp")
# inRstName<-paste0(dir.results,'S1_classif_v7/','S1_classif_',version,'_stratum_',stratum,'-overview1.tif')
# outRstName<-paste0(dir.results,'S1_classif_v7/','S1_classif_',version,'_stratum_',stratum,'-overview1_clipped_blend2.tif')
# 
# system(paste('gdalwarp',
#              '-cutline',shp.path,
#              '-cl',strsplit(basename(shp.path),'\\.')[[1]][1],
#              "-cwhere",'"',paste0("stratum=\'",aoi_all@data$stratum[as.integer(stratum)],"\'"),'"',
#              '-crop_to_cutline','-overwrite','-dstnodata "0"',
#              '-of','GTiff',"-co", "TILED=YES", "-co", "COMPRESS=LZW", "-co", "BIGTIFF=YES",'--config GDAL_CACHEMAX 4096 ',
#              "-cblend 2", #buffer to cutline
#              inRstName,
#              outRstName))



# crop to polygon (full scale)

version<-'v7'
stratum<-'1'

aoi_all<-shapefile(paste0(dir.data,"S1_GS/v7/land_mask/FinalStratum_astrid/checked_EU_Strata_1and2_Astrid_noHole.shp"))


shp.path<-paste0(dir.data,"S1_GS/v7/land_mask/FinalStratum_astrid/checked_EU_Strata_1and2_Astrid_noHole.shp")

inRstName<-paste0(dir.results,'pixac_v7/','pixac_',version,'_stratum_',stratum,'_byte.tif')
outRstName<-paste0(dir.results,'pixac_v7/','pixac_',version,'_stratum_',stratum,'_byte_clipped.tif')

system(paste('gdalwarp',
             '-cutline',shp.path,
             '-cl',strsplit(basename(shp.path),'\\.')[[1]][1],
             "-cwhere",'"',paste0("stratum=\'",aoi_all@data$stratum[as.integer(stratum)],"\'"),'"',
             '-crop_to_cutline','-overwrite','-dstnodata "255 0"',
             '-of','GTiff',"-co", "TILED=YES", "-co", "COMPRESS=LZW", "-co", "BIGTIFF=YES",'--config GDAL_CACHEMAX 4096 ',
             "-cblend 2", #buffer to cutline
             '--config GDALWARP_IGNORE_BAD_CUTLINE YES',
             inRstName,
             outRstName))



version<-'v7'
stratum<-'2'

aoi_all<-shapefile(paste0(dir.data,"S1_GS/v7/land_mask/FinalStratum_astrid/checked_EU_Strata_1and2_Astrid_noHole.shp"))

shp.path<-paste0(dir.data,"S1_GS/v7/land_mask/FinalStratum_astrid/checked_EU_Strata_1and2_Astrid_noHole.shp")
inRstName<-paste0(dir.results,'pixac_v7/','pixac_',version,'_stratum_',stratum,'_byte.tif')
outRstName<-paste0(dir.results,'pixac_v7/','pixac_',version,'_stratum_',stratum,'_byte_clipped.tif')

system(paste('gdalwarp',
             '-cutline',shp.path,
             '-cl',strsplit(basename(shp.path),'\\.')[[1]][1],
             "-cwhere",'"',paste0("stratum=\'",aoi_all@data$stratum[as.integer(stratum)],"\'"),'"',
             '-crop_to_cutline','-overwrite','-dstnodata "255 0"',
             '-of','GTiff',"-co", "TILED=YES", "-co", "COMPRESS=LZW", "-co", "BIGTIFF=YES",'--config GDAL_CACHEMAX 4096 ',
             "-cblend 2", #buffer to cutline
             '--config GDALWARP_IGNORE_BAD_CUTLINE YES',
             inRstName,
             outRstName))


###########################
# BUILD FULL MOSAIC
###########################

version<-'v7'

# BUILD VRT CLIPPED
outRstName<-paste0(dir.results,'pixac_v7/','pixac_',version,'_HR_clipped.vrt')
inRstName<-paste0(dir.results,'pixac_v7/','pixac_',version,'_stratum_','*','_byte_clipped.tif')

system(
  paste('gdalbuildvrt','-srcnodata 255 0',outRstName,inRstName))


# create the mosaic
inRstName<-paste0(dir.results,'pixac_v7/','pixac_',version,'_HR_clipped.vrt')
outRstName<-paste0(dir.results,'pixac_v7/','pixac_',version,'_HR_clipped.tif')

system(
  paste('gdal_translate',
        '-of','GTiff','-b 1',"-co", "TILED=YES", "-co", "COMPRESS=LZW", "-co", "BIGTIFF=YES",'--config GDAL_CACHEMAX 4096 ',
        inRstName,outRstName))

#add pyramid
system(
  paste('gdaladdo -r nearest -ro',outRstName,'2 4 8 16 32 64 128 256'))


###########################
#  MOSAIC OVERVIEW  5%
###########################

inRstName<-paste0(dir.results,'pixac_v7/','pixac_',version,'_HR_clipped.tif')
percent=5
outRstName<-paste0(dir.results,'pixac_v7/','pixac_',version,'_HR_clipped_overview.tif')

system(
  paste('gdal_translate',
        '-of','GTiff','-b 1',"-co", "TILED=YES", "-co", "COMPRESS=LZW", "-co", "BIGTIFF=YES",'--config GDAL_CACHEMAX 4096 ',
        '-outsize',paste(percent,'%',sep=''),paste(percent,'%',sep=''),
        inRstName,outRstName))


###########################
# MASK FULL MOSAIC WITH CLC
###########################
version<-'v7'

inRst1<-paste0(dir.results,'S1_classif_v7/mask_CLC2018.tif')
inRst2<-paste0(dir.results,'S1_classif_v7/','S1_classif_',version,'_HR_clipped.tif')
outRst<-paste0(dir.results,'S1_classif_v7/mask_CLC2018_',version,'warped.tif')

refRst <- raster(inRst2)
tr <- paste('-tr',xres(refRst), yres(refRst))
te <- paste('-te',extent(refRst)@xmin, extent(refRst)@ymin,extent(refRst)@xmax,extent(refRst)@ymax)

system(paste('gdalwarp','-t_srs EPSG:4326',tr,te,'-ot Byte',
             inRst1,'-wm 2048', "-co", "TILED=YES", "-co", "COMPRESS=LZW", "-co", "BIGTIFF=YES",'--config GDAL_CACHEMAX 4096 ',
             outRst))




inRst1<-paste0(dir.results,'S1_classif_v7/mask_CLC2018_',version,'warped.tif')
inRst2<-paste0(dir.results,'S1_classif_v7/','S1_classif_',version,'_HR_clipped.tif')
outRst<-paste0(dir.results,'S1_classif_v7/','S1_classif_',version,'_HR_clipped_CLC_masked.tif')

system(
  paste('gdal_calc.py',
        '-A',inRst1,
        '--A_band=1',
        '-B',inRst2,
        '--B_band=1',
        paste(' --outfile=', outRst,sep=''),
        paste(' --calc=','\"','(B*(A<1))','\"',sep=''),
        paste('--NoDataValue=',0,sep=''),
        paste('--type=','Int16',sep=''),
        paste('--co=','\"','COMPRESS=LZW','\"',sep=''),
        paste('--co=','\"','BIGTIFF=YES','\"',sep=''),
        paste('--co=','\"','TILED=YES','\"',sep='')
  )
)

###########################
# ADD PYRAMIDS
###########################

outRst<-paste0(dir.results,'S1_classif_v7/','S1_classif_',version,'_HR_clipped_CLC_masked.tif')

#add pyramid
system(
  paste('gdaladdo -r nearest -ro',outRst,'2 4 8 16 32 64 128 256'))


/eos/jeodpp/data/projects/BDA/REFOCUS/masked/eucropmap_masked.vrt


# create the mosaic
inRstName<-paste0('/eos/jeodpp/data/projects/BDA/REFOCUS/masked/eucropmap_masked.vrt')
outRstName<-paste0('/scratch/eucropmap_masked.tif')

system(
  paste('gdal_translate',
        '-of','GTiff','-b 1',"-co", "TILED=YES", "-co", "COMPRESS=LZW", "-co", "BIGTIFF=YES",'--config GDAL_CACHEMAX 4096 ',
        inRstName,outRstName))

#add pyramid
system(
  paste('gdaladdo -r nearest -ro',outRstName,'2 4 8 16 32 64 128 256'))


# # # # # # # # # # # # # # 

version<-'v7'

# BUILD VRT MASKED
outRstName<-paste0(dir.results,'pixac_v7/','pixac_',version,'_byte_masked.vrt')
inRstName<-paste0(dir.results,'pixac_v7/','masked/','eucropmap_pixac_','masked_','*','_1600.tif')

system(
  paste('gdalbuildvrt','-srcnodata 0',outRstName,inRstName))

# create the mosaic
inRstName<-paste0(dir.results,'pixac_v7/','pixac_',version,'_byte_masked.vrt')
outRstName<-paste0(dir.results,'pixac_v7/','pixac_',version,'_byte_masked.tif')

system(
  paste('gdal_translate',
        '-of','GTiff','-b 1',"-co", "TILED=YES", "-co", "COMPRESS=LZW", "-co", "BIGTIFF=YES",'--config GDAL_CACHEMAX 4096 ',
        inRstName,outRstName))

#add pyramid
system(
  paste('gdaladdo -r nearest -ro',outRstName,'2 4 8 16 32 64 128 256'))