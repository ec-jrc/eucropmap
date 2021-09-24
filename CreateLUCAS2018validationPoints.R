
# Set the plot dir ----
if(.Platform$OS.type == "unix") {dir.dropbox<-'/data/Dropbox'}
if (.Platform$OS.type == "windows") { dir.dropbox<-'C:/Users/rdand/Dropbox'}
plot.dir=paste0(dir.dropbox,'/JRC/LANDSENSE/CASE-STUDY-8-LUCAS-COPERNICUS-CLASSIF/figures/')
work.dir=paste0(dir.dropbox,'/JRC/LANDSENSE/CASE-STUDY-8-LUCAS-COPERNICUS-CLASSIF/')

work.dir<-'/data/Dropbox/JRC/LANDSENSE/CASE-STUDY-8-LUCAS-COPERNICUS-CLASSIF/data-published/'
plot.dir<-paste0(work.dir,'figure/')

library(raster)

###############################################
# 0/ Downloadnput files  ====
############################################### 

url_lucasharmo2018_survey<-'https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/LUCAS/LUCAS_harmonised/1_table/lucas_harmo_uf_2018.zip'
url_lucas_geometry<-'https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/LUCAS/LUCAS_harmonised/2_geometry/LUCAS_th_geom.zip'
url_lucas_copernicus<-'https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/LUCAS/LUCAS_2018_Copernicus.zip'

# Download LUCAS points geometries (i.e.the grid)
download.file(url_lucas_geometry,destfile=file.path(work.dir,basename(url_lucas_geometry)))

# Download LUCAS harmo 2018 table survey
download.file(url_lucasharmo2018_survey,destfile=file.path(work.dir,basename(url_lucasharmo2018_survey)))

# Download LUCAS Copernicus 2018 survey
download.file(url_lucas_copernicus,destfile=file.path(work.dir,basename(url_lucas_copernicus)))

###############################################
# 1/ Read the points  gemoetry ( it could take take about 5 min as the db is  big)  ====
############################################### 

# Read the shapefile file with the theroretical points locations
unzip(file.path(work.dir,basename(url_lucas_geometry)),exdir=file.path(work.dir))
lucas_geom<-shapefile(file.path(work.dir,'LUCAS_th_geom.shp')) # read the grid
lucas_geom<-lucas_geom[lucas_geom@data$YEAR==2018,] 

###############################################
# 2/ Read the LUCAS 2018 survey csv and merge it with the geometries  ====
############################################### 
# load the 2018 survey results 
unzip(file.path(work.dir,basename(url_lucasharmo2018_survey)),exdir=file.path(work.dir))
lucas_table<-read.csv(file.path(work.dir,'lucas_harmo_uf_2018.csv')) 

spdf<-merge(x = lucas_geom, y = lucas_table, by.x = "POINT_ID",by.y = "point_id") # merge them
colnames(spdf@data)<-toupper(colnames(spdf@data))

###############################################
# 3/ Read the LUCAS 2018 Copernicus  ====
############################################### 
# Load the Copernicus LUCAS
unzip(file.path(work.dir,basename(url_lucas_copernicus)),exdir=file.path(work.dir))
lucas_table<-read.csv(file.path(work.dir,'LUCAS_2018_Copernicus_attributes.csv')) 

lucas_copernicus<-read.csv(file.path(work.dir,'LUCAS_2018_Copernicus_attributes.csv'))
lucas_copernicus_point_id<-lucas_copernicus$POINT_ID[lucas_copernicus$COPERNICUS_CLEANED]

###############################################
# 3/ Remove the LUCAS points with Copernicus survey ====
############################################### 
# Drop the points used as Copernicus points (for training)
spdf<-spdf[!(spdf@data$POINT_ID %in% lucas_copernicus_point_id),]

###############################################
# 4/ Filter LUCAS point to have high quality points ====
############################################### 
condition1<-spdf@data$OBS_TYPE=='In situ < 100 m' # field survey, point visible
condition2<-spdf@data$OBS_DIRECT=='On the point' # observation neither to North or to West but on point
condition3<-((spdf@data$PARCEL_AREA_HA != '< 0.5 ha'  ) & (spdf@data$PARCEL_AREA_HA != 'Not relevant' ))# parcel area estimated > 0.5 hectare
condition4<-spdf@data$LC1_PERC== '> 75 %' # percentage of coverage of the LC class greater than 75%
condition12434<-condition1 & condition2 & condition3  & condition4

# Removing from validation the training points
validation<-spdf[condition12434,]

###############################################
# 5/ Add the stratum ====
############################################### 
biome<-shapefile(paste0(work.dir,"stratum/eucropmap_strata.shp"))

for (i_stratum in c(1,2)){
  stratum<-biome[biome@data$stratum==i_stratum,]
  validation_stratum<-validation[stratum,]
  validation_stratum@data$stratum<-rep(i_stratum,nrow(validation_stratum))
  
  if (i_stratum==1){
    vali<-validation_stratum
  }
  if (i_stratum>1){
    vali<- rbind(vali,validation_stratum)
  }
  
}

validation<-vali

###############################################
# 6/ Save the validation points as Rdata, shp and csv====
############################################### 

# Save as Rdata
save(validation, file=paste0(work.dir,'LUCAS_2018_points_validation.RData'))

dir.create(paste0(work.dir,'LUCAS_2018_Copernicus_validation'))

# Export the polygons alone
field.to.export<-c("POINT_ID")
export<-validation
export@data<-data.frame(POINT_ID= export@data$POINT_ID)
shapefile(export, filename=paste0(work.dir,'/LUCAS_2018_Copernicus_validation/LUCAS_2018_validations.shp'),overwrite=TRUE)

# Export the attributes as a csv dataset
write.csv(validation@data,
          paste0(work.dir,'/LUCAS_2018_Copernicus_validation/LUCAS_2018_Copernicus_validation_attributes.csv'), row.names = F)

# zip the shapefile and csv in one file
zip(zipfile =paste0(work.dir,'LUCAS_2018_Copernicus_validation.zip'), files = paste0(work.dir,'/LUCAS_2018_Copernicus_validation/'),extras = '-j')

