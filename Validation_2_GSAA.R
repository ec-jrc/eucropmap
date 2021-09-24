# Set the plot dir ----
if(.Platform$OS.type == "unix") {dir.dropbox<-'/data/Dropbox'}
if (.Platform$OS.type == "windows") { dir.dropbox<-'C:/Users/rdand/Dropbox'}
plot.dir=paste0(dir.dropbox,'/JRC/LANDSENSE/CASE-STUDY-8-LUCAS-COPERNICUS-CLASSIF/figures/')
work.dir=paste0(dir.dropbox,'/JRC/LANDSENSE/CASE-STUDY-8-LUCAS-COPERNICUS-CLASSIF/')

# Set the plot dir ----
if (Sys.info()["nodename"] =="d01ri1701915") {
  dir.data<-'/data/Dropbox/JRC/LANDSENSE/CASE-STUDY-8-LUCAS-COPERNICUS-CLASSIF/data/'
  dir.cropclassif<-'/data/Dropbox/JRC/LANDSENSE/CASE-STUDY-8-LUCAS-COPERNICUS-CLASSIF/data/S1_GS/v7/cropclassif/'
  dir.results<-'/data/'
}
if (Sys.info()["nodename"] =="jeodpp-terminal-151p-01") { 
  dir.data<-'/eos/jeodpp/data/projects/REFOCUS/data/'
  dir.cropclassif<-'/eos/jeodpp/data/projects/REFOCUS/cropclassif/v7/'
  dir.results<-'/eos/jeodpp/data/projects/REFOCUS/data/S1_GS/'}




###########################
##### CROP VRT OVER REGIONS WITH GSAA TABLES#####
###########################

version<-'v7'
library(raster)
shp<-shapefile('/data/Dropbox/JRC/LANDSENSE/CASE-STUDY-8-LUCAS-COPERNICUS-CLASSIF/data/S1_GS/v7/LPIS/lucas_extract.shp')

regions<-shp@data$GSAAtable
regions<-c("bevl2018","dk2018" , "frcv2018","nld2018", "nrw2018" ,   "si2018")


region<-"frcv2018"

for (region in regions){
  print(region)
  poolyg<-shp[shp@data$GSAAtable==region,]
  
  xmin<-poolyg@bbox[1,][1]
  xmax<-poolyg@bbox[1,][2]
  ymin<-poolyg@bbox[2,][1]
  ymax<-poolyg@bbox[2,][2]
  
  proj_win<-paste(xmin, ymax, xmax, ymin)
  
  
  
  inRstName<-paste0('.',dir.results,'S1_classif_v7/lpis/','S1_classif_',version,'_HR_clipped.tif')
  outRstName<-paste0('.',dir.results,'S1_classif_v7/lpis/s1-classif_',version,'_',region,'.vrt')
  
  system(
    paste('docker run -i --rm -v /:/usr/src/app glemoine62/dias_py',
          'gdal_translate','-projwin',proj_win,'-of VRT',inRstName,outRstName))
  
  
}

###########################
##### CROP VRT OVER REGIONS WITH GSAA TABLES#####
###########################



setwd('/data/S1_classif_v7/lpis')

region<-c("bevl2018")
for (region in regions){
  print(region)
  
  poolyg<-shp[shp@data$GSAAtable==region,]
  crop_code<-poolyg$cropcode
  
  system(
    paste('docker run -i --rm -v`pwd`:/usr/src/app glemoine62/dias_py',
          'python postgisGeoTiffExtract.py',
          paste0('s1-classif_v7_',region,'.vrt'),
          region,
          crop_code,
          '0.01 >',paste0(region,'_7th.txt')))
}


# system(
#   paste('docker run -it --rm -v`pwd`:/usr/src/app glemoine62/dias_py',
#         'python postgisGeoTiffExtract.py s1-classif_v7_bevl2018.vrt bevl2018 hfdtlt 0.01 > bevl2018_7th.txt'))


###########################
##### CROP VRT OVER REGIONS WITH GSAA TABLES#####
###########################

setwd('/data/Dropbox/JRC/LANDSENSE/CASE-STUDY-8-LUCAS-COPERNICUS-CLASSIF/data/LPIS_extract')

regions<-c("bevl2018","dk2018" , "frcv2018","nld2018", "nrw2018" ,   "si2018")



for (region in regions){
  system(
    paste('python3',
          '/data/Dropbox/JRC/LANDSENSE/CASE-STUDY-8-LUCAS-COPERNICUS-CLASSIF/scripts/Guido_DIAS/LPIS_lucas_conf_matrix/confusionMatrix.py',
          paste0('/data/Dropbox/JRC/LANDSENSE/CASE-STUDY-8-LUCAS-COPERNICUS-CLASSIF/data/LPIS_extract/',region,'_7th.txt'))
  )
}




###########################
##### EXTRACT PA, UA FOR CONF MATRIX
###########################

work.dir<-'/data/Dropbox/JRC/LANDSENSE/CASE-STUDY-8-LUCAS-COPERNICUS-CLASSIF/data-published/'

lucas_cropcodes<-read.csv('LPIS_extract/lucas_cropcodes_To_PG.csv')
regions<-unique(lucas_cropcodes$parceltable)



regions<-c("dk2018" , "nld2018", "nrw2018" , "bevl2018",  "si2018","frcv2018")

df_all<-data.frame()
for (region in regions){
  cm_count<-read.csv(paste0(work.dir,'LPIS_extract/',region,'_7th_abscount.csv'))
  #remove 300 and 500
  colnames(cm_count)[2:ncol(cm_count)]<-data.frame(strsplit(colnames(cm_count)[2:ncol(cm_count)],'X'))[2,]
  
  cm_count<-cm_count[cm_count$X<300,]
  label_crop<-cm_count$X
  cm_count<-cm_count[,-1]
  
  cm_count<-cm_count[,colnames(cm_count)<300]
  
  
  cm.with.sum<- cbind(cm_count,TOTAL=rowSums(cm_count))
  cm.with.sums.1<- rbind(cm.with.sum,TOTAL=colSums(cm.with.sum))
  
  
  # calculate PA and UA
  UA<-c()
  PA<-c()
  for (i in seq(1,ncol(cm_count))){
    j<-i
    UA[i]<-cm_count[i,j]/cm.with.sums.1[i,ncol(cm.with.sums.1)]
    PA[j]<-cm_count[i,j]/cm.with.sums.1[nrow(cm.with.sums.1),j]
  }
  
  df<-data.frame(label_crop,UA, PA,region=rep(region,length(PA)))
  df_all <-bind_rows(df_all,
                     df)
  # CM wit UA and PA
  cm.with.sums.1.UA<-cbind(cm.with.sums.1,UA=c(round(UA,digit=4),''))
  cm.with.sums.1.PA<-rbind(cm.with.sums.1.UA,PA=c(round(PA,digit=4),'',''))
  rownames(cm.with.sums.1.PA)<-c(colnames(cm.with.sums.1.PA)[1:(ncol(cm.with.sums.1.PA)-1)],'PA')
  # write.csv(cm.with.sums.1.PA,paste0(work.dir,'table/LPIS_CM/',region,'_7th_CM_PA_UA.csv'))
  # print(xtable(cm.with.sums.1.PA,include.rownames=TRUE,caption=paste('Confusion matrix for', region)))
  
}

save(df_all,file=paste0(work.dir,'LPIS_Extract','GSAA_UA_PA.Rdata'))





