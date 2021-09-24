
library(raster)
library(stringr)
library(ggplot2)
library(reshape2)


# Set the plot dir ----

work.dir<-'/data/Dropbox/JRC/LANDSENSE/CASE-STUDY-8-LUCAS-COPERNICUS-CLASSIF/data-published/'

###########################
# EXTRACTION EN EUCROPMAP DONE IN QGIS LEVEL 2
###########################

# Fix the geometries before running Zonal histogram
# classif<-raster(file.path(work.dir,'/result/v7/S1_classif_v7_HR_clipped_CLC_masked.tif'))
# url_nuts1<-'https://ec.europa.eu/eurostat/cache/GISCO/distribution/v2/nuts/shp/NUTS_RG_01M_2016_4326_LEVL_2.shp.zip'

# # load the LUCAS copernicus polygon
# url_nuts1<-'https://ec.europa.eu/eurostat/cache/GISCO/distribution/v2/nuts/shp/NUTS_RG_01M_2016_4326_LEVL_2.shp.zip'
# 
# # Download nuts1 delineation
# if(!file.exists(file.path(work.dir,'NUTS_RG_01M_2016_4326_LEVL_2.shp'))){
#   
#   download.file(url_nuts1,destfile=file.path(work.dir,basename(url_nuts1)))
#   unzip(file.path(work.dir,basename(url_nuts1)),exdir=file.path(work.dir))
# }
# nuts2<-shapefile(file.path(work.dir,'NUTS_RG_01M_2016_4326_LEVL_2.shp'))

# NUTS_RG_01M_2016_3035_corrected_fixedGeom_EUCROPMAPv7.shp

#eucropmap_hist_nuts2<-shapefile(file.path(work.dir,'nuts2/extractionNuts2/NUTS_RG_01M_2016_3035_corrected_fixedGeom_EUCROPMAPv7.shp'))

eucropmap_hist_nuts2<-shapefile(file.path(work.dir,'/nuts2/NUTS_RG_01M_2016_3035_corrected_fixedGeom_EUCROPMAPv7_masked.shp'))



df_eucropmap<-melt(eucropmap_hist_nuts2@data,
         # ID variables - all the variables to keep but not split apart on
         id.vars=c( "NUTS_ID"),
         # The source columns
         measure.vars=colnames(eucropmap_hist_nuts2@data)[10:(ncol(eucropmap_hist_nuts2@data))],
         # Name of the destination column that will identify the original
         # column that the measurement came from
         variable.name="crop",
         value.name="area_estimated"
)


df_eucropmap$crop<-str_split_fixed(df_eucropmap$crop, "_", 2)[,2]


###########################
# EXTRACTION EN EUCROPMAP DONE IN QGIS LEVEL 0,1,2,3
###########################
# 
# 'https://ec.europa.eu/eurostat/cache/GISCO/distribution/v2/nuts/shp/NUTS_RG_01M_2016_4326.shp.zip'
# { 'COLUMN_PREFIX' : 'HISTO_', 'INPUT_RASTER' : '/data/Dropbox/JRC/LANDSENSE/CASE-STUDY-8-LUCAS-COPERNICUS-CLASSIF/result/v7/S1_classif_v7_HR_clipped_CLC_masked.tif', 'INPUT_VECTOR' : '/home/andrrap/Downloads/NUTS_RG_01M_2016_4326.shp/NUTS_RG_01M_2016_4326.shp', 'OUTPUT' : '/data/Dropbox/JRC/LANDSENSE/CASE-STUDY-8-LUCAS-COPERNICUS-CLASSIF/data/nuts2/extractionNuts2/NUTS123_eucropmapV7_histogram.shp', 'RASTER_BAND' : 1 }

# Load nuts1 delineation (url_nuts1<-'https://ec.europa.eu/eurostat/cache/GISCO/distribution/v2/nuts/shp/NUTS_RG_01M_2016_3035.shp.zip'

nuts_all<-shapefile(file.path(work.dir,'NUTS_RG_01M_2016_3035.shp','NUTS_RG_01M_2016_3035.shp'))

# DE and UK report at the nuts 1 level and need to be removed
nuts2_withoutDE_UK<-nuts_all[(nuts_all$LEVL_CODE==2 & nuts_all$CNTR_CODE!='DE' & nuts_all$CNTR_CODE!='UK'),]
nuts1_withDE_UK<-nuts_all[(nuts_all$LEVL_CODE==1 & (nuts_all$CNTR_CODE=='DE' | nuts_all$CNTR_CODE=='UK')),]
spdf<-rbind(nuts2_withoutDE_UK,nuts1_withDE_UK,makeUniqueIDs = TRUE)

EU28<-c("AT", "BE", "BG",  "CY" ,"CZ" ,"DE" ,"DK", "EE" ,"EL" ,"ES" ,"FI" ,"FR" , "HU","HR", "IE" ,"IT" , "LT" ,"LU", "LV", "MT" ,"NL", "PL", "PT", "RO", "SE", "SI", "SK", "UK")

spdf<-spdf[spdf$CNTR_CODE %in% EU28,]
shapefile(spdf,filename=file.path(work.dir,'NUTS_RG_01M_2016_3035_corrected.shp'))


spdf<-shapefile(file.path(work.dir,'NUTS_RG_01M_2016_3035_corrected.shp'))
###########################
# EXTRACTION STAT ESTAT
###########################
# Extraction via databrowser
# https://ec.europa.eu/eurostat/databrowser/view/APRO_CPSHR__custom_575301/settings_1/table?lang=en
# table downloaded : file:///data/Dropbox/JRC/LANDSENSE/CASE-STUDY-8-LUCAS-COPERNICUS-CLASSIF/data/ComparisonWithEstat/APRO_CPSHR1613486268829.xlsx

estat<-read.csv(file.path(work.dir,'nuts2/ESTAT_2018_AREA_cleaned_noformat.csv'),encoding='utf-8')

# remove C1220 Winter cereal mixtures (maslin)
estat<-subset(estat,select= -c(C1220))
estat$C1500[2:nrow(estat)]<-as.numeric(estat$C1500[2:nrow(estat)])+ as.numeric(estat$G3000[2:nrow(estat)])
estat<-subset(estat,select= -c(G3000))

estat_legend<-data.frame(estat_crop_class=as.character(colnames(estat)[3:length(estat)]),
                         estat_crop_label=as.character(estat[1,3:length(estat)]) )

write.csv(estat_legend,file = file.path(work.dir,'nuts2/estat_legend.csv'))

estat_legend_NUTS2<-data.frame(estat_NUTS2_class=estat$GEO..Codes.,
                               estat_NUTS2_label=estat$GEO..Labels. )
write.csv(estat_legend_NUTS2,file = file.path(work.dir,'nuts2/estat_legend_NUTS2.csv'))

#remove second column
estat<-subset(estat,select= -c(GEO..Labels.))

#remove second line
estat<-estat[2:nrow(estat),]
estat[estat==""]<-NA

# remove empty lines or lines with only 0 or NA
lines_with_only_NA<-rowSums(is.na(estat[,2:ncol(estat)]) | estat[,2:ncol(estat)]==0 )==13
estat<-estat[!lines_with_only_NA,]

#keep only NUTS2 for which we have data in estat
estat<-estat[estat$GEO..Codes. %in% unique(eucropmap_hist_nuts2@data$NUTS_ID) ,]


# read convergence file to tranlstae legend
legend_convergence<-read.csv(file.path(work.dir,'nuts2/estat_legend_convergence_eucropmap.csv'))

estat_colnames_new<-colnames(estat)

for (i_crop in seq(2,length(estat_colnames_new))){
  estat_colnames_new[i_crop]<-legend_convergence$eucropmap_class[legend_convergence$estat_crop_class==estat_colnames_new[i_crop]]
}

colnames(estat)<-estat_colnames_new
colnames(estat)[1]<-'NUTS_ID' 



df_estat<-melt(estat,
         # ID variables - all the variables to keep but not split apart on
         id.vars=c( "NUTS_ID"),
         # The source columns
         measure.vars=colnames(estat)[2:(ncol(estat))],
         # Name of the destination column that will identify the original
         # column that the measurement came from
         variable.name="crop",
         value.name="area_reported"
)


df_estat$index<-rep('area_reported',nrow(df_estat))


# convert pixel in hectare : 1 pixel= 0.01 ha - 
df_eucropmap$area_estimated<-df_eucropmap$area_estimated*(0.01/1000)


df<-merge(df_estat,df_eucropmap,by.x=c('NUTS_ID','crop'),by.y=c('NUTS_ID','crop'))
save(df,file=paste0(work.dir,'nuts2/MergedDfEstatEucropmap.Rdata'))








