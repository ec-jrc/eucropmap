# This script extract the values of the EU crop map for the LUCAS core points preselected for the validation and generate F-SCORE per country needed to generate  FIG 8 and Supplementary Table 9

library(raster)
library(caret)
library(dplyr)


# Set  dir ----
work.dir<-'/data/Dropbox/JRC/LANDSENSE/CASE-STUDY-8-LUCAS-COPERNICUS-CLASSIF/data-published/'
url_ftp_eucropmap<-'https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/EUCROPMAP/2018/EUCROPMAP_2018.tif'  # file is 6.3 Gb

###########################
# LOAD VALIDATION POINTS AND EXTRACT THE CLASSIF VALUE ON THOSE POINTS (It could be long!)
###########################
# load the validation points
load(paste0(work.dir,'LUCAS_2018_points_validation.RData'))

# download the EU cropmap raster
options(timeout=1000) # avoid time out error when downloading as the file is big
download.file(url_ftp_eucropmap,destfile=file.path(work.dir,basename(url_ftp_eucropmap))) # file is 6.3 Gb

# read the raster
classif<-raster(file.path(work.dir,'EUCROPMAP_2018.tif'))

# extract the ratser on the point locations
validation_extract<-extract(classif,validation)

# save the extract
validation@data$prediction<-validation_extract
save(validation,file = paste0(work.dir,'s1-classif_v7_lucas-validation_all.RData'))


###########################
# CONVERT LUCAS LC1 TO EUCROPMAP LEGEND LABEL
###########################

# load the legend with label
conversion_LUCAS_eucropmap <-read.csv(paste0(work.dir,'EuroCropMap_v7_LUCAStoEUcropmaplabel.csv'))
df<-merge(validation,conversion_LUCAS_eucropmap, by.x ='LC1',by.y='LC1',all.x = TRUE,sort=FALSE)

# LOOP FOR 'F40U0'->600'  'F40U1'->'290' 
u1<-c('U111','U112','U113')
u0<-c('U120','U420','U415','U361','U362','U350','U370','U330','U150','U411','U412','U140','U312','U341','U314','U414','U210','U319','U318','U222','U321','U322','U413','U311','U317','U342','U313','U224','U226','U223','U316','U315','U225','U221','U227','U130','U228')

df@data$LU1_simplified<-rep('NA', nrow(df))

df@data$LU1_simplified[df@data$LU1 %in% u0]<-rep('u0', length(df[df@data$LU1 %in% u0,]))
df@data$LU1_simplified[df@data$LU1 %in% u1]<-rep('u1', length(df[df@data$LU1 %in% u1,]))

df$level_2[df$LU1_simplified=='u0' & df$LC1=='F40']<-rep(600, length(df$level_2[df$LU1_simplified=='u0' & df$LC1=='F40']))
df$level_2[df$LU1_simplified=='u1' & df$LC1=='F40']<-rep(290, length(df$level_2[df$LU1_simplified=='u1' & df$LC1=='F40']))
validation<-df

###########################
# ADD EXPLICIT LABEL AND FORMAT
###########################

legend_eucropmap<-read.csv(paste0(work.dir,'2020_EuroCropMap_v5.csv'))

df<-merge(validation,legend_eucropmap, by.x ='level_2',by.y='code',all.x = TRUE,sort=FALSE)
validation$level_2_label<-df$label

validation@data$classification<-validation@data$prediction
df<-merge(validation,legend_eucropmap, by.x ='classification',by.y='code',all.x = TRUE,sort=FALSE)
validation$classification_label<-df$label

# Conf Mat with all classes
validation$EUCROPMAP_class<-factor(validation$level_2)
validation$classification<-factor(validation$classification,levels=levels(validation$EUCROPMAP_class))


###########################
# SAVE THE TABLE TO CSV
###########################
write.csv(validation@data[,c("POINT_ID", "LC1" ,"LU1" ,"level_2","level_2_label","classification", "classification_label","stratum","NUTS0","NUTS1", "NUTS2","NUTS3","TH_LAT","TH_LONG")],
          file=paste0(work.dir,'s1-classif_v7_lucas-validation_explicit_label.csv'))


###########################
# CONFUSION MATRIX ALL POINTS
###########################
# cm = table(validation$classification, validation$EUCROPMAP_class,useNA = "ifany")
# cm.with.sum<- cbind(cm,TOTAL=margin.table(cm,margin=1))
# cm.with.sums.1<- rbind(cm.with.sum,TOTAL=margin.table(cm.with.sum,margin=2))
# confusionMatrix(validation$classification, validation$EUCROPMAP_class)
# write.csv(cm.with.sums.1,file=paste0(work.dir,'ConfMatrix_classif_v7_all.csv'))
# 

###########################
# CONFUSION MATRIX ONLY CROPS
###########################

# Conf Mat without 100 and 600 only level 2 classes ----
level2<-c(210,211,212,213,214,215,216,217,218,219,220,221,222, 223,230,231,232,233,240,250,290,300,500)

validation_crop<-validation[validation$EUCROPMAP_class %in% level2,]

validation_crop$EUCROPMAP_class<-factor(validation_crop$EUCROPMAP_class)
validation_crop$classification<-factor(validation_crop$classification,levels=levels(validation_crop$EUCROPMAP_class))



############################
# GENERATE THE DATA FOR FIG 8 and TABLE S9 (FSCORE PER COUNTRY)
############################
# FSCORE per country ----
df<-data.frame()
nuts_all<-unique(validation_crop$NUTS0)
for (i in seq(1,length(nuts_all))){
  
  nuts<-nuts_all[i]
  print(nuts)
  #\nuts<-'SI'
  validation_crop_stratum<-validation_crop[validation_crop$NUTS0 %in% nuts,]
  validation_crop_stratum$EUCROPMAP_class<-factor(validation_crop_stratum$EUCROPMAP_class)
  validation_crop_stratum$classification<-factor(validation_crop_stratum$classification,levels=levels(validation_crop_stratum$EUCROPMAP_class))
  
  cm = table(validation_crop_stratum$classification, validation_crop_stratum$EUCROPMAP_class,useNA = "ifany")
  cm.with.sum<- cbind(cm,TOTAL=margin.table(cm,margin=1))
  cm.with.sums.1<- rbind(cm.with.sum,TOTAL=margin.table(cm.with.sum,margin=2))
  
  cm_caret<-confusionMatrix(validation_crop_stratum$classification, validation_crop_stratum$EUCROPMAP_class)
  print(nuts)
  print(nrow(validation_crop_stratum))
  print(cm_caret$overall[1])
  
  if(i==1){ df <-data.frame(as.character(nuts),
                            nrow(validation_crop_stratum),
                            cm_caret$overall[1],
                            t(cm_caret$byClass[,c('F1')])
                            )}
  if(i>1){ df <-bind_rows(
    df,
    data.frame(as.character(nuts),
               nrow(validation_crop_stratum),
               cm_caret$overall[1],
               t(cm_caret$byClass[,c('F1')])
               ))}
}



colnames(df)<-c('Country','N','Accuracy',as.character(data.frame(strsplit( colnames(df)[4:ncol(df)],"Class.."))[2,]))
save(df,file=paste0(work.dir,'FSCORE_v7_by_country.Rdata'))

# ###########################
# # CONFUSION MATRIX, OA, OA per country , OA per stratum 
# ###########################
# cm = table(validation_crop$classification, validation_crop$EUCROPMAP_class,useNA = "ifany")
# cm.with.sum<- cbind(cm,TOTAL=margin.table(cm,margin=1))
# cm.with.sums.1<- rbind(cm.with.sum,TOTAL=margin.table(cm.with.sum,margin=2))
# 
# # cm_caret<-confusionMatrix(validation_crop$classification, validation_crop$EUCROPMAP_class)
# # write.csv(cm_caret$byClass,file=paste0(work.dir,'data/validation/ConfMatrix_classif_v7_crop_by Class.csv'))
# # write.csv(round(cm_caret$byClass,digits = 4),file=paste0(work.dir,'data/validation/ConfMatrix_classif_v7_crop_by Class_rounded.csv'))
# # 
# # cm_caret$overall
# 
# 
# # OA per country ----
# df<-data.frame()
# nuts_all<-unique(validation_crop$NUTS0)
# for (i in seq(1,length(nuts_all))){
#   
#   nuts<-nuts_all[i]
#   print(nuts)
#   #\nuts<-'SI'
#   validation_crop_stratum<-validation_crop[validation_crop$NUTS0 %in% nuts,]
#   validation_crop_stratum$EUCROPMAP_class<-factor(validation_crop_stratum$EUCROPMAP_class)
#   validation_crop_stratum$classification<-factor(validation_crop_stratum$classification,levels=levels(validation_crop_stratum$EUCROPMAP_class))
#   
#   cm = table(validation_crop_stratum$classification, validation_crop_stratum$EUCROPMAP_class,useNA = "ifany")
#   cm.with.sum<- cbind(cm,TOTAL=margin.table(cm,margin=1))
#   cm.with.sums.1<- rbind(cm.with.sum,TOTAL=margin.table(cm.with.sum,margin=2))
#   
#   cm_caret<-confusionMatrix(validation_crop_stratum$classification, validation_crop_stratum$EUCROPMAP_class)
#   print(nuts)
#   print(nrow(validation_crop_stratum))
#   print(cm_caret$overall[1])
#   if(i==1){ df <-data.frame(as.character(nuts),nrow(validation_crop_stratum),cm_caret$overall[1])}
#   if(i>1){ df <-rbind(df,data.frame(as.character(nuts),nrow(validation_crop_stratum),cm_caret$overall[1]) )}
# }
# 
# colnames(df)<-c('Country','N','Accuracy')
# df$Accuracy<-round(df$Accuracy,2)  
# df$Country<-as.character(df$Country)
# df <- df[order(df$Country),]
# 
# df<-data.frame(df)
# df$N<-as.numeric(df$N)
# df$Accuracy<-as.numeric(df$Accuracy)
# 
# write.csv(df,file=paste0(work.dir,'data/validation/ConfMatrix_classif_v7_crop_by_country.csv'),row.names=F)
# 
# 
# # OA per stratum ----
# df<-data.frame()
# nuts_all<-unique(validation_crop$stratum)
# for (i in seq(1,length(nuts_all))){
#   print(nuts)
#   nuts<-nuts_all[i]
#   
#   validation_crop_stratum<-validation_crop[validation_crop$stratum %in% nuts,]
#   validation_crop_stratum$EUCROPMAP_class<-factor(validation_crop_stratum$EUCROPMAP_class)
#   validation_crop_stratum$classification<-factor(validation_crop_stratum$classification,levels=levels(validation_crop_stratum$EUCROPMAP_class))
#   
#   cm = table(validation_crop_stratum$classification, validation_crop_stratum$EUCROPMAP_class,useNA = "ifany")
#   cm.with.sum<- cbind(cm,TOTAL=margin.table(cm,margin=1))
#   cm.with.sums.1<- rbind(cm.with.sum,TOTAL=margin.table(cm.with.sum,margin=2))
#   
#   cm_caret<-confusionMatrix(validation_crop_stratum$classification, validation_crop_stratum$EUCROPMAP_class)
#   print(nuts)
#   print(nrow(validation_crop_stratum))
#   print(cm_caret$overall[1])
#   if(i==1){ df <-data.frame(as.character(nuts),nrow(validation_crop_stratum),cm_caret$overall[1])}
#   if(i>1){ df <-rbind(df,data.frame(as.character(nuts),nrow(validation_crop_stratum),cm_caret$overall[1]) )}
# }
# 
# 
# colnames(df)<-c('Stratum','N','Accuracy')
# df$Accuracy<-round(df$Accuracy,2)  
# newdata <- df[order(df$Stratum),]
# 
# df<-data.frame(df)
# write.csv(df,file=paste0(work.dir,'data/validation/ConfMatrix_classif_v7_crop_by_stratum.csv'),row.names=F)




