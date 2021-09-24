# Code to generate  FIG 3, 4, 5, 8 , 9 and  10 ; Supp Fig S1, S5, S6. S7
# Data are loaded for each figure but not generated in this script
# For Data generation, see the following script:
# -CREATE LUCAS 2018 VALIDATION POINTS
# -Validation LUCAS POINT
#

require(dplyr)
require(reshape2)
library(tidyr)
require(ggplot2)
library(viridis)
library(xtable)


Sys.setlocale("LC_TIME", "C")

work.dir<-'/data/Dropbox/JRC/LANDSENSE/CASE-STUDY-8-LUCAS-COPERNICUS-CLASSIF/data-published/'
plot.dir<-paste0(work.dir,'figure/')


###########################
##### TABLE 2 ABOUT GSAA : TABLESUMMARY ##### 
###########################

library(tibble)
library(dplyr)
library(xtable)

load(file=file.path(work.dir,'LPIS_ClassGrt001.Rdata'))


df_summary <- df_all%>% 
  group_by(parceltable,crop_code_label) %>%
  summarize(parcel_area = sum(sum),ratio_are=sum(ratio), parcel_count = sum(count))



count_per_crop<-table(df_all$parceltable,as.character(as.integer(as.character(df_all$lucas_code))))

df<-data.frame(cbind(df_summary,count_per_crop,check.names=FALSE))
df$parcel_area<-round(df$parcel_area*0.0001) # SQM to HA
df$parcel_area<-round(df$parcel_area/1000000, digits = 2)
df<-subset(df,select=-( check.names))

df<-data.frame(df)
df$parcel_count<-as.integer(df$parcel_count)
df$ratio_are<-round(100*df$ratio_are)
colnames(df)<-c('Region','Crop code','Area (Mha)','Area Ratio of the total (%)','Parcels (#)',"211", "212","213" ,"214", "215" ,"216" ,"218", "219" ,"221" ,"222" ,"223","231" ,"232" ,"250", "290" ,"300" ,"500")

df<-add_column(df, Description = c('Flanders (Belgium)','Denmark','Centre - Val de Loire (France)','Netherlands','North Rhine-Westphalia (Germany)', 'Slovenia'), .after = "Region")



print(xtable(subset(
  df,select=c('Region','Description','Area (Mha)','Area Ratio of the total (%)','Parcels (#)')
),caption='Summary of LPIS-GSAA data used in the study.'),include.rownames=FALSE)


#####################
# FIG 3 : PLOT THE S1 TIME SERIES ######
#####################
load(paste0(work.dir,'LUCAS_2018_training_data_GEE_merged_12crops_summary.Rdata'))
lucas_legend = read.csv(paste0(work.dir,'legend-lucas7.csv'))
crops_subset<-c(211,212,213,214,215,216,218,221,222,231,232,233)


lucas_legend<-lucas_legend[lucas_legend$class %in% crops_subset,]
lucas_legend<-lucas_legend[order(lucas_legend$class),]

df_summary<-merge(x=df_summary,y=lucas_legend,by.x='level_2',by.y='class')

# New facet label names for supp variable
supp.labs <- c('North EU (Stratum 1)', 'South EU (Stratum 2)')
names(supp.labs) <- c("1", "2")


p3<-ggplot(df_summary[(df_summary$date>as.Date('2018-03-15') & df_summary$date<as.Date('2018-07-31') ), ],
           aes(x=date,y=(MU),fill=factor(level_2),color=factor(level_2)))+ 
  theme(text = element_text(size=7))+ylab('Sentinel-1 Backscatter')+
  scale_fill_manual('',breaks=lucas_legend$class,
                    values=as.character(lucas_legend$colorRGB),
                    labels = as.character(lucas_legend$label))+
  scale_colour_manual('',breaks=lucas_legend$class,
                      values=as.character(lucas_legend$colorRGB),
                      labels = as.character(lucas_legend$label))+
  theme_minimal()+
  theme(legend.position="top")+xlab('Month')+
  scale_x_date(date_breaks = '1 month',date_labels = "%b ")+
  geom_line()+
  geom_ribbon(aes(ymax=MU+3*( SD/sqrt(N)), ymin=MU-3*(SD/sqrt(N)),colour=factor(level_2),fill=factor(level_2)),linetype = 0,alpha=0.2) +
  facet_grid(vars(band), vars(stratum),scales="free_y",labeller = labeller(stratum = supp.labs))




ggsave(paste0('FIG3_','s1_ts_','12crops_smooth.pdf'),plot=p3,
       scale = 1, width = 150, height = 150, units = "mm",
       dpi = 300,path=plot.dir)


#####################
# FIG 4 : PLOT THE OA ######
#####################
library(ggrepel)
manip<-'DATE-BIOME-STRATIFY-CROP_pix'
biome<-'_regroup_remap'
load(paste0(work.dir,manip,biome,'DATE_OA','_summary.Rdata'))


# Accuracy plot
p2<-ggplot(data= accuracy.dates,aes(x=month,y=oa))+geom_point()+
  geom_line()+  theme(legend.position="top")+xlab('Month')+ylab('Overall accuracy [%]')+
  scale_x_date(date_breaks = '1 month',date_labels = "%b ")+ylim(0,100)+  theme(axis.text.x = element_text(angle = 90))+
  geom_label_repel(data =accuracy.dates,aes(label = oa),
                   nudge_x = 0,nudge_y = 3,size=2)+
  theme(text = element_text(size=7))

ggsave(filename=paste0('FIG4_',manip,biome,'DATE_OA.png'),plot = p2,
       scale = 1, width =140, height =  80 ,units = "mm",
       dpi = 300,path=paste0(work.dir,'figure/'))





#####################
# FIG 5 : PLOT THE FSCORE ######
#####################
library(ggrepel)


plot_title<-c('','Stratum 1','Stratum 2')

for ( manip in c('DATE-BIOME-STRATIFY-CROP_pix')) { #'DATE-BIOME-STRATIFY-CROP',
  i<-0
  for (biome in c('_regroup_remap','biome1','biome2')) #biome crop
  { 
    i<-i+1
    #load(paste0(work.dir,'result/',manip,biome,'DATES_F1score','_summary.Rdata'))
    load(paste0(work.dir,manip,biome,'DATE_FSCORE','_summary.Rdata'))
    
    p<-ggplot(data = df.dates[df.dates$X %in% crops_subset,], aes(x = date, y = f1.score,group=as.factor(X),col=as.factor(X)))+
      geom_line() +theme(text = element_text(size=9))+ylab('F1 score')+
      scale_color_manual('',breaks=lucas_legend$class,
                         values=as.character(lucas_legend$colorRGB),
                         labels = as.character(lucas_legend$label))+ 
      theme(legend.position="top")+xlab('Month')+
      theme_minimal()+ylim(0,1)+
      scale_x_date(date_breaks = '1 month',date_labels = "%b ")+
      theme(axis.text.x = element_text(angle = 90))+
      theme(text = element_text(size=6))+ggtitle(plot_title[i])+theme(plot.title = element_text(hjust = 0.5))+
      coord_cartesian(xlim = c(min(df.dates$date), max(df.dates$date) )) +
      geom_label_repel(data = df.dates[df.dates$X %in% crops_subset & df.dates$month==10,],aes(label = label),size=2,
                       na.rm = TRUE,xlim = c(-Inf, Inf), ylim = c(-Inf, Inf),fill = "light grey")+theme(legend.position = "none")
    
    ggsave(filename=paste0('FIG5_',manip,biome,'DATES_F1score_12crops.png'),plot = p,
           scale = 1, width =140, height =  140 ,units = "mm",
           dpi = 300,path=paste0(work.dir,'figure/'))
    
  }}

ggplot(data = df.dates[df.dates$X %in% crops_subset,], aes(x = date, y = f1.score,group=as.factor(X),col=as.factor(X)))+
  geom_line() +theme(text = element_text(size=9))+ylab('F1 score')+
  scale_color_manual('',breaks=lucas_legend$class,
                     values=as.character(lucas_legend$colorRGB),
                     labels = as.character(lucas_legend$label))+ 
  theme(legend.position="top")+xlab('Month')+
  theme_minimal()+
  scale_x_date(date_breaks = '1 month',date_labels = "%b ")+
  theme(axis.text.x = element_text(angle = 90))+
  theme(text = element_text(size=6))+ggtitle(plot_title[i])+theme(plot.title = element_text(hjust = 0.5))+
  coord_cartesian(xlim = c(min(df.dates$date), max(df.dates$date) )) +
  geom_label_repel(data = df.dates[df.dates$X %in% crops_subset & df.dates$month==10,],aes(label = label),
                   na.rm = TRUE,xlim = c(-Inf, Inf), ylim = c(-Inf, Inf),fill = "light grey")+theme(legend.position = "none")



###########################
# FIG 8 - PLOT FSCORE MAP ####
###########################
load(paste0(work.dir,'FSCORE_v7_by_country.Rdata'))

library(reshape2)
library(ggplot2)
library(tidyr)
library(cowplot)

df_fscore<-melt(df,
                # ID variables - all the variables to keep but not split apart on
                id.vars=c( "Country"),
                # The source columns
                measure.vars=colnames(df)[4:(ncol(df))],
                # Name of the destination column that will identify the original
                # column that the measurement came from
                variable.name="class",
                value.name="fscore"
)


crop_to_keep<-c(211,213,216,232,231,222)

df_fscore<-df_fscore[df_fscore$class %in% crop_to_keep,]

lucas_legend = read.csv(paste0(work.dir,'legend-lucas7.csv'))

lucas_legend<-lucas_legend[lucas_legend$class %in% crop_to_keep,]
lucas_legend<-lucas_legend[order(lucas_legend$class),]

# load the nuts 
spdf<-shapefile(file.path(work.dir,'NUTS_RG_01M_2016_3035.shp','NUTS_RG_01M_2016_3035.shp'))
spdf<-spdf[spdf$LEVL_CODE==0,]
spdf<- subset(spdf,select=NUTS_ID)
spdf<-st_as_sf(spdf)
# crop oveer continetal EU
spdf=st_intersection(st_make_valid(spdf),st_transform(st_set_crs(st_as_sf(as(raster::extent(-10.4559,34.0809, 34.5702,70.0841), "SpatialPolygons")), 4326),
                                                      st_crs(spdf)))


df_fscore_reshaped<-reshape2::dcast(df_fscore,Country~class)
colnames(df_fscore_reshaped)[2:ncol(df_fscore_reshaped)]<-paste0('crop_',colnames(df_fscore_reshaped)[2:ncol(df_fscore_reshaped)])


spdf_merge<-st_as_sf( merge(spdf,df_fscore_reshaped,by.x='NUTS_ID',by.y='Country'))



theme_set(cowplot::theme_map() +
            theme(panel.grid.major=element_line(colour="transparent")))
scale_fill_continuous <- function(...) ggplot2::scale_fill_continuous(..., type = "viridis")

spdf_merge2<-spdf_merge %>%
  gather("culture", "fscore", crop_211:crop_231)


# add the label for the crops
spdf_merge2<-spdf_merge2 %>% separate(culture, c("A","crop"), sep = "_")
lucas_legend = read.csv(paste0(work.dir,'legend-lucas7.csv'))
lucas_legend<-lucas_legend[lucas_legend$class %in% unique(spdf_merge2$crop),]
lucas_legend<-lucas_legend[order(lucas_legend$class),]

# reorder the crop label
spdf_merge2<-merge(spdf_merge2,lucas_legend,by.x='crop','class')
spdf_merge2$label_rorder<-factor(spdf_merge2$label, levels=unique(spdf_merge2$label[order(spdf_merge2$crop,spdf_merge2$label)]), ordered=TRUE)

# define the breaks and labels
breaks_diff<-c(0,0.3, 0.4 ,0.5 ,0.6 ,0.7 ,0.8 ,0.9 ,1) 
labels_diff<-c(' 0.3 <' ,'0.3 - 0.4', '0.4 - 0.5' ,'0.5 - 0.6','0.6 - 0.7' ,'0.7 - 0.8','0.8 - 0.9' ,'> 0.9')

# plot
p<-spdf_merge2 %>%
  mutate(value = cut(fscore,
                     breaks = breaks_diff, right = FALSE,
                     labels = labels_diff)) %>%
  ggplot(., aes(fill = value)) +
  geom_sf() +
  facet_wrap(~ label_rorder,nrow=2) + scale_fill_brewer(palette = "YlGn", direction = 1) +
  theme_map() + 
  theme(legend.position = "bottom")+
  guides(fill = guide_legend(nrow = 1))+
  labs(fill="Fscore")


plot.dir<-paste0(work.dir,'figure/')
ggsave(paste0('FIG8_DATES_F1score_6class_per_country_MAPPED.png'),plot=p,
       scale = 1, width = 250, height = 150, units = "mm",
       dpi = 300,path=plot.dir)



###########################
##### FIG 9 : GSAA PA, UA ##### 
###########################

library(ggplot2)

df_all<-load(paste0(work.dir,'LPIS_Extract','GSAA_UA_PA.Rdata'))

df_all<-drop_na(df_all)
lucas_legend = read.csv(paste0(work.dir,'legend-lucas7.csv'))
lucas_legend<-lucas_legend[lucas_legend$class %in% unique(df_all$label_crop),]
lucas_legend<-lucas_legend[order(lucas_legend$class),]


df_all$label_crop<-as.factor(df_all$label_crop)

df_all<-df_all[!df_all$UA==0 ,] 


# FIG 9 a) PA
p<-ggplot(df_all, aes(fill=label_crop, y=PA, x= region, label=round(PA,digits=2) )) + 
  geom_bar(position=position_dodge2(width=1,preserve = "single"), stat="identity")+
  geom_text(position = position_dodge2(width = 0.9,preserve = "single"), size = 1.5, hjust=-0.02,vjust=0.5,angle=90)+
  # geom_text(aes(y=PA,x= region, label= round(PA,digits=2)), position=position_dodge(width=1,preserve = "single"))+
  ylab('PA')+theme_minimal()+ xlab('Region')+
  theme(legend.position="top")+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual('',breaks=lucas_legend$class,
                    values=as.character(lucas_legend$colorRGB),
                    labels = as.character(lucas_legend$label))+
  theme(text = element_text(size=7))



ggsave(filename=paste0('LPIS_PA.png'),plot = p,
       scale = 1, width =150, height =  100 ,units = "mm",
       dpi = 300,path=paste0(work.dir,'figure/'))

# FIG 9 b) UA

p<-ggplot(df_all, aes(fill=label_crop, y=UA, x= region, label=round(UA,digits=2) )) + 
  geom_bar(position=position_dodge2(width=1,preserve = "single"), stat="identity")+
  geom_text(position = position_dodge2(width = 0.9,preserve = "single"), size =1.5, hjust=-0.02,vjust=0.5,angle=90)+
  # geom_text(aes(y=UA,x= region, label= round(UA,digits=2)), position=position_dodge(width=1,preserve = "single"))+
  ylab('UA')+theme_minimal()+ xlab('Region')+
  theme(legend.position="top")+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual('',breaks=lucas_legend$class,
                    values=as.character(lucas_legend$colorRGB),
                    labels = as.character(lucas_legend$label))+
  theme(text = element_text(size=7))



ggsave(filename=paste0('LPIS_UA.png'),plot = p,
       scale = 1, width =150, height =  100 ,units = "mm",
       dpi = 300,path=paste0(work.dir,'figure/'))



###########################
# TABLE S9 ####
###########################

load(paste0(work.dir,'FSCORE_v7_by_country.Rdata'))
crop_to_keep<-c('Country','N','211','213','216','222','231','232')

df_table<-df[,crop_to_keep]
colnames(df_table)<-c('Country','# points','Fscore 211','Fscore 213','Fscore 216','Fscore 222','Fscore 231','Fscore 232')

df_table[,3:ncol(df_table)]<-round(df_table[,3:ncol(df_table)],2)  
df_table$Country<-as.character(df_table$Country)
df_table <- df_table[order(df_table$Country),]

library(xtable)
print(xtable(df_table,caption='F-score by country for the 6 main crops'),include.rownames=FALSE)



###########################
# FIG 10 - Scatter of area reported compared with are estimated ####
###########################
library(ggpubr)

FacetEqualWrap <- ggproto(
  "FacetEqualWrap", FacetWrap,
  
  train_scales = function(self, x_scales, y_scales, layout, data, params) {
    
    # doesn't make sense if there is not an x *and* y scale
    if (is.null(x_scales) || is.null(x_scales)) {
      stop("X and Y scales required for facet_equal_wrap")
    }
    
    # regular training of scales
    ggproto_parent(FacetWrap, self)$train_scales(x_scales, y_scales, layout, data, params)
    
    # switched training of scales (x and y and y on x)
    for (layer_data in data) {
      match_id <- match(layer_data$PANEL, layout$PANEL)
      
      x_vars <- intersect(x_scales[[1]]$aesthetics, names(layer_data))
      y_vars <- intersect(y_scales[[1]]$aesthetics, names(layer_data))
      
      SCALE_X <- layout$SCALE_X[match_id]
      ggplot2:::scale_apply(layer_data, y_vars, "train", SCALE_X, x_scales)
      
      SCALE_Y <- layout$SCALE_Y[match_id]
      ggplot2:::scale_apply(layer_data, x_vars, "train", SCALE_Y, y_scales)
    }
    
  }
)

facet_wrap_equal <- function(...) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)
  
  ggproto(NULL, FacetEqualWrap,
          shrink = facet_super$shrink,
          params = facet_super$params
  )
}



# load the matched data
load(paste0(work.dir,'nuts2/MergedDfEstatEucropmap.Rdata'))


# read convergence file to tranlstae legend
legend_convergence<-read.csv(file.path(work.dir,'nuts2/estat_legend_convergence_eucropmap.csv'))
legend_convergence<-legend_convergence[legend_convergence$estat_crop_class!='G3000',]

# keep only two main crops
crops_toremove<-c('NODA',100, 219, 223, 230, 240, 250, 290, 300, 500, 600)
df<-df[!df$crop %in% crops_toremove,]
# convert pixel in hectare : 1 pixel= 0.01 ha
crops_toremove<-c(215, 217, 218, 223,233)
df<-df[!df$crop %in% crops_toremove,]

df$area_reported<-as.numeric(as.character(df$area_reported))
df$area_estimated<-as.numeric(df$area_estimated)
df<-merge(df,legend_convergence,by.x='crop','eucropmap_class')

table(df$eucropmap_label,df$estat_crop_label)

theme_set(theme_grey())


p<-ggscatter(df,x='area_reported', y='area_estimated', add = "reg.line",add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
             conf.int = TRUE ,
             label = "NUTS_ID", repel = TRUE) +
  #stat_cor(method = "pearson", label.x = 3,label.x.npc='right')+
  stat_cor(aes(label = paste(..r.label.., ..rr.label.., sep = "~`,`~")),label.x = 3,color="red")+
  facet_wrap_equal(~eucropmap_label,nrow=3, scales = "free") +
  geom_abline(slope=1,intercept = 0,linetype='solid',color="black")+
  theme(aspect.ratio = 1)+xlab("Area reported by Eurostrat (x1000 Ha)") + ylab("Area estimated by EUcropmap (x1000 Ha)")


plot.dir<-paste0(work.dir,'figure/')
ggsave(paste0('FIG10_Scatter_EstatVsEucropmapByNuts2_crop_all_withmaizeG3000.png'),plot=p,
       scale = 1, width = 300, height = 300, units = "mm",
       dpi = 300,path=plot.dir)



###########################
##### TABLE S10  ABOUT GSAA : TABLE MATCHING LEGEND ##### 
###########################
df_all<-load(paste0(work.dir,'LPIS_Extract','GSAA_UA_PA.Rdata'))

df_all<-drop_na(df_all)
lucas_legend = read.csv(paste0(work.dir,'legend-lucas7.csv'))
lucas_legend<-lucas_legend[lucas_legend$class %in% unique(df_all$label_crop),]
lucas_legend<-lucas_legend[order(lucas_legend$class),]


df_all$label_crop<-as.factor(df_all$label_crop)

df_all<-df_all[!df_all$UA==0 ,] 
lucas_crop<-c("211", "212","213" ,"214", "215" ,"216" ,"218", "219" ,"221" ,"222" ,"223","231" ,"232" ,"250", "290" ,"300" ,"500")
parceltable_list<-unique( df_all$parceltable)


df_test<-subset(df,select =c("211","212" ,"213" ,"214", "215" ,"216" ,"218", "219" ,"221" ,"222" ,"223","231" ,"232" ,"250", "290" ,"300" ,"500"))

for (parceltable_list_i in seq(1,length(parceltable_list)) ) {
  for (lucas_crop_i in seq(1,length(lucas_crop)) ){
    
    code_list<-c(df_all$crop_code[df_all$lucas_code==lucas_crop[lucas_crop_i] & df_all$parceltable==parceltable_list[parceltable_list_i]])
    
    if (length(code_list)!=0) {
      # code_list<-as.character(list(code_list))
      code_list<-paste(code_list,collapse = ',')
      df_test[parceltable_list_i,lucas_crop_i]<-code_list }
    
  }}

print(xtable(
  cbind(subset(df,select =c('Region','Crop code')),df_test),
  caption='Summary of LPIS-GSAA legend matching data used in the study.'),include.rownames=FALSE)



###########################
#  TABLE S17 ####
###########################

paste(sort(unique(df$NUTS_ID)),collapse = ', ')
library(xtable)
print(xtable(legend_convergence[,2:5],caption='Legend convergence between Eurostat and EU crop map classes'),include.rownames=FALSE)


###########################
# FIG S5 TO MAP MAP DIFF  ####
###########################
library(sf)
library(reshape2)
library(tidyr)
# library(dplyr)
library(magrittr) # viva la %<>% 
library(tidyverse)
library(cowplot)
library(RColorBrewer)


spdf<-shapefile(file.path(work.dir,'NUTS_RG_01M_2016_3035_corrected.shp'))
spdf$area<-as.numeric((st_area(st_as_sf(spdf))*0.0001)/1000)

df$diff<-df$area_reported-df$area_estimated
df1<-merge(df,spdf@data,by.x='NUTS_ID',by.y='NUTS_ID')
df$diff_percent<-round(df1$diff/df1$area,digits=4)*100

df2<-subset(df,select=c(crop, NUTS_ID, diff))

df3<-reshape2::dcast(df2,NUTS_ID~crop)
colnames(df3)[2:ncol(df3)]<-paste0('crop_',colnames(df3)[2:ncol(df3)])

spdf_sf<-st_as_sf(spdf) 
spdf_sf<- subset(spdf_sf,select=NUTS_ID)

# crop oveer continetal EU
spdf_sf=st_intersection(st_make_valid(spdf_sf),st_transform(st_set_crs(st_as_sf(as(raster::extent(-10.4559,34.0809, 34.5702,70.0841), "SpatialPolygons")), 4326),
                                                            st_crs(spdf_sf)
)
)

spdf_merge<-merge(spdf_sf,df3,by.x='NUTS_ID',by.y='NUTS_ID')


theme_set(cowplot::theme_map() +
            theme(panel.grid.major=element_line(colour="transparent")))
scale_fill_continuous <- function(...) ggplot2::scale_fill_continuous(..., type = "viridis")

spdf_merge2<-spdf_merge %>%
  gather("culture", "diff", crop_211:crop_231)

# add the label for the crops
spdf_merge2<-spdf_merge2 %>% separate(culture, c("A","crop"), sep = "_")
lucas_legend = read.csv(paste0(work.dir,'table/legend-lucas7.csv'))
lucas_legend<-lucas_legend[lucas_legend$class %in% unique(spdf_merge2$crop),]
lucas_legend<-lucas_legend[order(lucas_legend$class),]

# reorder the crop label
spdf_merge2<-merge(spdf_merge2,lucas_legend,by.x='crop','class')
spdf_merge2$label_rorder<-factor(spdf_merge2$label, levels=unique(spdf_merge2$label[order(spdf_merge2$crop,spdf_merge2$label)]), ordered=TRUE)

# define the breaks and labels
breaks_diff<-c(-1000,-15, -10 ,-5 ,0 ,5 ,10 ,15 ,1000) 
labels_diff<-c(' -15 % <' ,'-15 % : -10 %', '-10 % : -5 %' ,'-5 % : 0 %','0 : +5 %' ,'+5 % : +10 %','+10 % : +15 %' ,'> 15 %')

# plot
p<-spdf_merge2 %>%
  mutate(value = cut(diff,
                     breaks = breaks_diff, right = FALSE,
                     labels = labels_diff)) %>%
  ggplot(., aes(fill = value)) +
  geom_sf() +
  facet_wrap(~ label_rorder) + scale_fill_brewer(palette = "BrBG", direction = 1) +
  theme_map() + 
  theme(legend.position = "bottom")+
  guides(fill = guide_legend(nrow = 1))+
  labs(fill="Difference [%]")

plot.dir<-paste0(work.dir,'figure/EstatVsEucropmap/')
ggsave(paste0('FIG_S5_Map_EstatVsEucropmapByNuts2_crop_all_withG3000.png'),plot=p,
       scale = 1, width = 300, height = 300, units = "mm",
       dpi = 300,path=plot.dir)

###########################
# FIGURE S6  HIST DIFF  ####
###########################


p<-spdf_merge2 %>%
  mutate(value = cut(diff,
                     breaks = breaks_diff, right = FALSE,
                     labels = labels_diff)) %>%
  ggplot(., aes(fill = value,x=value)) +
  geom_bar(stat='count') +
  facet_wrap(~ label_rorder) + scale_fill_brewer(palette = "BrBG", direction = 1) +
  theme_minimal() + 
  theme(legend.position = "bottom")+
  geom_vline(xintercept =0, linetype="dotted", color = "red", size=1)+
  guides(fill = guide_legend(nrow = 1))+
  labs(fill="Difference [%]")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs( y='Count of NUTS regions [#]')

plot.dir<-paste0(work.dir,'figure/')
ggsave(paste0('FIG_S6_Histogram_EstatVsEucropmapByNuts2_crop_all_withG300.png'),plot=p,
       scale = 1, width = 300, height = 300, units = "mm",
       dpi = 300,path=plot.dir)



#####################
# FIG S7 PLOT THE S1 TIME SERIES for crop groups ######
#####################
load(paste0(work.dir,'/LUCAS_2018_training_data_GEE_merged_4cropgroups_summary.Rdata'))


supp.labs <- c('North EU (Stratum 1)', 'South EU (Stratum 2)')
names(supp.labs) <- c("1", "2")

p5<-ggplot(df_summary[(df_summary$date>as.Date('2018-03-15') & df_summary$date<as.Date('2018-07-31') ), ],
           aes(x=date,y=(MU),fill=factor(crop_group),color=factor(crop_group),linetype =factor(stratum)))+ 
  theme(text = element_text(size=7))+ylab('Sentinel-1 Backscatter')+
  theme(legend.position="top")+xlab('Month')+
  scale_linetype_manual("Stratum",values=c("1"=1,"2"=2))+
  scale_x_date(date_breaks = '1 month',date_labels = "%b ")+
  geom_line()+ scale_colour_discrete(name = "Crop group")+
  #geom_ribbon(aes(ymax=upper_limit, ymin=lower_limit,colour=factor(crop_group),fill=factor(crop_group)),linetype = 0,alpha=0.2) +
  facet_grid(vars(band),scales="free_y",labeller = labeller(stratum = supp.labs))+theme_bw()

plot.dir<-paste0(work.dir,'figure/')
ggsave(paste0('s1_ts_','4crop_group_smooth_together.pdf'),plot=p5,
       scale = 1, width = 150, height = 150, units = "mm",
       dpi = 300,path=plot.dir)

