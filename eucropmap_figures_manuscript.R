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
# FIG 8 - PLOT FSCORE MAP
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






