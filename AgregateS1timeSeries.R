
#agreggate the S1_point_allV7_10days_10m_1Jan-31Dec_EU_ratio-db.csv to plot


##############################
### 1 .Aggregate by crop and by stratum to PREPARE FIG 3 ######
##############################
# read only VV and VH
s1<-read.csv(pipe("cut -f5,42-115 -d','/S1_point_allV7_10days_10m_1Jan-31Dec_EU_ratio-db.csv"))

# keep only two main crops
crops_subset<-c(211,212,213,214,215,216,218,221,222,231,232,233)
s1<-s1[s1$level_2 %in% crops_subset,]


df<-melt(s1,
         # ID variables - all the variables to keep but not split apart on
         id.vars=c( "POINT_ID",'stratum' ,'level_2'),
         # The source columns
         measure.vars=colnames(s1)[2:(ncol(s1)-2)],
         # Name of the destination column that will identify the original
         # column that the measurement came from
         variable.name="index",
         value.name="value"
)

df<-separate(df, col=c("index"),into = c("band","date"), sep="_")
df$date<-as.Date(df$date,format="%Y%m%d")

df_summary <- df%>% 
  group_by(date,band,stratum,level_2) %>%
  summarise(
    MU = mean(value),
    SD = sd(value, na.rm = TRUE),
    N = sum(!is.na(value)),
    upper_limit = MU + SD/sqrt(N),
    lower_limit = MU - SD/sqrt(N)
  )

save(df_summary,file=paste0(work.dir,'LUCAS_2018_training_data_GEE_merged_12crops_summary.Rdata'))




##############################
### 2 .Aggregate by crop GROUP  and by stratum to PREPARE FIG S7 ######
##############################
# read only VV and VH
s1<-read.csv(pipe("cut -f5,42-115 -d','S1_point_allV7_10days_10m_1Jan-31Dec_EU_ratio-db.csv"))

# keep only two main crops
crops_subset<-c(211,212,213,214,215,216,218,221,222,231,232)
s1<-s1[s1$level_2 %in% crops_subset,]

#re aggregate by group
s1$crop_group[s1$level_2 %in% c(211,212,213,214,215,218)]<-rep('cereals',sum(s1$level_2 %in% c(211,212,213,214,215,218)))
s1$crop_group[s1$level_2 %in% c(221,222,231)]<-rep('broad leaf crops',sum(s1$level_2 %in% c(221,222,231)))
s1$crop_group[s1$level_2 %in% c(216)]<-rep('maize',sum(s1$level_2 %in% c(216)))
s1$crop_group[s1$level_2 %in% c(232)]<-rep('rapeseed',sum(s1$level_2 %in% c(232)))



df<-melt(s1,
         # ID variables - all the variables to keep but not split apart on
         id.vars=c( "POINT_ID",'stratum' ,'crop_group'),
         # The source columns
         measure.vars=colnames(s1)[2:(ncol(s1)-3)],
         # Name of the destination column that will identify the original
         # column that the measurement came from
         variable.name="index",
         value.name="value"
)

df<-separate(df, col=c("index"),into = c("band","date"), sep="_")
df$date<-as.Date(df$date,format="%Y%m%d")


df_summary <- df%>% 
  group_by(date,band,stratum,crop_group) %>%
  summarise(
    MU = mean(value),
    SD = sd(value, na.rm = TRUE),
    N = sum(!is.na(value)),
    upper_limit = MU + SD/sqrt(N),
    lower_limit = MU - SD/sqrt(N)
  )
save(df_summary,file=paste0(work.dir,'LUCAS_2018_training_data_GEE_merged_4cropgroups_summary.Rdata'))

