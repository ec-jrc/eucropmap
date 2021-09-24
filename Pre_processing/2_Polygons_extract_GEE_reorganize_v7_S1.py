# directories
# JEODPP

# The S1 times series extracted for the COPERNICUS polygons from GEE
data_path='/eos/jeodpp/data/projects/REFOCUS/data/S1_GS/all-10days/extract_GEE_dec2020/'
project_path='/eos/jeodpp/data/projects/REFOCUS/data/S1_GS/all-10days/Map_v7/'
csv_path='/eos/jeodpp/data/projects/REFOCUS/data/polygons/v7/'

path_pol = '/eos/jeodpp/data/projects/REFOCUS/data/polygons/'

#working directory
pwd = project_path

#import modules
import pandas as pd
from pandas import Series,DataFrame
import csv
import numpy as np
import glob
import os

#reorganize as well the polygon csv file for the legend


# Add the strata and the legend from the polygones attributes info

#load polygons
lucas_polygons=pd.read_csv(os.path.join(csv_path,'LUCAS_2018_Copernicus_attributes_cropmap_level1-2.csv'),usecols=['POINT_ID','stratum','level_1','level_2'],dtype={'level_1':int,'level_2':int})

lucas_polygons.info()
print(lucas_polygons.groupby('POINT_ID').apply(min).shape)

#all the extracted polygones in GEE
#pd_lucas= pd.concat(map(pd.read_csv, glob.glob(os.path.join(data_path,'*.csv'))),sort=False)
pd_lucas = pd.concat(map(lambda file: pd.read_csv(file), glob.glob(os.path.join(data_path,'*.csv'))),sort=False)
#remove the stratum and column geo
pd_lucas=pd_lucas.drop(columns=['.geo'])

pd_lucas.head()
print(pd_lucas.columns)
print(pd_lucas.info())


#link the POINT_ID and select the ones available in the polygones - merge in right - the polygons
pd_lucas_b=pd.merge(pd_lucas,lucas_polygons.loc[:,['POINT_ID','stratum','level_2','level_1']],left_on=['POINT_ID'], right_on =['POINT_ID'], how = 'right')#,validate="many_to_one")

#remove nan - not for S2 as everything will disappear
pd_lucas_b=pd_lucas_b.dropna()
print(pd_lucas_b.groupby('POINT_ID').apply(min).shape)

print(pd_lucas_b.groupby('POINT_ID').apply(min).level_1.value_counts())
print(pd_lucas_b.groupby('POINT_ID').apply(min).level_2.value_counts())
#save v7 dataset
pd_lucas_b.to_csv(os.path.join(project_path,'S1_point_allV7_10days_10m_1Jan-31Dec_EU_ratio-db.csv'))