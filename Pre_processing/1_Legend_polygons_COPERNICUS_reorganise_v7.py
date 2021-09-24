# Preparing the data for the classificaiton - LUCAS polygons

# 1) Selecting the polygons used according to the legend 
# 2) Recoding the legend of the LUCAS COPERNICUS POLYGONS to fit the EU CROP MAP Legend

# # # # # # # # # # # # # 
# directories JEODPP# # # 
# # # # # # # # # # # # # 

project_path='/eos/jeodpp/data/projects/REFOCUS/data/polygons/v7/'
path_pol = '/eos/jeodpp/data/projects/REFOCUS/data/polygons/'
 
#import modules
import pandas as pd
from pandas import Series,DataFrame
import csv
import numpy as np
import time
import sklearn
import scipy
import matplotlib.pyplot as plt
import os
from datetime import datetime
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix

import glob
import os

#####################################
##########prepare the data ##########
#####################################

#1)load csv with the COPERNIUS polygons with the strata information
lucas_polygons = pd.read_csv(os.path.join(path_pol,'LUCAS_2018_Copernicus_polygons_attributes_strata.csv'))

print(lucas_polygons.shape)
print(lucas_polygons.head())

#2) drop the classes not relevant for the EU crop map
drop=lucas_polygons[lucas_polygons.LC1.isin(['A30','H11','H12','H21','H22','H23','G21','G30','G11'])]
print('number of dropped polygons',drop.shape)
                        
lucas_polygons=lucas_polygons[~lucas_polygons.LC1.isin(['A30','H11','H12','H21','H22','H23','G21','G30','G11'])]

#3) Reclassify the land use classes of LUCAS in two categories

lucas_polygons['LU1s']=lucas_polygons['LU1']

lucas_polygons['LU1s']=lucas_polygons['LU1s'].replace(['U111','U112','U113'],'U1')
lucas_polygons['LU1s']=lucas_polygons['LU1s'].replace(['U120','U420','U415','U361','U362','U350','U370','U330','U150','U411','U412','U140','U312','U341','U314','U414','U210'
                                           ,'U319','U318','U222','U321','U322','U413','U311','U317','U342','U313','U224','U226','U223','U316','U315','U225','U221','U227'],'U0')

print(lucas_polygons['LU1s'].value_counts())

#joint field - #create a field for the legend mixing LC1 and LU1
lucas_polygons['LC1_LU1']=lucas_polygons['LC1'].astype(str) + lucas_polygons['LU1s'].astype(str)

###########################
##########level_2##########
###########################

lucas_polygons['level_2']=lucas_polygons['LC1_LU1']

#reclassify in 100
lucas_polygons.head()
lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace='A[12][0123]U[01]',value='100',regex=True)

#croplands
lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace=('B11U[01]'),value='211',regex=True)
lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace=('B12U[01]'),value='212',regex=True)
lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace=('B13U[01]'),value='213',regex=True)
lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace=('B14U[01]'),value='214',regex=True)
lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace=('B15U[01]'),value='215',regex=True)
lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace=('B16U[01]'),value='216',regex=True)
lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace=('B17U[01]'),value='217',regex=True)
lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace=('B18U[01]'),value='218',regex=True)
lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace=('B19U[01]'),value='219',regex=True)

lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace=('B21U[01]'),value='221',regex=True)
lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace=('B22U[01]'),value='222',regex=True)
lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace=('B23U[01]'),value='223',regex=True)

lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace=('B3[4-7]U[01]'),value='230',regex=True)
lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace=('B31U[01]'),value='231',regex=True)
lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace=('B32U[01]'),value='232',regex=True)
lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace=('B33U[01]'),value='233',regex=True)

lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace=('B4[1-5]U[01]'),value='240',regex=True)
lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace=('B5[1-4]U[01]'),value='250',regex=True)
lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace=('B55U[01]'),value='500',regex=True)

lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace=('B7[1-7]U[01]'),value='300',regex=True)
lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace=('B8[1-4]U[01]'),value='300',regex=True)

#woodland - no land use choice
lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace='C[123][0123]U[01]',value='300',regex=True)

#shrubland - no land use choice
lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace='D[123][0123]U[01]',value='300',regex=True)

#grassland - 

#agricultural grasslands
lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace='E[123][0123]U[1]',value='500',regex=True)
#non agricultural grasslands
lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace='E[1][0123]U[0]',value='500',regex=True)
lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace='E[2][0123]U[0]',value='500',regex=True)
lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace='E[3][0123]U[0]',value='500',regex=True)

#agricultural bare lands
lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace='F40U1',value='290',regex=True)
#non agricultural grasslands
lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace='F[123][0]U[01]',value='600',regex=True)
lucas_polygons['level_2']=lucas_polygons['level_2'].replace(to_replace='F40U0',value='600',regex=True)

print(lucas_polygons['level_2'].value_counts())

###########################
##########level_1##########
###########################

lucas_polygons['level_1']=lucas_polygons['level_2']
#pixels counts per class
#print(lucas_polygons['level_1'].value_counts())

#arable land + grassland
lucas_polygons['level_1']=lucas_polygons['level_1'].replace(to_replace=('2[1-6][0-9]'),value='200',regex=True)

#permament crops
lucas_polygons['level_1']=lucas_polygons['level_1'].replace(to_replace=('2[7][0-9]'),value='300',regex=True)
lucas_polygons['level_1']=lucas_polygons['level_1'].replace(to_replace=('2[8][0-9]'),value='300',regex=True)

#reclassify agricultural bareland in agricultural 200
lucas_polygons['level_1']=lucas_polygons['level_1'].replace(to_replace=('290'),value='200',regex=True)

#pixels counts per class
print(lucas_polygons['level_1'].value_counts())

#drop columns
lucas_polygons=lucas_polygons.drop(columns=['LC1_LU1','LU1s','stratum.1'])
print(lucas_polygons.shape)

################################
##########save as csv ##########
################################

lucas_polygons.to_csv(os.path.join(project_path,'LUCAS_2018_Copernicus_attributes_cropmap_level1-2.csv'))
lucas_polygons.head()