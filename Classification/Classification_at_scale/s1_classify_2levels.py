#!/usr/bin/env python3

from datetime import datetime
import numpy as np
import pickle
import os
from osgeo import gdal, gdal_array, gdalconst

import sys



input_tif = sys.argv[1]
working_dir= sys.argv[2]
model_path_level1 = sys.argv[3]
model_path_level2 = sys.argv[4]


#'/eos/jeodpp/data/projects/REFOCUS/cropclassif/list_rasters_eu_stratum_1.lst'
#input_tif = '/eos/jeodpp/data/projects/REFOCUS/classification/s1_stack_EU_NW210000000000-0000102400.tif'
#without crop
#input_tif = '/eos/jeodpp/data/projects/REFOCUS/data/S1_GS/lucas-astrid/s1_stack_EU_NW220000045056-0000151552.tif'
#with nodata
#input_tif = '/eos/jeodpp/data/projects/REFOCUS/data/S1_GS/lucas-astrid/s1_stack_EU_NW110000069632-0000114688.tif'
#working_dir= '/eos/jeodpp/data/projects/REFOCUS/data/S1_GS/S1_classif_v2_test/'
#model_path_level1 ='/eos/jeodpp/data/projects/REFOCUS/classification/result/models/level1-mask/RFmodel_LUCAS_[1]_level1-mask_all-polygons_janv-jul2018_24042020'
#model_path_level2 ='/eos/jeodpp/data/projects/REFOCUS/classification/result/models/level2-crop/RFmodel_LUCAS_[1,3,4]_level2-crop_all-polygons_janv-jul2018_24042020'
#model_path = '/eos/jeodpp/data/projects/REFOCUS/classification/RFmodel_LUCAS2_janv-jul2018'

#working_dir='/eos/jeodpp/data/projects/REFOCUS/data/S1_GS/S1_classif_v2/'
#output_tif_level1 = os.path.basename(input_tif)[:-4] + '_predict_level1.tif'
output_tif_level2 = os.path.basename(input_tif)[:-4] + '_predict_level2.tif'

decades = 22  # use 22  decades for prediction

print("Date and time start: " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
# Tell GDAL to throw Python exceptions, and register all drivers
gdal.UseExceptions()
gdal.AllRegister()

print('load the model : ' + model_path_level1)
model1 = pickle.load(open(model_path_level1, 'rb'))
print('load the model : ' + model_path_level2)
model2 = pickle.load(open(model_path_level2, 'rb'))

print('Read input : ' + input_tif)
img_ds = gdal.Open(input_tif, gdal.GA_ReadOnly)
RasterXSize=img_ds.RasterXSize
RasterYSize=img_ds.RasterYSize
RasterCount=img_ds.RasterCount
Projection=img_ds.GetProjection()
GeoTransform=img_ds.GetGeoTransform()

print('Reorder bands')
img = np.zeros((RasterYSize,RasterXSize, decades * 2), np.float32)

# for VH
c = img_ds.RasterCount
for b in range(decades):
    img[:, :, b] = img_ds.GetRasterBand(c).ReadAsArray()
    c = c - 2

# for VV
c = img_ds.RasterCount - 1
for b in range(decades, (decades * 2)):
    img[:, :, b] = img_ds.GetRasterBand(c).ReadAsArray()
    c = c - 2

del c, img_ds

# Take our full image, ignore the bands after July and reshape into long 2d array (nrow*ncol, nband)for classification
new_shape = (img.shape[0] * img.shape[1], img.shape[2])
img_as_array = img[:, :, :(decades * 2)].reshape(new_shape)

# replace no data value by 0 -also not ideal way
nodata_mask = np.sum(np.isnan(img_as_array),axis=1)>0
img_as_array[np.isnan(img_as_array)] = 0


# predict for each pixel level 1
print('Predict level 1')
class_prediction1 = model1.predict(img_as_array)
del model1
# Reshape the classification map
#class_prediction1[nodata_mask]=0

#class_prediction1_reshaped = class_prediction1.reshape(img[:, :, 1].shape)

# write classif level 1
#os.chdir(working_dir)
#print('Save output: ' + output_tif_level1)
#gtiff_driver = gdal.GetDriverByName('GTiff')
#out_ds = gtiff_driver.Create(output_tif_level1, RasterXSize,RasterYSize,  1, gdalconst.GDT_Int16,['compress=LZW','tiled=YES'])
#out_ds.SetProjection(Projection)
#out_ds.SetGeoTransform(GeoTransform)

#out_band = out_ds.GetRasterBand(1)
#out_band.WriteArray(class_prediction1_reshaped)

#out_ds = None

# predict for crop pixel level 2
#print('Predict level 2')

class_prediction2=class_prediction1

# if 200 in the tile run the second model
if 200 in np.unique(class_prediction1):
	print('Predict level 2')
	class_prediction2[class_prediction1==200]=model2.predict(img_as_array[class_prediction1==200])


class_prediction2[nodata_mask]=0

del img_as_array, model2, nodata_mask

# Reshape the classification map
class_prediction2_reshaped = class_prediction2.reshape(img[:, :, 1].shape)
del class_prediction1, class_prediction2, img


# write classif level 2
os.chdir(working_dir)
print('Save output: ' + output_tif_level2)
gtiff_driver = gdal.GetDriverByName('GTiff')
out_ds = gtiff_driver.Create(output_tif_level2, RasterXSize,RasterYSize,  1, gdalconst.GDT_Int16,['compress=LZW','tiled=YES'])
out_ds.SetProjection(Projection)
out_ds.SetGeoTransform(GeoTransform)

out_band = out_ds.GetRasterBand(1)
out_band.WriteArray(class_prediction2_reshaped)

out_ds = None

print("Date and time end: " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
