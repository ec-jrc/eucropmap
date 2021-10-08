import rasterio
import numpy as np
from rasterio.plot import reshape_as_raster, reshape_as_image
import pickle
import matplotlib.pyplot as plt
import os

from pixac.pixelaccuracy import *

def save_model_as_pickle(filename, mdl):
    with open(filename, 'wb') as handle:
        pickle.dump(mdl, handle, protocol=4)
    return

def load_model_from_pickle(filename):
    with open(filename, 'rb') as f:
        mdl = pickle.load(f)
    return mdl

# Read the training data
with open('./example/example_data.pickle', 'rb') as handle:
    example_data = pickle.load(handle)

# my_X: features; lbls: reference labels; prds: predicted labels
my_X = example_data['X']
lbls = example_data['labels']
prds = example_data['preds']


# read the image stack to use for inference
with rasterio.open(r'./example/example_image.tif') as src:
    ras = src.read()
    img = reshape_as_image(ras).astype(np.int)
    profile = src.profile

img = np.nan_to_num(img)

# Building the rf forest models
mdl_fn = r'./example/pa_models.pickle'
if not os.path.exists(mdl_fn):
    pixel_accuracy = PixelAccuracy()
    pixel_accuracy.fit(my_X, lbls, prds)
    save_model_as_pickle(mdl_fn, pixel_accuracy)
else:
    pixel_accuracy = load_model_from_pickle(mdl_fn)

# Model inference
print('First case: no map available')
my_pamaps = pixel_accuracy.inference(img)

# save output to raster
profile.update(
    dtype=rasterio.float64,
    count=my_pamaps.shape[0],
    compress='lzw')
with rasterio.open('./example/pixel_based_accuracy_nomask.tif', 'w', **profile) as dst:
    dst.write(my_pamaps)

print('Second case: map available for masking')
with rasterio.open(r'./example/example_map.tif') as src:
    _ras = src.read(1)
    map = _ras.astype(np.int)

# apply inference
my_pamaps = pixel_accuracy.inference(img, map)

profile.update(
    dtype=rasterio.float64,
    count=1,
    compress='lzw')
with rasterio.open('./example/pixel_based_accuracy_mask.tif', 'w', **profile) as dst:
    my_pamaps[my_pamaps==-10000]=np.nan
    dst.write(np.expand_dims(my_pamaps, 0))


def pixel_based_accuracy_class_i(clf):
    pclf = PixelClassifier(clf)
    probs = pclf.image_predict_proba(np.expand_dims(_img, 0))
    return probs



# EOF
