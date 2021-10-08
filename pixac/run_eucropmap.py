import rasterio
import os, glob
import numpy as np
from rasterio.plot import reshape_as_raster, reshape_as_image
import pickle
import matplotlib.pyplot as plt
import time
import re

from pixac.pixelaccuracy import *

def save_model_as_pickle(filename, mdl):
    with open(filename, 'wb') as handle:
        pickle.dump(mdl, handle, protocol=4)
    return

def load_model_from_pickle(filename):
    with open(filename, 'rb') as f:
        mdl = pickle.load(f)
    return mdl
 
# choose the stratum
stratum=2

# Read the training data
root_dir = r'/eos/jeodpp/data/projects/REFOCUS/data/S1_GS/v7/pixaccuracy/pixac/'
rdata_fn = os.path.join(root_dir, r'LUCAS_validation_map_v7_ref_S1_1x1_pixac.csv')
mdl_fn=os.path.join(root_dir+'model/', r'pa_models_str'+str(stratum)+'.pickle') 
mp_fns   = glob.glob(os.path.join(root_dir, 's1_stack_EU*level2.tif'))
im_fns   = [x for x in glob.glob(os.path.join(root_dir, 's1_stack_EU*.tif')) if 'level2' not in x]

my_nan = -10000
nnodes = 10

# Read validation data and subset to stratum 1 or 2
df = pd.read_csv(rdata_fn)
df = df[df['stratum'].isin([stratum])]

# classes with 0 correct prediction in stratum 1
if stratum == 1:
    df = df[df['level_2']!=217]
    df = df[df['level_2']!=218]
    df = df[df['level_2']!=219]
    df = df[df['level_2']!=223]
    df = df[df['level_2']!=233]
    df = df[df['level_2']!=600]
elif stratum == 2 :
    df = df[df['level_2']!=214]
    df = df[df['level_2']!=215]
    df = df[df['level_2']!=217]
    df = df[df['level_2']!=218]
    df = df[df['level_2']!=219]
    df = df[df['level_2']!=221]
    df = df[df['level_2']!=222]
    df = df[df['level_2']!=223]

df = df.dropna()

# Generate feature set, labels, and predictions
lbls = df['level_2'].values
prds = df['classification'].values
regex='(((?<![\w\d])VH_)|((?<![\w\d])VV_))(20180[1-7])'
my_X = df[[x for x in list(df.columns) if re.findall(regex,x)]].values
#my_X = df[[x for x in list(df.columns) if x.startswith('V')]].values

# Building the rf forest models
if not os.path.exists(mdl_fn):
    pixel_accuracy = PixelAccuracy()
    pixel_accuracy.fit(my_X, lbls, prds)
    save_model_as_pickle(mdl_fn, pixel_accuracy)
else:
    pixel_accuracy = load_model_from_pickle(mdl_fn)

# Model inference on Sentinel-2 tiles
for mp_fn, im_fn in zip(mp_fns, im_fns):
    start_time = time.time()
    pa_fn = mp_fn.split('.tif')[0]+'_pa.tif'

    # read the image stack to use for inference
    with rasterio.open(im_fn) as src:
        #ras = src.read()
        #select and re-order bands
        VH=list(range(72,28,-2))
        VV=list(range(71,27,-2))
        ras= src.read((*VH,*VV))
        img = reshape_as_image(ras).astype(np.int)
        profile = src.profile

    img = np.nan_to_num(img)

    with rasterio.open(mp_fn) as src:
        _ras = src.read(1)
        map = _ras.astype(np.int)

    # apply inference
    my_pamaps = pixel_accuracy.inference(img, map, my_nan)

    # save output to raster
    profile.update(
        dtype=rasterio.float64,
        count=1,
        compress='lzw')

    with rasterio.open(pa_fn, 'w', **profile) as dst:
        my_pamaps[my_pamaps == my_nan] = np.nan
        dst.write(np.expand_dims(my_pamaps, 0))

    end_time = time.time()
    print(f'Tile processed in {(end_time - start_time)/60} minutes')

# EOF