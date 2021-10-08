import rasterio
import os, glob
import numpy as np
from rasterio.plot import reshape_as_raster, reshape_as_image
import pickle
import matplotlib.pyplot as plt
import time
import re
import sys
import pandas as pd

from pixac.pixelaccuracy import *

def save_model_as_pickle(filename, mdl):
    with open(filename, 'wb') as handle:
        pickle.dump(mdl, handle, protocol=4)
    return

def load_model_from_pickle(filename):
    with open(filename, 'rb') as f:
        mdl = pickle.load(f)
    return mdl

input_tif = sys.argv[1]
working_dir= sys.argv[2]
classif_dir= sys.argv[3]
model_path = sys.argv[4]
input_csv= sys.argv[5]
stratum = sys.argv[6]

#sys.path.append('/eos/jeodpp/home/users/verheas/pixac/')

# Read the training data
#root_dir = r'/eos/jeodpp/home/users/verheas/pixac'
#rdata_fn = os.path.join(root_dir, r'LUCAS_validation_map_v7_ref_S1_1x1.csv')
#mdl_fn = os.path.join(root_dir, r'pa_models.pickle')
#mp_fns   = glob.glob(os.path.join(root_dir, 's1_stack_EU*level2.tif'))
#im_fns   = [x for x in glob.glob(os.path.join(root_dir, 's1_stack_EU*.tif')) if 'level2' not in x]
#mdl_fn = os.path.join(model_path, r'pa_models_str'+str(stratum)+'.pickle')



# data imported for the rest of the script
mdl_fn = os.path.join(model_path, r'pa_models_str'+str(stratum)+'.pickle')
rdata_fn = input_csv
im_fns=input_tif
mp_fns=os.path.join(classif_dir, os.path.basename(input_tif)[:-4] + '_predict_level2.tif')

my_nan = -10000
nnodes = 10

test_dir='/eos/jeodpp/data/projects/REFOCUS/data/S1_GS/pixac_v7_stratum_'+stratum+'/'
pa_fn_test = os.path.join(test_dir, os.path.basename(mp_fns)[:-4]+'_pa.tif')

if not os.path.exists(pa_fn_test):
    print("process"+pa_fn_test)

    # Read validation data and subset to stratum 1
    df = pd.read_csv(rdata_fn,engine='python')
    df = df[df['stratum'].isin([stratum])]
    print(stratum)
    df = df.dropna()

    # classes with 0 correct prediction in stratum 1 and 2
    if stratum == 1:
        df = df[df['level_2']!=217]
        df = df[df['level_2']!=218]
        df = df[df['level_2']!=219]
        df = df[df['level_2']!=223]
        df = df[df['level_2']!=233]
        df = df[df['level_2']!=600]
    elif stratum == 2:
        df = df[df['level_2']!=214]
        df = df[df['level_2']!=215]
        df = df[df['level_2']!=217]
        df = df[df['level_2']!=218]
        df = df[df['level_2']!=219]
        df = df[df['level_2']!=221]
        df = df[df['level_2']!=222]
        df = df[df['level_2']!=223]

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
        print('create model')
    else:
        pixel_accuracy = load_model_from_pickle(mdl_fn)
        print('load model')

    # Model inference on Sentinel-2 tiles
    # for mp_fn, im_fn in zip(mp_fns, im_fns):
    start_time = time.time()
    # pa_fn = mp_fn.split('.tif')[0]+'_pa.tif'
    pa_fn = os.path.join(working_dir, os.path.basename(mp_fns)[:-4]+'_pa.tif')

    # read the image stack to use for inference
    with rasterio.open(im_fns) as src:
        #ras = src.read()
        #select and re-order bands
        VH=list(range(72,28,-2))
        VV=list(range(71,27,-2))
        ras= src.read((*VH,*VV))
        img = reshape_as_image(ras).astype(np.int)
        profile = src.profile

    img = np.nan_to_num(img)

    with rasterio.open(mp_fns) as src:
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
else:
    print(pa_fn_test+" already exists")
# EOF