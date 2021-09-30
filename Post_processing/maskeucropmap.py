from pathlib import Path
from scipy import ndimage
import argparse
import pyjeo as pj

parser=argparse.ArgumentParser()

parser.add_argument("-input","--input",help="input path for crop map",dest="input",required=True,type=str)
parser.add_argument("-output","--output",help="output path for masked crop map",dest="output",required=True,type=str)
parser.add_argument("-ulx","--ulx",help="upper left X coordinate in decimal degrees",dest="ulx",required=True,type=float)
parser.add_argument("-uly","--uly",help="upper left Y coordinate in decimal degrees",dest="uly",required=True,type=float)
parser.add_argument("-lrx","--lrx",help="lower right X coordinate in decimal degrees",dest="lrx",required=True,type=float)
parser.add_argument("-lry","--lry",help="lower right Y coordinate in decimal degrees",dest="lry",required=True,type=float)
parser.add_argument("-maxslope","--maxslope",help="threshold for maximum slope (in degrees)",dest="maxslope",required=False,type=float, default=10)
parser.add_argument("-maxdem","--maxdem",help="threshold for maximum dem (in m)",dest="maxdem",required=False,type=float, default=800)
parser.add_argument("-dem","--dem",help="path of DEM (can be in VRT format)",dest="dem",required=True,type=str)
parser.add_argument("-gsw","--gsw",help="path of Global Surface Water (can be in VRT format)",dest="gsw",required=True,type=str)
parser.add_argument("-clc","--clc",help="path of Corine Land Cover (can be in VRT format)",dest="clc",required=True,type=str)
parser.add_argument("-ghsl","--ghsl",help="path of Global Human Settlement Layer (directory with files in GeoTIFF format)",dest="ghsl",required=True,type=str)
parser.add_argument("-srcnodata","--srcnodata",help="nodata value of eu cropmap",dest="srcnodata",required=False,type=float,default=0)

args = parser.parse_args()

maskclc = [123,124,322, 331, 332, 333, 335, 411, 412, 421, 422, 423]
thresholdgsw = 2# mask values >= 1
thresholdslope = 100*args.maxslope
thresholddem = args.maxdem
thresholdghsl = 249 #250: built-up
#GHSL 2m
# 0 = no data
# 1 = land
# 2 = water
# 255 = built-up area
#GHSL 10m
# 0 = no data
# 1 = land
# 250 = non-residential built-up area
# 255 = residential built-up area

bbox = [args.ulx, args.uly, args.lrx, args.lry]

print("bbox requested: {}".format(bbox))
jimcrop = pj.Jim(args.input, bbox = bbox)
jimcrop.geometry.warp('epsg:3035', dx = 10, dy = 10)

bboxlaea = jimcrop.properties.getBBox()

print("masking slope based on DEM")
jimdem = pj.Jim(args.dem, bbox = bbox)
jimdem.geometry.warp('epsg:3035',
                    ulx = bboxlaea[0],
                    uly = bboxlaea[1],
                    lrx = bboxlaea[2],
                    lry = bboxlaea[3],
                    dx = jimcrop.properties.getDeltaX(),
                    dy = jimcrop.properties.getDeltaY())

jimslope = pj.demops.slope(jimdem)*100.0
jimslope[jimdem<thresholddem] = 0
jimslope.pixops.convert('GDT_UInt16')
jimslope.np()[:] = ndimage.median_filter(jimslope.np(), size = 3)
jimcrop[(jimcrop >0) & (jimslope > thresholdslope)] = args.srcnodata

print("masking global surface water")
jimgsw = pj.Jim(args.gsw, bbox = bbox)
jimgsw.geometry.warp('epsg:3035',
                    ulx = bboxlaea[0],
                    uly = bboxlaea[1],
                    lrx = bboxlaea[2],
                    lry = bboxlaea[3],
                    dx = jimcrop.properties.getDeltaX(),
                    dy = jimcrop.properties.getDeltaY())

jimcrop[(jimcrop >0) & (jimgsw > thresholdgsw)] = args.srcnodata

print("masking global human settlements layer")
for f in args.ghsl.iterdir():
    if f.match("N*.tif"):
        jimghsl = pj.Jim(f, noread = True)
        if jimghsl._jipjim.covers(bboxlaea[0], bboxlaea[1], bboxlaea[2], bboxlaea[3]):
            jimghsl = pj.Jim(f, bbox = jimcrop.properties.getBBox(),
                            dx = jimcrop.properties.getDeltaX(),
                            dy = jimcrop.properties.getDeltaY())
            if jimghsl.properties.getBBox() != jimcrop.properties.getBBox():
                jimghsl.geometry.crop(bbox = jimcrop.properties.getBBox(),
                                    dx = jimcrop.properties.getDeltaX(),
                                    dy = jimcrop.properties.getDeltaY())

            jimcrop[(jimcrop >0) & (jimghsl > thresholdghsl)] = args.srcnodata

print("masking Corine land cover")
jimclc = pj.Jim(args.clc, bbox = bboxlaea,
                dx = jimcrop.properties.getDeltaX(),
                dy = jimcrop.properties.getDeltaY())
for c in maskclc:
    jimcrop[(jimcrop >0) & (jimclc == c)] = args.srcnodata

jimcrop.io.write(args.output,co = ['COMPRESS=LZW', 'TILED=YES'])
