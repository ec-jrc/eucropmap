from pathlib import Path
from scipy import ndimage
import math
import argparse
import osgeo
import osgeo.osr
import osgeo.ogr
import pyjeo as pj

def transformBBox(bbox, s_srs=4326, t_srs=None):
    # print("*** Function bb2wkt was called")
    bbox_copy = bbox.copy()
    ulx = bbox_copy[0]
    uly = bbox_copy[1]
    lrx = bbox_copy[2]
    lry = bbox_copy[3]
    if t_srs:
        # create coordinate transformation
        inSpatialRef = osgeo.osr.SpatialReference()
        outSpatialRef = osgeo.osr.SpatialReference()
        if int(osgeo.__version__[0]) >= 3:
            # GDAL 3 changes axis order: https://github.com/OSGeo/gdal/issues/1546
            inSpatialRef.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
            outSpatialRef.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
        if isinstance(s_srs, int):
            inSpatialRef.ImportFromEPSG(s_srs)
        else:
            inSpatialRef.ImportFromWkt(s_srs)
        if isinstance(t_srs, int):
            outSpatialRef.ImportFromEPSG(t_srs)
        else:
            outSpatialRef.ImportFromWkt(t_srs)
        coordTransform = osgeo.osr.CoordinateTransformation(inSpatialRef, outSpatialRef)
        # transform points
        point = osgeo.ogr.Geometry(osgeo.ogr.wkbPoint)
        point.AddPoint(ulx, uly)
        point.Transform(coordTransform)
        ulx = point.GetX()
        uly = point.GetY()
        point = osgeo.ogr.Geometry(osgeo.ogr.wkbPoint)
        point.AddPoint(lrx, lry)
        point.Transform(coordTransform)
        lrx = point.GetX()
        lry = point.GetY()
    if ulx > lrx:
        print("Error: bounding box not correct (ulx>lrx)")
        return (None)
    if uly < lry:
        print("Error: bounding box not correct (uly<lry)")
        return (None)
    return [ulx,uly,lrx,lry]


parser=argparse.ArgumentParser()

parser.add_argument("-output","--output",help="output path for masked crop map",dest="output",required=True,type=str)
parser.add_argument("-outputdem","--outputdem",help="outputdem path for masked crop map",dest="outputdem",required=False,type=str, default=None)
parser.add_argument("-outputslope","--outputslope",help="outputslope path for masked crop map",dest="outputslope",required=False,type=str, default = None)
parser.add_argument("-ulx","--ulx",help="upper left X coordinate in decimal degrees",dest="ulx",required=True,type=float)
parser.add_argument("-uly","--uly",help="upper left Y coordinate in decimal degrees",dest="uly",required=True,type=float)
parser.add_argument("-lrx","--lrx",help="lower right X coordinate in decimal degrees",dest="lrx",required=True,type=float)
parser.add_argument("-lry","--lry",help="lower right Y coordinate in decimal degrees",dest="lry",required=True,type=float)
parser.add_argument("-maxslope","--maxslope",help="threshold for maximum slope (in degrees)",dest="maxslope",required=False,type=float, default=10)
parser.add_argument("-maxdem","--maxdem",help="threshold for maximum dem (in m)",dest="maxdem",required=False,type=float, default=800)
parser.add_argument("-noghsl","--noghsl",help="do not mask global human settlements",dest="noghsl",required=False,type=bool, default=False)
parser.add_argument("-noslope","--noslope",help="do not mask slopes",dest="noslope",required=False,type=bool, default=False)
parser.add_argument("-dem","--dem",help="base on dem 30 m or 10 m",dest="dem",required=False,type=int, default=30)
parser.add_argument("-noclc","--noclc",help="do not mask corine land cover",dest="noclc",required=False,type=bool, default=False)
parser.add_argument("-nogsw","--nogsw",help="do not mask global surface water",dest="nogsw",required=False,type=bool,default=False)
parser.add_argument("-srcnodata","--srcnodata",help="nodata value of eu cropmap",dest="srcnodata",required=False,type=float,default=0)

args = parser.parse_args()

eucropmapfn = Path('/eos/jeodpp/data/projects/REFOCUS/data/S1_GS/S1_classif_v7/S1_classif_v7_HR_clipped.tif')
slopefn = Path('/eos/jeodpp/data/base/Elevation/GLOBAL/AW3D30/VER2-1/Data/GeoTIFF/slopeUtmZones')
dem30mfn = Path('/eos/jeodpp/data/base/Elevation/GLOBAL/AW3D30/VER2-1/Data/GeoTIFF/tiles10deg')
dem10mfn = Path('/eos/jeodpp/data/base/Elevation/GLOBAL/COP-DEM/RESTRICTED/GLO-30-DGED/VER1-0/Data/VRT/Copernicus_DSM_10.vrt')
gswfn = Path('/eos/jeodpp/data/base/Hydrography/GLOBAL/JRC/GlobalSurfaceWater/YearlyClassification/VER3-0/Data/GeoTIFF/yearlyClassification2019.tif')
ghsl10mfn = Path('/eos/jeodpp/data/products/Landcover/EUROPE/GHSL/BuiltUp/DWH2-VHR2015/VER-2019_1-0/Data/EPSG-3035/10m/GeoTIFF')
clcfn = Path('/eos/jeodpp/data/base/Landcover/EUROPE/CorineLandCover/CLC2018/VER20-b2/Data/GeoTIFF/100m/clc2018_Version_20_b2.tif')

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

nodataclc = -32768
nodataslope = -9999
nodatagsw = 0

bbox = [args.ulx, args.uly, args.lrx, args.lry]

print("bbox requested: {}".format(bbox))
jimcrop = pj.Jim(eucropmapfn, bbox = bbox)
jimcrop.geometry.warp('epsg:3035', dx = 10, dy = 10)

bboxlaea = jimcrop.properties.getBBox()

if not args.noslope:
    if args.dem == 10:
        print("masking slope based on 10 m DEM")
        jimdem = pj.Jim(dem10mfn, bbox = bbox)
        jimdem.geometry.warp('epsg:3035',
                            ulx = bboxlaea[0],
                            uly = bboxlaea[1],
                            lrx = bboxlaea[2],
                            lry = bboxlaea[3],
                            dx = jimcrop.properties.getDeltaX(),
                            dy = jimcrop.properties.getDeltaY())
        if args.outputdem is not None:
            jimdem.io.write(args.outputdem, co = ['COMPRESS=LZW', 'TILED=YES'])
        jimslope = pj.demops.slope(jimdem)*100.0
        jimslope[jimdem<thresholddem] = 0
        jimslope.pixops.convert('GDT_UInt16')
        jimslope.np()[:] = ndimage.median_filter(jimslope.np(), size = 3)
        if args.outputslope is not None:
            jimslope.io.write(args.outputslope, co = ['COMPRESS=LZW', 'TILED=YES'])
        if not args.noslope:
            jimcrop[(jimcrop >0) & (jimslope > thresholdslope)] = args.srcnodata
    else:# args.dem == 30
        mosaicdem = None
        mosaicslope = None
        print("masking slope based on 30 m DEM")

        lonulx = 10 * (int(bbox[0]) // 10)
        lonlrx = 10 * (int(bbox[2]) // 10)
        latuly = - 10 * (int(-bbox[1]) // 10)
        latlry = - 10 * (int(-bbox[3]) // 10)
        for lon in range(lonulx,lonlrx + 1):
            for lat in range(latlry,latuly + 1):
                if lon < 0:
                    orient = 'W'
                else:
                    orient = 'E'
                if lat < 0:
                    hemisphere = 'S'
                else:
                    hemisphere = 'N'
                f = dem30mfn / Path('aw3d30_dsm_' + orient + str(abs(lon)).zfill(3) + hemisphere + str(abs(lat)).zfill(3)).with_suffix('.tif')
                print('filename is: {}'.format(f))
                jimdem = pj.Jim(f, bbox = bbox)
                jimdem.geometry.warp('epsg:3035',
                                     ulx = bboxlaea[0],
                                     uly = bboxlaea[1],
                                     lrx = bboxlaea[2],
                                     lry = bboxlaea[3],
                                     dx = jimcrop.properties.getDeltaX(),
                                     dy = jimcrop.properties.getDeltaY())
                if mosaicdem is None:
                    mosaicdem = pj.geometry.crop(
                        jimdem,
                        bbox = jimcrop.properties.getBBox(),
                        dx = jimcrop.properties.getDeltaX(),
                        dy = jimcrop.properties.getDeltaY())
                else:
                    mosaicdem.geometry.stackPlane(
                        pj.geometry.crop(jimdem,
                                         bbox = jimcrop.properties.getBBox(),
                                         dx = jimcrop.properties.getDeltaX(),
                                         dy = jimcrop.properties.getDeltaY()))
                    mosaicdem.geometry.reducePlane(nodata = -32768, refband = 0)
        if args.outputdem is not None:
            mosaicdem.io.write(args.outputdem, co = ['COMPRESS=LZW', 'TILED=YES'])
        mosaicdem[mosaicdem < 0] = 0
        jimslope = pj.demops.slope(mosaicdem)*100.0
        jimslope[mosaicdem<thresholddem] = 0
        jimslope.pixops.convert('GDT_UInt16')
        jimslope.np()[:] = ndimage.median_filter(jimslope.np(), size = 9)
        if args.outputslope is not None:
            jimslope.io.write(args.outputslope, co = ['COMPRESS=LZW', 'TILED=YES'])
        if not args.noslope:
            jimcrop[(jimcrop >0) & (jimslope > thresholdslope)] = args.srcnodata

        # utmzone_left = math.ceil((bbox[0]+180)/6)
        # utmzone_right = math.ceil((bbox[2]+180)/6)
        # for utmzone in range(utmzone_left,utmzone_right+1):
        #     f = slopefn / str("slope_utm" + str(utmzone) + "N.tif")
        # utmzone_left = math.ceil((bbox[0]+180)/6)
        # utmzone_right = math.ceil((bbox[2]+180)/6)
        # for utmzone in range(utmzone_left,utmzone_right+1):
        #     f = slopefn / str("slope_utm" + str(utmzone) + "N.tif")
        #     jimslope = pj.Jim(f, noread = True)
        #     bboxutm = transformBBox(bbox, 4326, jimslope.properties.getProjection())
        #     if jimslope._jipjim.covers(bboxutm[0], bboxutm[1], bboxutm[2], bboxutm[3]):
        #         jimslope = pj.Jim(f, bbox = bboxutm)
        #         # print("smoothnodata")
        #         # jimslope.ngbops.filter2d('smoothnodata', nodata = 0)
        #         # print("median filter")
        #         # jimslope.ngbops.filter2d('median',nodata = 0)
        #         jimslope.geometry.warp('epsg:3035',
        #                                 ulx = bboxlaea[0],
        #                                 uly = bboxlaea[1],
        #                                 lrx = bboxlaea[2],
        #                                 lry = bboxlaea[3],
        #                                 dx = jimcrop.properties.getDeltaX(),
        #                                 dy = jimcrop.properties.getDeltaY())

        #         if mosaicslope is None:
        #             mosaicslope = pj.geometry.crop(
        #                 jimslope,
        #                 bbox = jimcrop.properties.getBBox(),
        #                 dx = jimcrop.properties.getDeltaX(),
        #                 dy = jimcrop.properties.getDeltaY())
        #         else:
        #             mosaicslope.geometry.stackPlane(
        #                 pj.geometry.crop(jimslope,
        #                                 bbox = jimcrop.properties.getBBox(),
        #                                 dx = jimcrop.properties.getDeltaX(),
        #                                 dy = jimcrop.properties.getDeltaY()))
        #             mosaicslope.geometry.reducePlane(nodata = 0, refband = 0)
        #         jimcrop[(jimcrop >0) & (jimslope > thresholdslope)] = args.srcnodata

        # if args.outputslope is not None:
        #     jimslope.io.write(args.outputslope, co = ['COMPRESS=LZW', 'TILED=YES'])

#mask global surface water
if not args.nogsw:
    print("masking global surface water")
    jimgsw = pj.Jim(gswfn, bbox = bbox)
    jimgsw.geometry.warp('epsg:3035',
                        ulx = bboxlaea[0],
                        uly = bboxlaea[1],
                        lrx = bboxlaea[2],
                        lry = bboxlaea[3],
                        dx = jimcrop.properties.getDeltaX(),
                        dy = jimcrop.properties.getDeltaY())

    jimcrop[(jimcrop >0) & (jimgsw > thresholdgsw)] = args.srcnodata

#mask human settlements
if not args.noghsl:
    print("masking global human settlements layer")
    for f in ghsl10mfn.iterdir():
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

#mask Corine classes
if not args.noghsl:
    print("masking Corine land cover")
    jimclc = pj.Jim(clcfn, bbox = bboxlaea,
                    dx = jimcrop.properties.getDeltaX(),
                    dy = jimcrop.properties.getDeltaY())
    for c in maskclc:
        jimcrop[(jimcrop >0) & (jimclc == c)] = args.srcnodata

jimcrop.io.write(args.output,co = ['COMPRESS=LZW', 'TILED=YES'])