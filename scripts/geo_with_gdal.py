#!/usr/bin/env python3

import os
import sys
import gdal
import shutil
import pickle
import argparse
import xml.etree.ElementTree as ET

import isce
import isceobj
from isceobj.Image import createImage

from crlpac import runCmd

def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser( description='geocode images with gdal')
    parser.add_argument('-input', dest='input', type=str, required=True,
            help = 'input file to be geocoded')
    parser.add_argument('-output', dest='output', type=str, required=True,
            help = 'output file to be geocoded')
    parser.add_argument('-lat', dest='lat', type=str, required=True,
            help = 'latitude file.')
    parser.add_argument('-lon', dest='lon', type=str, required=True,
            help = 'longitude file')
    parser.add_argument('-bbox', dest='bbox', type=str, required=True,
            help='geocode bounding box (format: S/N/W/E). or you can input master frame pickle file that contains a bounding box')
    parser.add_argument('-ssize', dest='ssize', type=float, default=1.0,
            help = 'output sample size. default: 1.0 arcsec')
    parser.add_argument('-rmethod', dest='rmethod', type=int, default=0,
            help = 'resampling method. 0: near (default), 1: bilinear, 2: cubic, 3: cubicspline, 4: lanczos, 5: average, 6: mode, 7: max, 8: min, 9: med, 10: Q1, 11: Q3.')
    parser.add_argument('-topshift', dest='topshift', type=int, default=0,
            help = 'top shift of image to be geocoded. default: 0')
    parser.add_argument('-leftshift', dest='leftshift', type=int, default=0,
            help = 'left shift of image to be geocoded. default: 0')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':

    inps = cmdLineParse()

    #get bounding box
    if os.path.isfile(inps.bbox):
        with open(inps.bbox, 'rb') as f:
            frame = pickle.load(f)
        bbox = frame.snwe
    else:
        bbox = [float(val) for val in inps.bbox.split('/')]
        if len(bbox) != 4:
            raise Exception('bbox should contain 4 floating point values!')
    print("geocode bounding box:")
    print("south: {}".format(bbox[0]))
    print("north: {}".format(bbox[1]))
    print("west: {}".format(bbox[2]))
    print("east: {}".format(bbox[3]))


    #afer comparing with geo_with_ll.py, I find that it's better to
    inps.topshift += -1

    #many things only work in the current directory because of many levels of function calls
    #I do the following to overcome that
##################################################################################################
    input_cd = 0
    lat_cd = 0
    lon_cd = 0
    if not os.path.isfile(os.path.basename(inps.input)):
        input_cd = 1
        os.symlink(inps.input, os.path.basename(inps.input))
        shutil.copyfile(inps.input+'.xml', os.path.basename(inps.input)+'.xml')
        inps.input = os.path.basename(inps.input)

    #if not os.path.isfile(os.path.basename(inps.output)):
    #    output_cd = 1
    #    inps.output = os.path.basename(inps.output)
    #    os.symlink(inps.output, os.path.basename(inps.output))
    #    shutil.copyfile(inps.output+'.xml', os.path.basename(inps.output)+'.xml')
    if not os.path.isfile(os.path.basename(inps.lat)):
        lat_cd = 1
        os.symlink(inps.lat, os.path.basename(inps.lat))
        shutil.copyfile(inps.lat+'.xml', os.path.basename(inps.lat)+'.xml')
        inps.lat = os.path.basename(inps.lat)
    if not os.path.isfile(os.path.basename(inps.lon)):
        lon_cd = 1
        os.symlink(inps.lon, os.path.basename(inps.lon))
        shutil.copyfile(inps.lon+'.xml', os.path.basename(inps.lon)+'.xml')
        inps.lon = os.path.basename(inps.lon)
##################################################################################################

    cmd = 'isce2gis.py vrt -i ' + inps.input
    runCmd(cmd)
    cmd = 'isce2gis.py vrt -i ' + inps.lat
    runCmd(cmd)
    cmd = 'isce2gis.py vrt -i ' + inps.lon
    runCmd(cmd)

    inImage = createImage()
    inImage.load(inps.input+'.xml')
    latImage = createImage()
    latImage.load(inps.lat+'.xml')
    lonImage = createImage()
    lonImage.load(inps.lon+'.xml')

    width = inImage.width
    length = inImage.length
    if width != latImage.width or width != lonImage.width:
        raise Exception('file widths are different!')
    if length != latImage.length or length != latImage.length:
        raise Exception('file lengths are different!')

    lat_shift = 'lat.tmp'
    lon_shift = 'lon.tmp'
    cmd = "gdal_translate -of VRT -srcwin {} {} {} {} -outsize {} {} -a_nodata 0 {} {}".format(
        inps.leftshift,
        inps.topshift,
        width,
        length,
        width,
        length,
        inps.lat + '.vrt',
        lat_shift + '.vrt')
    runCmd(cmd)
    cmd = "gdal_translate -of VRT -srcwin {} {} {} {} -outsize {} {} -a_nodata 0 {} {}".format(
        inps.leftshift,
        inps.topshift,
        width,
        length,
        width,
        length,
        inps.lon + '.vrt',
        lon_shift + '.vrt')
    runCmd(cmd)

    tree = ET.parse(inps.input + '.vrt')
    root = tree.getroot()
    meta = ET.SubElement(root, 'metadata')
    meta.attrib['domain'] = "GEOLOCATION"
    meta.tail = '\n'
    meta.text = '\n    '
    rdict = {'Y_DATASET' : lat_shift + '.vrt',
            'X_DATASET' :  lon_shift + '.vrt',
            'X_BAND' : "1",
            'Y_BAND' : "1",
            'PIXEL_OFFSET': "0",
            'LINE_OFFSET' : "0",
            'LINE_STEP' : "1",
            'PIXEL_STEP' : "1"}
    for key, val in rdict.items():
        data = ET.SubElement(meta, 'mdi')
        data.text = val
        data.attrib['key'] = key
        data.tail = '\n    '
    data.tail = '\n'
    tree.write(inps.input+'.vrt')

    if os.path.isfile(inps.output):
        os.remove(inps.output)
    wsen = "{} {} {} {}".format(bbox[2], bbox[0], bbox[3], bbox[1])
    rmethod = ['near', 'bilinear', 'cubic', 'cubicspline', 'lanczos', 'average', 'mode', 'max', 'min', 'med', 'Q1', 'Q3']
    cmd = "gdalwarp -of ENVI -geoloc -te {} -tr {} {} -srcnodata 0 -dstnodata 0 -r {} {} {}".format(
        wsen,
        inps.ssize / 3600.0,
        inps.ssize / 3600.0,
        rmethod[inps.rmethod],
        inps.input + '.vrt',
        inps.output)
    runCmd(cmd)

    ds=gdal.Open(inps.output)
    b=ds.GetRasterBand(1)

    #outImage = createImage()
##########################################################
    outImage = inImage
    inImage.scheme = 'BSQ'
##########################################################
    outImage.setFilename(inps.output)
    outImage.setWidth(b.XSize)
    outImage.setLength(b.YSize)
    outImage.coord2.coordDescription = 'Latitude'
    outImage.coord2.coordUnits = 'degree'
    #outImage.coord2.coordStart = bbox[1]
    #outImage.coord2.coordDelta = -sample_size_lat
    outImage.coord2.coordStart = ds.GetGeoTransform()[3]
    outImage.coord2.coordDelta = ds.GetGeoTransform()[5]
    outImage.coord1.coordDescription = 'Longitude'
    outImage.coord1.coordUnits = 'degree'
    #outImage.coord1.coordStart = bbox[2]
    #outImage.coord1.coordDelta = sample_size_lon
    outImage.coord1.coordStart = ds.GetGeoTransform()[0]
    outImage.coord1.coordDelta = ds.GetGeoTransform()[1]
    outImage.renderHdr()


    os.remove(os.path.splitext(inps.output)[0]+'.hdr')
    os.remove(lat_shift + '.vrt')
    os.remove(lon_shift + '.vrt')
    #many things only work in the current directory because of many levels of function calls
    #I do the following to overcome that
##################################################################################################
    if input_cd == 1:
        os.remove(os.path.basename(inps.input))
        os.remove(os.path.basename(inps.input)+'.xml')
    if lat_cd == 1:
        os.remove(os.path.basename(inps.lat))
        os.remove(os.path.basename(inps.lat)+'.xml')
    if lon_cd == 1:
        os.remove(os.path.basename(inps.lon))
        os.remove(os.path.basename(inps.lon)+'.xml')
##################################################################################################

