#!/usr/bin/env python3

#Cunren Liang, JPL/Caltech


import os
import sys
import glob
import shutil
import ntpath
import pickle
import argparse
import datetime
import numpy as np
from scipy.interpolate import griddata

import isce
import isceobj
from isceobj.Image import createImage,createDemImage


def read_bands(filename, length, width, scheme, nbands, datatype):

    if datatype.upper() == 'FLOAT':
        datatype1 = np.float32
    elif datatype.upper() == 'CFLOAT':
        datatype1 = np.complex64
    elif datatype.upper() == 'DOUBLE':
        datatype1 = np.float64
    elif datatype.upper() == 'BYTE':
        datatype1 = np.int8
    elif datatype.upper() == 'SHORT':
        datatype1 = np.int16
    else:
        raise Exception('data type not supported, please contact crl for your data type!')

    bands = []
    if scheme.upper() == 'BIP':
        data = np.fromfile(filename, dtype=datatype1).reshape(length, width*nbands)
        for i in range(nbands):
            bands.append(data[:, i:width*nbands:nbands])
    elif scheme.upper() == 'BIL':
        data = np.fromfile(filename, dtype=datatype1).reshape(length*nbands, width)
        for i in range(nbands):
            bands.append(data[i:length*nbands:nbands, :])
    elif scheme.upper() == 'BSQ':
        data = np.fromfile(filename, dtype=datatype1).reshape(length*nbands, width)
        for i in range(nbands):
            offset = length * i
            bands.append(data[0+offset:length+offset, :])        
    else:
        raise Exception('unknown band scheme!')

    return bands


def write_bands(filename, length, width, scheme, nbands, datatype, bands):

    if datatype.upper() == 'FLOAT':
        datatype1 = np.float32
    elif datatype.upper() == 'CFLOAT':
        datatype1 = np.complex64
    elif datatype.upper() == 'DOUBLE':
        datatype1 = np.float64
    elif datatype.upper() == 'BYTE':
        datatype1 = np.int8
    elif datatype.upper() == 'SHORT':
        datatype1 = np.int16
    else:
        raise Exception('data type not supported, please contact crl for your data type!')

    if scheme.upper() == 'BIP':
        data = np.zeros((length, width*nbands), dtype=datatype1)
        for i in range(nbands):
            data[:, i:width*nbands:nbands] = bands[i]
    elif scheme.upper() == 'BIL':
        data = np.zeros((length*nbands, width), dtype=datatype1)
        for i in range(nbands):
            data[i:length*nbands:nbands, :] = bands[i]
    elif scheme.upper() == 'BSQ':
        data = np.zeros((length*nbands, width), dtype=datatype1)
        for i in range(nbands):
            offset = length * i
            data[0+offset:length+offset, :] = bands[i]
    else:
        raise Exception('unknown band scheme!')    

    #output result
    data.astype(datatype1).tofile(filename)


def cmdLineParse():
    '''
    Command line parser.
    '''
    parser = argparse.ArgumentParser( description='geocode file using lat and lon files',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="")
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
    parser.add_argument('-rmethod', dest='rmethod', type=int, default=1,
            help = 'resampling method. 0: nearest. 1: linear (default). 2: cubic.')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':

    inps = cmdLineParse()

    inImage = createImage()
    inImage.load(inps.input+'.xml')
    inImage.filename = inps.input

    width = inImage.width
    length = inImage.length
    scheme = inImage.scheme
    nbands = inImage.bands
    datatype = inImage.dataType

    #check file width and length
    latImage = createImage()
    latImage.load(inps.lat+'.xml')
    latImage.filename = inps.lat

    lonImage = createImage()
    lonImage.load(inps.lon+'.xml')
    lonImage.filename = inps.lon

    if width != latImage.width or width != lonImage.width:
        raise Exception('file width are different!')
    if length != latImage.length or length != latImage.length:
        raise Exception('file length are different!')

    #convert to degrees
    sample_size = inps.ssize/3600.0
    sample_size_lat = sample_size
    sample_size_lon = sample_size
    print("geocoding sample size: {} [arcsec]".format(inps.ssize))
    print("geocoding sample size: {} [degree]".format(sample_size))

    if inps.rmethod == 0:
        rmethod = 'nearest'
    elif inps.rmethod == 1:
        rmethod = 'linear'
    elif inps.rmethod == 2:
        rmethod = 'cubic'
    else:
        raise Exception('unknow resampling method!')
    print("interpolation method: {}".format(rmethod))

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

    #get input bands saved in a list
    print("reading data to be geocoded")
    bands = read_bands(inps.input, length, width, scheme, nbands, datatype)

    #get latitude
    print("reading latitude")
    lat = read_bands(inps.lat, length, width, latImage.scheme, latImage.bands, latImage.dataType)
    lat = lat[0]
    #get longitude
    print("reading longitude")
    lon = read_bands(inps.lon, length, width, lonImage.scheme, lonImage.bands, lonImage.dataType)
    lon = lon[0]

    latlon = np.zeros((length*width, 2), dtype=np.float64)
    latlon[:, 0] = lat.reshape(length*width)
    latlon[:, 1] = lon.reshape(length*width)
    #                             n       s                         w       e
    grid_lat, grid_lon = np.mgrid[bbox[1]:bbox[0]:-sample_size_lat, bbox[2]:bbox[3]:sample_size_lon]

    print("interpolate input data")
    bands_geo = []
    for i in range(nbands):
        geoband = griddata(latlon, (bands[i]).reshape(length*width), (grid_lat, grid_lon), method=rmethod, fill_value = 0.0)
        bands_geo.append(geoband)

    print("write result")
    (length_geo, width_geo) = geoband.shape
    write_bands(inps.output, length_geo, width_geo, scheme, nbands, datatype, bands_geo)


    outImage = inImage
    outImage.setFilename(inps.output)
    outImage.setWidth(width_geo)
    outImage.setLength(length_geo)
    outImage.coord2.coordDescription = 'Latitude'
    outImage.coord2.coordUnits = 'degree'
    outImage.coord2.coordStart = bbox[1]
    outImage.coord2.coordDelta = -sample_size_lat
    outImage.coord1.coordDescription = 'Longitude'
    outImage.coord1.coordUnits = 'degree'
    outImage.coord1.coordStart = bbox[2]
    outImage.coord1.coordDelta = sample_size_lon
    outImage.renderHdr()



#./geo.py -input filt_diff_150909-160824_8rlks_16alks_msk.unw -output filt_diff_150909-160824_8rlks_16alks_msk.unw.geo -lat 150909-160824_8rlks_16alks.lat -lon 150909-160824_8rlks_16alks.lon -bbox 150909.slc.pck -ssize 1.0


