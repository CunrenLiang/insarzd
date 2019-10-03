#!/usr/bin/env python3

import os
import sys
import pickle
import datetime
import argparse
import numpy as np

import isce
import isceobj
from isceobj.Constants import SPEED_OF_LIGHT


def runGeo(frame, demImage, inImage, bbox, rlks, alks, method, outputName, topShift, leftShift):
    from zerodop.geozero import createGeozero
    from isceobj.Planet.Planet import Planet

    #####Run Topo
    planet = Planet(pname='Earth')
    topo = createGeozero()
    topo.configure()

    topo.slantRangePixelSpacing = 0.5 * SPEED_OF_LIGHT / frame.rangeSamplingRate
    topo.prf = frame.PRF #should be changed to reciprocal of azimuth time interval for burst mode!!!
    topo.radarWavelength = frame.radarWavelegth
    topo.orbit = frame.getOrbit()
    topo.width = frame.getNumberOfSamples()
    topo.length = frame.getNumberOfLines()
    topo.wireInputPort(name='dem', object=demImage)
    topo.wireInputPort(name='planet', object=planet)
    topo.wireInputPort(name='tobegeocoded', object=inImage)
    topo.numberRangeLooks = rlks
    topo.numberAzimuthLooks = alks
    #topo.lookSide = -1
    topo.lookSide = frame.getInstrument().getPlatform().pointingDirection
    topo.setSensingStart(frame.getSensingStart())
    topo.rangeFirstSample = frame.startingRange
    topo.method=method
    topo.demCropFilename = 'crop.dem'
    topo.geoFilename = outputName
    topo.dopplerCentroidCoeffs = [0.]
    topo.snwe = bbox

    ##consider the shifts
    topo.rangeFirstSample += leftShift * topo.slantRangePixelSpacing
    topo.setSensingStart(frame.getSensingStart()+datetime.timedelta(seconds=(topShift/topo.prf)))



    topo.geocode()

    print('South: ', topo.minimumGeoLatitude)
    print('North: ', topo.maximumGeoLatitude)
    print('West:  ', topo.minimumGeoLongitude)
    print('East:  ', topo.maximumGeoLongitude)
    
    return

def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser( description='geocode InSAR products')
    parser.add_argument('-m', '--mframe', dest='mframe', type=str, required=True,
            help = 'master frame (picke file)')
    parser.add_argument('-d', '--dem', dest='dem', type=str, required=True,
            help = 'input DEM to use')
    parser.add_argument('-i', '--input', dest='input', type=str, required=True,
            help='input file to be geocoded')
    parser.add_argument('-o', '--output', dest='output', type=str, required=True,
            help='output geocoded file')
    parser.add_argument('-b', '--bbox', dest='bbox', type=str, default='useMasterDefaultBbox',
            help='Bounding box (SNWE). Default: read from master frame')
    parser.add_argument('-r','--rlks', dest='rlks', type=int, default=1,
            help = 'Number of range looks')
    parser.add_argument('-a','--alks', dest='alks', type=int, default=1,
            help = 'Number of azimuth looks')
    parser.add_argument('-ts','--topshift', dest='topshift', type=int, default=0,
            help = 'top shift in original SLC lines')
    parser.add_argument('-ls','--leftshift', dest='leftshift', type=int, default=0,
            help = 'left shift in original SLC samples')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':

    inps = cmdLineParse()

    #get info from master frame
    with open(inps.mframe, 'rb') as fid:
        frame = pickle.load(fid)

    #check bbox
    if inps.bbox != 'useMasterDefaultBbox':
        inps.bbox = [float(val) for val in inps.bbox.split('/')]
        if len(inps.bbox) != 4:
            raise Exception('Bbox should contain 4 floating point values')
    else:
        inps.bbox = frame.snwe

    #setup dem
    demImage = isceobj.createDemImage()
    demImage.load(inps.dem + '.xml')
    demImage.setAccessMode('read')

    #setup input file
    inImage = isceobj.createImage()
    inImage.load(inps.input + '.xml')
    inImage.setAccessMode('read')

    #choose method here!!!
    if inps.input.endswith('.int'):
        method = 'sinc'
    else:
        method = 'bilinear'

    runGeo(frame, demImage, inImage, inps.bbox, inps.rlks, inps.alks, method, inps.output, inps.topshift, inps.leftshift)


#./geo.py -m 20130927.slc.pck -d demLat_N25_N29_Lon_E060_E067.dem.wgs84 -i filt_diff_20130927-20141211_16rlks_16alks.unw -o filt_diff_20130927-20141211_16rlks_16alks.unw.geo -b 26.12/26.62/64.87/65.28 -r 16 -a 16
#./geo.py -m 20130927.slc.pck -d demLat_N25_N29_Lon_E060_E067.dem.wgs84 -i diff_20130927-20141211_16rlks_16alks.int -o diff_20130927-20141211_16rlks_16alks.int.geo -b 26.12/26.62/64.87/65.28 -r 16 -a 16

