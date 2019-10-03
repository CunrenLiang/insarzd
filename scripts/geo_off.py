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

from crlpac import Dummy

def runGeo(frame, demImage, inImage, bbox, rlks, alks, method, outputName, aff, topShift, leftShift):
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

    ##############
    topo.ATCa = aff.a
    topo.ATCb = aff.b
    topo.ATCc = aff.c
    topo.ATCd = aff.d
    topo.ATCe = aff.e
    topo.ATCf = aff.f
    topo.ATCrlks = aff.rlks
    topo.ATCalks = aff.alks
    ##############

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
    parser.add_argument('-t', '--aff', dest='aff', type=str, default=None,
            help='affine transformation file')
    parser.add_argument('-x','--xrlks', dest='xrlks', type=int, default=1,
            help = 'Number of range looks used in the affine transformation')
    parser.add_argument('-y','--yalks', dest='yalks', type=int, default=1,
            help = 'Number of azimuth looks used in the affine transformation')
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

##################################################
    with open(inps.aff) as f:
        lines = f.readlines()
        f.close

    i = 0
    for linex in lines:
        if 'Affine Matrix ' in linex:
            a = float(lines[i + 2].split()[0])
            b = float(lines[i + 2].split()[1])
            c = float(lines[i + 3].split()[0])
            d = float(lines[i + 3].split()[1])
            e  = float(lines[i + 7].split()[0])
            f  = float(lines[i + 7].split()[1])
            break
        i += 1    

    #inverse affine transformation
    affi = Dummy()
    p = a * d - b * c
    affi.a = d / p
    affi.b = -b / p
    affi.c = -c / p
    affi.d = a / p
    affi.e = (b * f - d * e) / p
    affi.f = (c * e - a * f) / p

    affi.rlks = inps.xrlks
    affi.alks = inps.yalks

    print(a, b)
    print(c, d)
    print(e, f)
    print('***********')
    print(affi.a, affi.b)
    print(affi.c, affi.d)
    print(affi.e, affi.f)
    print(affi.rlks, affi.alks)
##################################################




    runGeo(frame, demImage, inImage, inps.bbox, inps.rlks, inps.alks, method, inps.output, affi, inps.topshift, inps.leftshift)


#./geo_off.py -m 141018.slc.pck -d 141018-141123.dem.wgs84 -i filt_diff_141018-141123_16rlks_16alks.unw -o filt_diff_141018-141123_16rlks_16alks.unw.geo -b useMasterDefaultBbox -r 16 -a 16 -t ampsim_16rlks_16alks.aff -x 16 -y 16

#example after adding top and left shifts. 28-AUG-2015
#./geo_off.py -m ../150412.slc.pck -d ./demLat_N22_N33_Lon_E078_E092.dem.wgs84 -i offsets.bil -o offsets.bil.geo -b useMasterDefaultBbox -r 64 -a 64 -t ../ampsim.aff -x 16 -y 16 -ts 512 -ls 512
#top and left shifts can be used when geocoding pixel offset results because edges are cut off when doing pixel offset.


