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

def runGeo2rdr(frame, latImage, lonImage, demImage, rgoffName='range.off', azoffName='azimuth.off', rlks=1, alks=1):

    from zerodop.geo2rdr import createGeo2rdr
    from isceobj.Planet.Planet import Planet

    #create topo
    planet = Planet(pname='Earth')
    topo = createGeo2rdr()
    topo.configure()

    #set parameters
    topo.slantRangePixelSpacing = rlks * 0.5 * SPEED_OF_LIGHT / frame.rangeSamplingRate
    topo.prf = frame.PRF / alks #!!!should be changed to azimuth time interval for burst mode
    topo.radarWavelength = frame.radarWavelegth
    topo.orbit = frame.getOrbit()
    topo.width = int(frame.getNumberOfSamples()/rlks)
    topo.length = int(frame.getNumberOfLines()/alks)
    topo.demLength = demImage.length
    topo.demWidth = demImage.width
    topo.wireInputPort(name='planet', object=planet)
    topo.numberRangeLooks = 1 #
    topo.numberAzimuthLooks = 1 # must be set to be 1
    #topo.lookSide = -1
    topo.lookSide = frame.getInstrument().getPlatform().pointingDirection
    topo.setSensingStart(frame.getSensingStart() + datetime.timedelta(seconds=(alks-1.0)/2.0 * (1.0/frame.PRF)))
    topo.rangeFirstSample = frame.startingRange + (rlks - 1.0)/2.0 * (0.5 * SPEED_OF_LIGHT / frame.rangeSamplingRate)
    topo.dopplerCentroidCoeffs = [0.] # we are using zero doppler geometry

    #set files
    topo.demImage = demImage
    topo.latImage = latImage
    topo.lonImage = lonImage
    topo.rangeOffsetImageName = rgoffName
    topo.azimuthOffsetImageName = azoffName

    #run it
    topo.geo2rdr()

    return


def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser( description='offset estimation for the InSAR pair')
    parser.add_argument('-s', '--sframe', dest='sframe', type=str, required=True,
            help = 'slave frame (picke file)')
    parser.add_argument('-a','--lat', dest='lat', type=str, required=True,
            help = 'lattitude file')
    parser.add_argument('-o','--lon', dest='lon', type=str, required=True,
            help = 'longtitude file')
    parser.add_argument('-z','--hgt', dest='hgt', type=str, required=True,
            help = 'height file (DEM in radar coordinates)')
    parser.add_argument('-r','--rgoff', dest='rgoff', type=str, default='range.off',
            help = '(output) range offset file')
    parser.add_argument('-i','--azoff', dest='azoff', type=str, default='azimuth.off',
            help = '(output) azimuth offset file')
    parser.add_argument('-rlks', dest='rlks', type=int, default=1,
            help = 'number of range looks')
    parser.add_argument('-alks', dest='alks', type=int, default=1,
            help = 'number of azimuth looks')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':

    inps = cmdLineParse()

    with open(inps.sframe, 'rb') as fid:
        slaveFrame = pickle.load(fid)

    latImage = isceobj.createImage()
    latImage.load(inps.lat + '.xml')
    latImage.setAccessMode('read')

    lonImage = isceobj.createImage()
    lonImage.load(inps.lon + '.xml')
    lonImage.setAccessMode('read')

    demImage = isceobj.createDemImage()
    demImage.load(inps.hgt + '.xml')
    demImage.setAccessMode('read')
    
    #run it
    runGeo2rdr(slaveFrame, latImage, lonImage, demImage, rgoffName=inps.rgoff, azoffName=inps.azoff, rlks=inps.rlks, alks=inps.alks)


#./geo2rdr.py -s 20141211.slc.pck -a 20130927-20141211.lat -o 20130927-20141211.lon -z 20130927-20141211.hgt -r 20130927-20141211_rg.off -i 20130927-20141211_az.off




