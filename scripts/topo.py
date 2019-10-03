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

from crlpac import getWidth
from crlpac import getLength
from crlpac import runCmd
from crlpac import create_xml


def runTopo(frame, demImage, latName='lat.rdr', lonName='lon.rdr', hgtName='z.rdr', losName='los.rdr', incName = 'inc.rdr', mskName = 'msk.rdr', rlks=1, alks=1):
    from zerodop.topozero import createTopozero
    from isceobj.Planet.Planet import Planet

    #####Run Topo
    planet = Planet(pname='Earth')
    topo = createTopozero()
    topo.slantRangePixelSpacing = rlks * 0.5 * SPEED_OF_LIGHT / frame.rangeSamplingRate
    topo.prf = frame.PRF / alks #!!!should be changed to azimuth time interval for burst mode
    topo.radarWavelength = frame.radarWavelegth
    topo.orbit = frame.getOrbit()
    topo.width = int(frame.getNumberOfSamples()/rlks)
    topo.length = int(frame.getNumberOfLines()/alks)
    topo.wireInputPort(name='dem', object=demImage)
    topo.wireInputPort(name='planet', object=planet)
    topo.numberRangeLooks = 1 #must be set as 1
    topo.numberAzimuthLooks = 1 #must be set as 1 Cunren
    #topo.lookSide = -1
    topo.lookSide = frame.getInstrument().getPlatform().pointingDirection
    topo.sensingStart = frame.getSensingStart() + datetime.timedelta(seconds=(alks-1.0)/2.0 * (1.0/frame.PRF))
    topo.rangeFirstSample = frame.startingRange + (rlks - 1.0)/2.0 * (0.5 * SPEED_OF_LIGHT / frame.rangeSamplingRate)
    topo.demInterpolationMethod='BIQUINTIC'

    topo.latFilename = latName
    topo.lonFilename = lonName
    topo.heightFilename = hgtName
    topo.losFilename = losName
    topo.incFilename = incName
    topo.maskFilename = mskName

    topo.topo()

    return list(topo.snwe)

def runSimAmp(hgtName='z.rdr', simName='simamp.rdr'):

    #get width from hgt file
    hgtImage = isceobj.createImage()
    hgtImage.load(hgtName + '.xml')
    hgtImage.setAccessMode('read')
    hgtImage.createImage()

    hgtWidth = hgtImage.getWidth()
    hgtLength = hgtImage.getLength()
    hgtImage.finalizeImage()

    #run simamp
    cmd = "$INSAR_ZERODOP_BIN/simamp {} {} {} 3.0 100.0".format(hgtName, simName, hgtWidth)
    #print("{}".format(cmd))
    runCmd(cmd)

    #create xml file of simamp.rdr
    create_xml(simName, hgtWidth, hgtLength, 'float')


def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser( description='calculate lattitude, longtitude, height, los, simulated radar image, incidence angle, and layover mask file in radar geometry')
    parser.add_argument('-m', '--mframe', dest='mframe', type=str, required=True,
            help = 'master frame (picke file)')
    parser.add_argument('-d', '--dem', dest='dem', type=str, required=True,
            help = 'input DEM to use')
    parser.add_argument('-a','--lat', dest='lat', type=str, default='lat.rdr',
            help = '(output) lattitude file')
    parser.add_argument('-o','--lon', dest='lon', type=str, default='lon.rdr',
            help = '(output) longtitude file')
    parser.add_argument('-z','--hgt', dest='hgt', type=str, default='z.rdr',
            help = '(output) height file (DEM in radar coordinates)')
    parser.add_argument('-l','--los', dest='los', type=str, default='los.rdr',
            help = '(output) line of sight file')
    parser.add_argument('-s','--sim', dest='sim', type=str, default='simamp.rdr',
            help = '(output) simulated radar images')
    parser.add_argument('-i','--inc', dest='inc', type=str, default='inc.rdr',
            help = '(output) incidence angle file')
    parser.add_argument('-k','--msk', dest='msk', type=str, default='msk.rdr',
            help = '(output) layover mask file')    
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

    #get radar parameters from master frame
    with open(inps.mframe, 'rb') as fid:
        masterFrame = pickle.load(fid)

    #Setup dem
    demImage = isceobj.createDemImage()
    demImage.load(inps.dem + '.xml')
    demImage.setAccessMode('read')

    #1. produce lat, lon, hgt, and los
    snwe = runTopo(masterFrame, demImage, latName=inps.lat, lonName=inps.lon, hgtName=inps.hgt, losName=inps.los, incName = inps.inc, mskName = inps.msk, rlks=inps.rlks, alks=inps.alks)

    #add snwe to masterFrame, which is used in geo.py
    masterFrame.snwe = snwe
    os.remove(inps.mframe)
    with open(inps.mframe, 'wb') as f:
        pickle.dump(masterFrame, f)

    #2. produce simulated radar image
    runSimAmp(hgtName=inps.hgt, simName=inps.sim)

#./topo.py -m 20130927.slc.pck -d demLat_N25_N29_Lon_E060_E067.dem.wgs84 -a 20130927-20141211.lat -o 20130927-20141211.lon -z 20130927-20141211.hgt -l 20130927-20141211.los -s 20130927-20141211.sim -i 20130927-20141211.inc -k 20130927-20141211.msk



