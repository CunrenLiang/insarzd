#!/usr/bin/env python3


# 1. For ScanSAR-ScanSAR interferometry, this program calculates the burst synchronization.
# 2. For ScanSAR-stripmap interferometry, this program calculates the corresponding burst
#    start lines and times for stripmap, and add them to stripmap frame pickle file. Burst
#    length and burst cycle length from ScanSAR are also added to stripmap frame pickle file.
#    In this case, the burst synchronization is 100%.
# 3. In both cases, the burst synchronization information is reported and written to a txt
#    file.

#Cunren Liang, JPL/Caltech



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


def runTopo(frame, demImage, latName='lat.rdr', lonName='lon.rdr', hgtName='z.rdr', losName='los.rdr', incName = 'inc.rdr', mskName = 'msk.rdr', offsetFromStart=0, numOfLines=100):
    from zerodop.topozero import createTopozero
    from isceobj.Planet.Planet import Planet

    #####Run Topo
    planet = Planet(pname='Earth')
    topo = createTopozero()
    topo.slantRangePixelSpacing = 0.5 * SPEED_OF_LIGHT / frame.rangeSamplingRate
    topo.prf = frame.PRF #!!!should be changed to azimuth time interval for burst mode
    topo.radarWavelength = frame.radarWavelegth
    topo.orbit = frame.getOrbit()
    topo.width = frame.getNumberOfSamples()
    #topo.length = frame.getNumberOfLines()
    topo.length = numOfLines
    topo.wireInputPort(name='dem', object=demImage)
    topo.wireInputPort(name='planet', object=planet)
    topo.numberRangeLooks = 1 #must be set as 1
    topo.numberAzimuthLooks = 1 #must be set as 1 Cunren
    #topo.lookSide = -1
    topo.lookSide = frame.getInstrument().getPlatform().pointingDirection
    topo.sensingStart = frame.getSensingStart() + datetime.timedelta(seconds=offsetFromStart * (1.0/topo.prf))
    topo.rangeFirstSample = frame.startingRange
    topo.demInterpolationMethod='BIQUINTIC'

    topo.latFilename = latName
    topo.lonFilename = lonName
    topo.heightFilename = hgtName
    topo.losFilename = losName
    topo.incFilename = incName
    topo.maskFilename = mskName

    topo.topo()

    return list(topo.snwe)


def runGeo2rdr(frame, latImage, lonImage, demImage, rgoffName='range.off', azoffName='azimuth.off', extentionLines=0):

    from zerodop.geo2rdr import createGeo2rdr
    from isceobj.Planet.Planet import Planet

    #create topo
    planet = Planet(pname='Earth')
    topo = createGeo2rdr()
    topo.configure()

    #set parameters
    topo.slantRangePixelSpacing = 0.5 * SPEED_OF_LIGHT / frame.rangeSamplingRate
    topo.prf = frame.PRF #!!!should be changed to azimuth time interval for burst mode
    topo.radarWavelength = frame.radarWavelegth
    topo.orbit = frame.getOrbit()
    topo.width = frame.getNumberOfSamples()
    topo.length = frame.getNumberOfLines() + extentionLines
    topo.demLength = demImage.length
    topo.demWidth = demImage.width
    topo.wireInputPort(name='planet', object=planet)
    topo.numberRangeLooks = 1 #
    topo.numberAzimuthLooks = 1 # must be set to be 1
    #topo.lookSide = -1
    topo.lookSide = frame.getInstrument().getPlatform().pointingDirection
    topo.setSensingStart(frame.getSensingStart() - datetime.timedelta(seconds=extentionLines * (1.0/topo.prf))          )
    topo.rangeFirstSample = frame.startingRange
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




def cal_offset(mframeFile, sframeFile, DEMfile):

    #master frame
    with open(mframeFile, 'rb') as f:
        mframe = pickle.load(f)

    #slave frame
    with open(sframeFile, 'rb') as f:
        sframe = pickle.load(f)

    #dem
    demImage = isceobj.createDemImage()
    demImage.load(DEMfile + '.xml')
    demImage.setAccessMode('read')

    latFile = 'tmp.lat'
    lonFile = 'tmp.lon'
    hgtFile = 'tmp.hgt'
    losFile = 'tmp.los'
    incFile = 'tmp.inc'
    mskFile = 'tmp.msk'

    rgoffFile = 'range.off'
    azoffFile = 'azimuth.off'

    ##########################################################################
    #offset calculation parameters
    masterOffsetFromStart=int(mframe.getNumberOfLines()/2) #number of lines
    masterNumberOfLines=100 #number of lines
    slaveExtentionLines=100000 #number of lines
    ##########################################################################

    #run topo
    runTopo(mframe, demImage, latName=latFile, lonName=lonFile, hgtName=hgtFile, losName=losFile, incName=incFile, mskName=mskFile, offsetFromStart=masterOffsetFromStart, numOfLines=masterNumberOfLines)


    latImage = isceobj.createImage()
    latImage.load(latFile + '.xml')
    latImage.setAccessMode('read')

    lonImage = isceobj.createImage()
    lonImage.load(lonFile + '.xml')
    lonImage.setAccessMode('read')

    hgtImage = isceobj.createDemImage()
    hgtImage.load(hgtFile + '.xml')
    hgtImage.setAccessMode('read')

    #run geo: extend file length by 2 * extentionLines, to make sure it overlaps with master
    runGeo2rdr(sframe, latImage, lonImage, hgtImage, rgoffName=rgoffFile, azoffName=azoffFile, extentionLines=slaveExtentionLines)

    #get range and azimuth offsets
    width = getWidth(rgoffFile+'.xml')
    length = getLength(rgoffFile+'.xml')
    rgoff = np.fromfile(rgoffFile, dtype=np.float32, count=length*width).reshape(length,width)
    azoff = np.fromfile(azoffFile, dtype=np.float32, count=length*width).reshape(length,width)

    #remove BAD_VALUE = -999999.0 as defined in geo2rdr.f90
    # http://docs.scipy.org/doc/numpy-1.10.0/reference/generated/numpy.mean.html
    # In single precision, mean can be inaccurate
    # np.mean(a, dtype=np.float64)
    rgoffm = np.mean(rgoff[np.nonzero(rgoff!=-999999.0)], dtype=np.float64)
    azoffm = np.mean(azoff[np.nonzero(azoff!=-999999.0)], dtype=np.float64) - masterOffsetFromStart - slaveExtentionLines

    print('\noffsets from geometrical calculation:')
    print('++++++++++++++++++++++++++++++++++++++')
    print('range offset: {}'.format(rgoffm))
    print('azimuth offset: {}\n'.format(azoffm))

    #tidy up
    os.remove(latFile)
    os.remove(latFile + '.xml')
    os.remove(latFile + '.vrt')
    os.remove(lonFile)
    os.remove(lonFile + '.xml')
    os.remove(lonFile + '.vrt')
    os.remove(hgtFile)
    os.remove(hgtFile + '.xml')
    os.remove(hgtFile + '.vrt')
    os.remove(losFile)
    os.remove(losFile + '.xml')
    os.remove(losFile + '.vrt')
    os.remove(incFile)
    os.remove(incFile + '.xml')
    os.remove(incFile + '.vrt')
    os.remove(mskFile)
    os.remove(mskFile + '.xml')
    os.remove(mskFile + '.vrt')

    os.remove(rgoffFile)
    os.remove(rgoffFile + '.xml')
    os.remove(rgoffFile + '.vrt')
    os.remove(azoffFile)
    os.remove(azoffFile + '.xml')
    os.remove(azoffFile + '.vrt')


    return [rgoffm, azoffm]


def cmdLineParse():
    '''
    Command line parser.
    '''
    parser = argparse.ArgumentParser( description='calculate busrt synchronization, if ScanSAR-stripmap interferometry, calculate the corresponding burst start lines and times for stripmap (100% synchronization)')
    parser.add_argument('-m', dest='mframe', type=str, required=True,
            help = 'master frame pickle file')
    parser.add_argument('-s', dest='sframe', type=str, required=True,
            help = 'slave frame pickle file')
    parser.add_argument('-d', dest='dem', type=str, required=True,
            help = 'DEM file')
    parser.add_argument('-c', dest='cmb', type=int, default=0,
            help = 'inteferometry combination. 0: ScanSAR-ScanSAR (default). 1: ScanSAR-stripmap')
    parser.add_argument('-i', dest='ind', type=int, default=0,
            help = 'if inteferometry combination is ScanSAR-stripmap, master is? 0: master is stripmap (default). 1: master is ScanSAR. This option has no effect on ScanSAR-ScanSAR.')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':

    inps = cmdLineParse()

    #master frame
    with open(inps.mframe, 'rb') as f:
        mframe = pickle.load(f)

    #slave frame
    with open(inps.sframe, 'rb') as f:
        sframe = pickle.load(f)


################################################################################################
    if not ((hasattr(mframe, 'nbraw') and hasattr(mframe, 'ncraw')) or (hasattr(sframe, 'nbraw') and hasattr(sframe, 'ncraw'))):
        print('nothing is done by this program.')
        sys.exit(0)
################################################################################################











    offset = cal_offset(inps.mframe, inps.sframe, inps.dem)

    if inps.cmb == 0:
        scburstStartLine = mframe.burstStartLines[0] + offset[1]
        nb = sframe.nbraw
        nc = sframe.ncraw
        #slave burst start times corresponding to master burst start times (100% synchronization)
        scburstStartLines = np.arange(scburstStartLine - 100000*nc, scburstStartLine + 100000*nc, nc)
        dscburstStartLines = -(sframe.burstStartLines[0] - scburstStartLines)
        #find the difference with minimum absolute value
        unsynLines = dscburstStartLines[np.argmin(np.absolute(dscburstStartLines))]
       
        if np.absolute(unsynLines) >= nb:
            synLines = 0
            if unsynLines > 0:
                unsynLines = nb
            else:
                unsynLines = -nb
        else:
            synLines = nb - np.absolute(unsynLines)

        ############################################################################################
        #illustration of the sign of the number of unsynchronized lines (unsynLines)     
        #The convention is the same as ampcor offset, that is,
        #              slaveLineNumber = masterLineNumber + unsynLines
        #
        # |-----------------------|     ------------
        # |                       |        ^
        # |                       |        |
        # |                       |        |   unsynLines < 0
        # |                       |        |
        # |                       |       \ /
        # |                       |    |-----------------------|
        # |                       |    |                       |
        # |                       |    |                       |
        # |-----------------------|    |                       |
        #        Master Burst          |                       |
        #                              |                       |
        #                              |                       |
        #                              |                       |
        #                              |                       |
        #                              |-----------------------|
        #                                     Slave Burst
        #
        #
        ############################################################################################
    else:
        if inps.ind == 0: #stripmap is master
            scburstStartLine = sframe.burstStartLines[0] - offset[1]
            nb = sframe.nbraw
            nc = sframe.ncraw
            prf = sframe.PRF

            #slave burst start times corresponding to master burst start times (100% synchronization)
            burstStartLines=[]
            burstStartTimes=[]
            for i in range(-100000, 100000):
                saz_burstx = scburstStartLine + nc * i
                st_burstx = mframe.getSensingStart() + datetime.timedelta(seconds=saz_burstx * (1.0/prf))
                if saz_burstx >= 0.0 and saz_burstx <= mframe.getNumberOfLines():
                    burstStartLines.append(saz_burstx)
                    burstStartTimes.append(st_burstx)

            mframe.burstStartLines = burstStartLines
            mframe.burstStartTimes = burstStartTimes
            mframe.nbraw = sframe.nbraw
            mframe.ncraw = sframe.ncraw

            #add the burst starting times to stripmap frame pickle
            os.remove(inps.mframe)
            with open(inps.mframe, 'wb') as f:
                pickle.dump(mframe, f)


        else: #scansar is master
            scburstStartLine = mframe.burstStartLines[0] + offset[1]
            nb = mframe.nbraw
            nc = mframe.ncraw
            prf = mframe.PRF

            #slave burst start times corresponding to master burst start times (100% synchronization)
            burstStartLines=[]
            burstStartTimes=[]
            for i in range(-100000, 100000):
                saz_burstx = scburstStartLine + nc * i
                st_burstx = sframe.getSensingStart() + datetime.timedelta(seconds=saz_burstx * (1.0/prf))
                if saz_burstx >= 0.0 and saz_burstx <= sframe.getNumberOfLines():
                    burstStartLines.append(saz_burstx)
                    burstStartTimes.append(st_burstx)

            sframe.burstStartLines = burstStartLines
            sframe.burstStartTimes = burstStartTimes
            sframe.nbraw = mframe.nbraw
            sframe.ncraw = mframe.ncraw

            #add the burst starting times to stripmap frame pickle
            os.remove(inps.sframe)
            with open(inps.sframe, 'wb') as f:
                pickle.dump(sframe, f)


        #in this case burst synchronization is 100%
        unsynLines = 0
        synLines = nb

        #report the burst start lines and times of stripmap data
        print('burst no     start line                 start time')
        print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        for i in range(len(burstStartLines)):
            print('{:4d} {:18.2f}          {}'.format(i, burstStartLines[i], burstStartTimes[i]))


    #report burst synchronization information
    synInfo = '''burst synchronization result
image pair range offset: {}
image pair azimuth offset: {}
number of lines in a burst: {}
number of lines in a burst cycle: {}
number of unsynchronized lines in a burst: {}
burst synchronization: {}%
'''.format(
    offset[0],
    offset[1],
    nb,
    nc,
    unsynLines,
    synLines/nb*100.0
    )

    with open('burst_synchronization.txt', 'w') as f:
        f.write(synInfo)

    print(synInfo)






