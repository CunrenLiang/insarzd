#!/usr/bin/env python3

#Cunren Liang, 22-MAY-2015
#JPL/Caltech

#./mosaicframe.py -f ../150221-150502/f500/150221.slc.pck -f ../150221-150502/f510/150221.slc.pck -f ../150221-150502/f520/150221.slc.pck -o 150221.slc.pck
#./mosaicframe.py -f ../150221-150502/f500/150502.slc.pck -f ../150221-150502/f510/150502.slc.pck -f ../150221-150502/f520/150502.slc.pck -o 150502.slc.pck

import os
import sys
import argparse
import datetime
import pickle
import ntpath
import shutil
from xml.etree.ElementTree import ElementTree

import isce
import isceobj
import mroipac
from mroipac.ampcor.Ampcor import Ampcor
from isceobj.Orbit.Orbit import Orbit
from isceobj.Constants import SPEED_OF_LIGHT


def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser( description='mosaic consecutive frames')
    parser.add_argument('-f', action='append', dest='frames', default=[],
            help = 'frames IN THE ORDER OF acquisition time. This should be frame pickle files.')
    parser.add_argument('-o', '--output', dest='output', type=str, required=True,
            help = '(output) pickle file of mosaicked frame. A parameter file is also generated at the same time.')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':

    inps = cmdLineParse()
    if len(inps.frames) < 2:
        raise Exception('At least 2 frames should be specified!')

    print(inps.frames)

    BIN_DIR='$INSAR_ZERODOP_BIN'

    nframe = len(inps.frames)
    frames = inps.frames

    ############################
    pcks = []
    
    dopCoeff = []
    
    width = []
    length = []
    startingRange = []
    farRange = []
    sensingStart = []
    sensingStop = []
    orbit = []

    prf = []
    rangePixelSize = []
    ############################

    for frame in frames:
        with open(frame, 'rb') as fid:
            pcks.append(pickle.load(fid))

    for pck in pcks:
        dopCoeff.append(pck.dopCoeff)
        
        width.append(pck.getNumberOfSamples())
        length.append(pck.getNumberOfLines())
        startingRange.append(pck.startingRange)
        farRange.append(pck.getFarRange())
        sensingStart.append(pck.getSensingStart())
        sensingStop.append(pck.getSensingStop())
        orbit.append(pck.getOrbit())
        
        prf.append(pck.PRF)
        rangePixelSize.append(0.5 * SPEED_OF_LIGHT / pck.rangeSamplingRate)
    
    ################################################
    #should add range/azimuth pixel size  check here
    #currently only handle the case that all frames
    #have the same range/azimuth pixle size.
    ################################################

    #get offset from gemoetrical parameters
    rangeOffset1 = []
    azimuthOffset1 = []
    rangeOffset1.append(0.0)
    azimuthOffset1.append(0.0)
    for i in range(1, nframe):
        rangeOffset1.append(-(startingRange[i] - startingRange[i-1]) / rangePixelSize[i]) 
        azTimeDiff = (sensingStart[i] - sensingStart[i-1]).total_seconds()
        azimuthOffset1.append(-azTimeDiff * prf[i])
    
    #change processing dir
    #os.mkdir('proc')
        

    rangeOffset = []
    azimuthOffset = []
    rangeOffset.append(0.0)
    azimuthOffset.append(0.0)
    #get offset relative to the first frame
    for i in range(1, nframe):
        rangeOffset.append(0.0)
        azimuthOffset.append(0.0)
        for j in range(1, i+1):
            rangeOffset[i] += rangeOffset1[j]
            azimuthOffset[i] += azimuthOffset1[j]


    xs = []
    xe = []
    ys = []
    ye = []
    #determine output width and length
    #actually no need to calculate in azimuth direction
    for i in range(nframe):
        if i == 0:
            xs.append(0)
            xe.append(width[i] - 1)
            ys.append(0)
            ye.append(length[i] - 1)
        else:
            xs.append(0 - int(rangeOffset[i]))
            xe.append(width[i] - 1 - int(rangeOffset[i]))
            ys.append(0 - int(azimuthOffset[i]))
            ye.append(length[i] - 1 - int(azimuthOffset[i]))

    (xmin, xminIndex) = min((v,i) for i,v in enumerate(xs))
    (xmax, xmaxIndex) = max((v,i) for i,v in enumerate(xe))
    (ymin, yminIndex) = min((v,i) for i,v in enumerate(ys))
    (ymax, ymaxIndex) = max((v,i) for i,v in enumerate(ye))

    outWidth = xmax - xmin + 1
    outLength = ymax - ymin + 1

    #update range offset here!!!

    #prepare offset for mosaicing
    rangeOffset3 = []
    azimuthOffset3 = []
    for i in range(nframe):
        rangeOffset3.append(int(rangeOffset[i]) - int(rangeOffset[xminIndex]))
        if i != 0:
            azimuthOffset3.append(int(azimuthOffset[i]) - int(azimuthOffset[i-1]))
        else:
            azimuthOffset3.append(0)


    #update pickle file
    outDopCoeff = [0, 0, 0, 0]
    for dopCoeffx in dopCoeff:
        for i in range(4):
            outDopCoeff[i] += dopCoeffx[i] / nframe

    #outStartingRange = startingRange[0] - rangeOffset3[0] * rangePixelSize[0]
    #fixed the uppler line's bug, CL, 03-JUN-2015
    outStartingRange = startingRange[0] + rangeOffset3[0] * rangePixelSize[0]
    outFarRange = outStartingRange + (outWidth - 1) * rangePixelSize[0]
    outSensingStart = sensingStart[0]
    outSensingStop = outSensingStart + datetime.timedelta(seconds=(outLength - 1) / prf[0])


    ##
    outpck = pcks[0]
    outpck.dopCoeff = outDopCoeff
    outpck.setNumberOfSamples(outWidth)
    outpck.setNumberOfLines(outLength)
    outpck.startingRange = outStartingRange
    outpck.setFarRange(outFarRange)
    outpck.setSensingStart(outSensingStart)
    outpck.setSensingStop(outSensingStop)
    ##still needs to determin which orbit to use:
    #outpck.setOrbit(orbit[nframe - 1])


    #########################################
    #trim orbit: remove the extended orbits when reading data, except the start of
    #first frame and the end of last frame. CL, 14-AUG-2015
    orbitTrim = []
    for i in range(nframe):
        sensingStart1 = sensingStart[i] - datetime.timedelta(seconds=1.0 / prf[i])
        sensingStop1  = sensingStop[i] + datetime.timedelta(seconds=1.0 / prf[i])

        if i == 0:
            print('trimming orbit {}'.format(i))
            sensingStart1 = pcks[i].getOrbit().minTime - datetime.timedelta(seconds=1.0 / prf[i])
        elif i == nframe - 1:
            print('trimming orbit {}'.format(i))
            sensingStop1  = pcks[i].getOrbit().maxTime + datetime.timedelta(seconds=1.0 / prf[i])
        else:
            print('trimming orbit {}'.format(i))

        orbitTrim.append(orbit[i].trimOrbit(sensingStart1, sensingStop1))

    ###########################################
    #mosaic orbits
    #update all the orbit. CL, 14-AUG-2015
    orb = orbitTrim[0]
    for orbitTrimx in orbitTrim[1:]:
        for sv in orbitTrimx:
            if sv.time > orb.maxTime:
                orb.addStateVector(sv)
    outpck.setOrbit(orb)
    ###########################################


    ###############################################################################################
    #with open('./alos2/new/150221-150502/f560/150221.slc.pck', 'rb') as fid:
    #    tmppck = pickle.load(fid)
    #    outpck.setOrbit(tmppck.getOrbit())
    #with open('./alos2/new/150221-150502/f560/150502.slc.pck', 'rb') as fid:
    #    tmppck = pickle.load(fid)
    #    outpck.setOrbit(tmppck.getOrbit())
    ###############################################################################################


    #sensing mid
    outpck.setSensingMid(outpck.getSensingStart() + datetime.timedelta(seconds=(outpck.getSensingStop()-outpck.getSensingStart()).total_seconds()/2.0))
    with open(inps.output, 'wb') as f:
        pickle.dump(outpck, f)



