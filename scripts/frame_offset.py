#!/usr/bin/env python3

#Cunren Liang, 22-MAY-2015
#JPL/Caltech

############################################################################
#program for calculating offsets of consecutive frames
#
#update history:
#10-apr-2017, crl, support burst-by-burst processing
############################################################################


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
from iscesys.StdOEL.StdOELPy import create_writer
from isceobj.Constants import SPEED_OF_LIGHT

from crlpac import getWidth
from crlpac import getLength
from crlpac import runCmd
from crlpac import writeOffset
from crlpac import meanOffset
from crlpac import cullOffsetbyMean
from crlpac import cullOffset
from crlpac import getOffset


def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="program for calculating offsets of consecutive frames\
                    \n\
                    ")
    parser.add_argument('-s', action='append', dest='slc', default=[],
            help = 'frame SLCs or amplitude images IN ORDER')
    parser.add_argument('-f', action='append', dest='frame', default=[],
            help = 'frame pck file IN THE ORDER OF acquisition time.')
    parser.add_argument('-o', dest='output', type=str, required=True,
            help = 'output txt file containing the offsets')
    parser.add_argument('-c', dest='cflag', type=int, default=0,
            help = 'offsets culling method. 0: fitoff in ROI_pac (default). 1: culling routine in isce. 2: use gemoetrical offsets.')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':

    BIN_DIR='$INSAR_ZERODOP_BIN'

    inps = cmdLineParse()

    slcs = inps.slc
    frames = inps.frame

    if not (len(slcs) == len(frames)):
        raise Exception('the numbers of files are not equal!')
    if len(slcs) < 2:
        raise Exception('at least 2 subswaths should be specified!')

    nframe = len(slcs)

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

    #get offset from matching
    rangeOffset2 = []
    azimuthOffset2 = []
    rangeOffset2.append(0.0)
    azimuthOffset2.append(0.0)
    #do the matching between two consecutive frames:
    for i in range(1, nframe):

        if os.path.splitext(slcs[0])[1] == '.slc':
            mSLC = isceobj.createSlcImage()
        elif os.path.splitext(slcs[0])[1] == '.amp':
            mSLC = isceobj.createImage()
            mSLC.dataType='FLOAT'
        else:
            raise Exception('file type not supported yet.')
        mSLC.setFilename(slcs[i-1])
        mSLC.setWidth(width[i-1])
        mSLC.setLength(length[i-1])
        mSLC.setAccessMode('read')
        mSLC.createImage()

        if os.path.splitext(slcs[0])[1] == '.slc':
            sSLC = isceobj.createSlcImage()
        elif os.path.splitext(slcs[0])[1] == '.amp':
            sSLC = isceobj.createImage()
            sSLC.dataType='FLOAT'
        else:
            raise Exception('file type not supported yet.')
        sSLC.setFilename(slcs[i])
        sSLC.setWidth(width[i])
        sSLC.setLength(length[i])
        sSLC.setAccessMode('read')
        sSLC.createImage()

        ampcor = Ampcor(name='insarapp_slcs_ampcor')
        ampcor.configure()

        #DATA TYPE
        if os.path.splitext(slcs[0])[1] == '.slc':
            ampcor.setImageDataType1('complex')
            ampcor.setImageDataType2('complex')
        elif os.path.splitext(slcs[0])[1] == '.amp':
            ampcor.setImageDataType1('real')
            ampcor.setImageDataType2('real')
        else:
            raise Exception('file type not supported yet.')

        #INPUT/OUTPUT FILES
        ampcor.setMasterSlcImage(mSLC)
        ampcor.setSlaveSlcImage(sSLC)

        #MATCH REGION
        ########################################
        numRange = 60;
        numAzimuth = 30;
        firstSample = 1
        firstLine = 1

        if rangeOffset1[i] < 0:
            firstSample = int(35 - rangeOffset1[i])
        if azimuthOffset1[i] < 0:
            firstLine = int(35 - azimuthOffset1[i])
        ########################################
        ampcor.setFirstSampleAcross(firstSample)
        ampcor.setLastSampleAcross(width[i-1])
        ampcor.setNumberLocationAcross(numRange)
        ampcor.setFirstSampleDown(firstLine)
        ampcor.setLastSampleDown(length[i-1])
        ampcor.setNumberLocationDown(numAzimuth)

        #MATCH PARAMETERS
        ########################################
        fftWidth = 64
        fftLength = 256
        searchMaxWidth = 20
        searchMaxLength = 20

        ########################################
        ampcor.setWindowSizeWidth(fftWidth)
        ampcor.setWindowSizeHeight(fftLength)
        ampcor.setSearchWindowSizeWidth(searchMaxWidth)
        ampcor.setSearchWindowSizeHeight(searchMaxLength)
        ampcor.setAcrossLooks(1)
        ampcor.setDownLooks(1)
        ampcor.setOversamplingFactor(64)
        ampcor.setZoomWindowSize(16)
        ampcor.setAcrossGrossOffset(rangeOffset1[i])
        ampcor.setDownGrossOffset(azimuthOffset1[i])
        #1. The following not set
        #Matching Scale for Sample/Line Directions                       (-)    = 1. 1.
        #should add the following in Ampcor.py?
        #if not set, in this case, Ampcor.py'value is also 1. 1.
        #ampcor.setScaleFactorX(1.)
        #ampcor.setScaleFactorY(1.)


        #MATCH THRESHOLDS AND DEBUG DATA
        #2. The following not set
        #in roi_pac the value is set to 0 1
        #in isce the value is set to 0.001 1000.0
        #SNR and Covariance Thresholds                                   (-)    =  {s1} {s2}
        #should add the following in Ampcor?
        #THIS SHOULD BE THE ONLY THING THAT IS DIFFERENT FROM THAT OF ROI_PAC
        #ampcor.setThresholdSNR(0)
        #ampcor.setThresholdCov(1)
        ampcor.setDebugFlag(False)
        ampcor.setDisplayFlag(False)

        #in summary, only two things not set which are indicated by 'The following not set' above.

        
        #run ampcor
        ampcor.ampcor()
        
        #get offsets
        offsets = ampcor.getOffsetField()
        offsetFile = 'ampcor_{}_{}.off'.format(i-1, i)
        cullOffsetFile = 'cull_{}_{}.off'.format(i-1, i)
        dumpFile = 'fitoff_{}_{}.out'.format(i-1, i)


        #first cull by mean
        threshold = 3.5
        meanOffsets = meanOffset(offsets)
        offsets = cullOffsetbyMean(offsets, meanOffsets[1], threshold)

        #cull offset
        if inps.cflag == 0:
            #cull offsets using ROI_pac's fitoff
            atParam = getOffset(offsets, offsetFile, cullOffsetFile, dumpFile)

            if atParam != None:
                rangeOffset2.append(atParam[4])
                azimuthOffset2.append(atParam[5])
            else:
                print('******************************************************************')
                print('WARNING: There are not enough offsets left, so we are forced to')
                print('         use geometrical offset for frame mosaicking')
                print('******************************************************************')
                rangeOffset2.append(rangeOffset1[i])
                azimuthOffset2.append(azimuthOffset1[i])

            #remove the dumpfile of fitoff
            os.remove(dumpFile)
            os.remove(offsetFile)
            os.remove(cullOffsetFile)
        elif inps.cflag == 1:
            #cull offsets using isce
            #record original offsets
            writeOffset(offsets, offsetFile)
            distances = (10,5,3)
            #numCullOffsetsLimits = (100, 75, 50)
            numCullOffsetsLimits = (50, 40, 30) #crl, apr. 25, 2017
            refinedOffsets = cullOffset(offsets, distances, numCullOffsetsLimits)

            if refinedOffsets != None:
                #record culled offsets
                writeOffset(refinedOffsets, cullOffsetFile)
                #use mean offsets
                rangeOffset2.append(meanOffset(refinedOffsets)[0])
                azimuthOffset2.append(meanOffset(refinedOffsets)[1])
                os.remove(cullOffsetFile)
            else:
                print('******************************************************************')
                print('WARNING: There are not enough offsets left, so we are forced to')
                print('         use geometrical offset for frame mosaicking')
                print('******************************************************************')
                rangeOffset2.append(rangeOffset1[i])
                azimuthOffset2.append(azimuthOffset1[i])
            os.remove(offsetFile)

        else:

            print('******************************************************************')
            print('WARNING: You are using geometrical offset for frame mosaicking')
            print('******************************************************************')


            #just use geometrical offsets
            rangeOffset2.append(rangeOffset1[i])
            azimuthOffset2.append(azimuthOffset1[i])


    for i in range(1, nframe):
        if abs(azimuthOffset2[i] - azimuthOffset1[i]) > 4.0:
            print('******************************************************************')
            print("WARNING: There may be an error with azimuth offset {} probably caused by offset culling, using geometrical azimuth offset instead!".format(i))
            print('******************************************************************')
            azimuthOffset2[i] = azimuthOffset1[i]


    offsetComp = "\n\ncomparision of offsets:\n\n"
    offsetComp += "offset type       i     geometrical           match      difference\n"
    offsetComp += "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
    for i, (offset1, offset2) in enumerate(zip(rangeOffset1, rangeOffset2)):
        offsetComp += "range offset     {:2d}   {:13.3f}   {:13.3f}   {:13.3f}\n".format(i, offset1, offset2, offset1 - offset2)
    for i, (offset1, offset2) in enumerate(zip(azimuthOffset1, azimuthOffset2)):
        offsetComp += "azimuth offset   {:2d}   {:13.3f}   {:13.3f}   {:13.3f}\n".format(i, offset1, offset2, offset1 - offset2)

    with open(inps.output, 'w') as f:
        f.write(offsetComp)

    print("{}".format(offsetComp))


