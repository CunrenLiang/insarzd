#!/usr/bin/env python3

#Cunren Liang, 26-MAY-2015
#JPL/Caltech

import os
import re
import sys
import pickle
import argparse
import numpy as np
from xml.etree.ElementTree import ElementTree

import isce
import stdproc
import isceobj
import mroipac
from mroipac.ampcor.Ampcor import Ampcor
from isceobj.Location.Offset import OffsetField,Offset
from isceobj.Constants import SPEED_OF_LIGHT
from iscesys.StdOEL.StdOELPy import create_writer
from zerodop.baseline.Baseline import Baseline

from crlpac import getWidth
from crlpac import getLength
from crlpac import runCmd
from crlpac import writeOffset
from crlpac import meanOffset
from crlpac import cullOffsetbyMean
from crlpac import cullOffset
from crlpac import getOffset


def readOffset(filename):

    with open(filename, 'r') as f:
        lines = f.readlines()
    #                                          0      1       2       3      4         5             6          7
    #retstr = "%s %s %s %s %s %s %s %s" % (self.x,self.dx,self.y,self.dy,self.snr, self.sigmax, self.sigmay, self.sigmaxy)

    offsets = OffsetField()
    for linex in lines:
        #linexl = re.split('\s+', linex)
        #detect blank lines with only spaces and tabs
        if linex.strip() == '':
            continue

        linexl = linex.split()
        offset = Offset()
        #offset.setCoordinate(int(linexl[0]),int(linexl[2]))
        offset.setCoordinate(float(linexl[0]),float(linexl[2]))
        offset.setOffset(float(linexl[1]),float(linexl[3]))
        offset.setSignalToNoise(float(linexl[4]))
        offset.setCovariance(float(linexl[5]),float(linexl[6]),float(linexl[7]))
        offsets.addOffset(offset)

    return offsets


def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser( description='match, cull offsets and form interferogram')
    parser.add_argument('-m', '--master', dest='master', type=str, required=True,
            help = 'master image for matching')
    parser.add_argument('-s', '--slave', dest='slave', type=str, required=True,
            help = 'slave image for matching')

    parser.add_argument('-mframe', '--masterframe', dest='masterframe', type=str, default=None,
            help = 'master frame. default: master.pck')
    parser.add_argument('-sframe', '--slaveframe', dest='slaveframe', type=str, default=None,
            help = 'slave frame. default: slave.pck')
    parser.add_argument('-m2', '--master2', dest='master2', type=str, default=None,
            help = 'master image for forming the interferogram. default: master')
    parser.add_argument('-s2', '--slave2', dest='slave2', type=str, default=None,
            help = 'slave image for forming the interferogram. default: slave')

    parser.add_argument('-offset', '--offsetfile', dest='offsetfile', type=str, default=None,
            help = 'offset file to use. if set, matching won\'t be done anymore. In this case, other parameters to be set are only -m, -s, -int, -amp, -sframe. default: None')

    parser.add_argument('-int', '--interferogram', dest='interferogram', type=str, required=True,
            help = '(output) interferogram')
    parser.add_argument('-amp', '--amplitude', dest='amplitude', type=str, required=True,
            help = '(output) amplitude image')
    parser.add_argument('-nr', '--nrgoff', dest='nrgoff', type=int, default=40,
            help = 'number of range offsets')
    parser.add_argument('-na', '--nazoff', dest='nazoff', type=int, default=40,
            help = 'number of azimuth offsets')
    parser.add_argument('-rfw', '--rgfftw', dest='rgfftw', type=int, default=32,
            help = 'range fft window size')
    parser.add_argument('-afw', '--azfftw', dest='azfftw', type=int, default=32,
            help = 'azimuth fft window size')
    parser.add_argument('-rsm', '--rgsmax', dest='rgsmax', type=int, default=32,
            help = 'range search samples')
    parser.add_argument('-asm', '--azsmax', dest='azsmax', type=int, default=32,
            help = 'azimuth search lines')
    parser.add_argument('-rlks', dest='rlks', type=int, default=1,
            help = 'number of range looks of the interferograms and amplitudes. default: 1')
    parser.add_argument('-alks', dest='alks', type=int, default=1,
            help = 'number of azimuth looks of the interferograms and amplitudes. default: 1')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':
    
    inps = cmdLineParse()

    masterWidth = getWidth(inps.master + '.xml')
    masterLength = getLength(inps.master + '.xml')

    slaveWidth = getWidth(inps.slave + '.xml')
    slaveLength = getLength(inps.slave + '.xml')


    if inps.offsetfile != None:
        refinedOffsets = readOffset(inps.offsetfile)

        if inps.slaveframe == None:
            slaveFrameName = inps.slave + '.pck'
        else:
            slaveFrameName = inps.slaveframe

        with open(slaveFrameName, 'rb') as f:
            slaveFrame = pickle.load(f)

    else:
#####################################################################
#   STEP 0. calculate mean offsets
#####################################################################
        if inps.masterframe == None:
            masterFrameName = inps.master + '.pck'
        else:
            masterFrameName = inps.masterframe
        if inps.slaveframe == None:
            slaveFrameName = inps.slave + '.pck'
        else:
            slaveFrameName = inps.slaveframe

        with open(masterFrameName, 'rb') as f:
            masterFrame = pickle.load(f)
        with open(slaveFrameName, 'rb') as f:
            slaveFrame = pickle.load(f)

        bObj = Baseline()
        bObj.configure()
        bObj.baselineLocation = 'top'
        bObj.wireInputPort(name='masterFrame', object=masterFrame)
        bObj.wireInputPort(name='slaveFrame', object=slaveFrame)
        azoff, rgoff = bObj.baseline()
        coarseAcross = int(np.mean(rgoff))
        coarseDown = int(np.mean(azoff))



#####################################################################
#   STEP 1. do the matching
#####################################################################
        mSLC = isceobj.createSlcImage()
        mSLC.setFilename(inps.master)
        mSLC.setWidth(masterWidth)
        mSLC.setLength(masterLength)
        mSLC.setAccessMode('read')
        mSLC.createImage()

        sSLC = isceobj.createSlcImage()
        sSLC.setFilename(inps.slave)
        sSLC.setWidth(slaveWidth)
        sSLC.setLength(slaveLength)
        sSLC.setAccessMode('read')
        sSLC.createImage()

        ampcor = Ampcor(name='insarapp_slcs_ampcor')
        ampcor.configure()

        #DATA TYPE
        ampcor.setImageDataType1('complex')
        ampcor.setImageDataType2('complex')

        #INPUT/OUTPUT FILES
        ampcor.setMasterSlcImage(mSLC)
        ampcor.setSlaveSlcImage(sSLC)

        #MATCH REGION
        ########################################
        numRange = inps.nrgoff
        numAzimuth = inps.nazoff
        meanOffsetSamples = coarseAcross
        meanOffsetLines = coarseDown
        #it seems that we cannot use 0, don't know why
        if meanOffsetSamples == 0:
            meanOffsetSamples = 1
        if meanOffsetLines  == 0:
            meanOffsetLines = 1


        firstSample = 1
        firstLine = 1
        if meanOffsetSamples < 0:
            firstSample = int(35 - meanOffsetSamples)
        if meanOffsetLines < 0:
            firstLine = int(35 - meanOffsetLines)
        ########################################
        ampcor.setFirstSampleAcross(firstSample)
        ampcor.setLastSampleAcross(masterWidth)
        ampcor.setNumberLocationAcross(numRange)
        ampcor.setFirstSampleDown(firstLine)
        ampcor.setLastSampleDown(masterLength)
        ampcor.setNumberLocationDown(numAzimuth)

        #MATCH PARAMETERS
        ########################################
        fftWidth = inps.rgfftw
        fftLength = inps.azfftw
        searchMaxWidth = inps.rgsmax
        searchMaxLength = inps.azsmax

        ########################################
        ampcor.setWindowSizeWidth(fftWidth)
        ampcor.setWindowSizeHeight(fftLength)
        ampcor.setSearchWindowSizeWidth(searchMaxWidth)
        ampcor.setSearchWindowSizeHeight(searchMaxLength)
        ampcor.setAcrossLooks(1)
        ampcor.setDownLooks(1)
        ampcor.setOversamplingFactor(64)
        ampcor.setZoomWindowSize(16)
        ampcor.setAcrossGrossOffset(meanOffsetSamples)
        ampcor.setDownGrossOffset(meanOffsetLines)
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

        ampcorOffsetFile = 'ampcor.off'
        writeOffset(offsets, ampcorOffsetFile)

        #finalize image, and re-create it
        #otherwise the file pointer is still at the end of the image
        mSLC.finalizeImage()
        sSLC.finalizeImage()


#####################################################################
#   STEP 2. cull offsets
#####################################################################

        #first cull by mean
        threshold = 3.5
        meanOffsets = meanOffset(offsets)
        offsets = cullOffsetbyMean(offsets, meanOffsets[1], threshold)


        distances = (10,5,3,3,3,3,3,3)
        numCullOffsetsLimits = (100, 75, 50, 50, 50, 50, 50, 50)
    
        refinedOffsets = offsets
        for i, (distance, numCullOffsetsLimit) in enumerate(zip(distances, numCullOffsetsLimits)):

            cullOff = isceobj.createOffoutliers()
            cullOff.wireInputPort(name='offsets', object=refinedOffsets)
            cullOff.setSNRThreshold(2.0)
            cullOff.setDistance(distance)
        
            #set the tag used in the outfile. each message is precided by this tag
            #is the writer is not of "file" type the call has no effect
            stdWriter = create_writer("log", "", True, filename="offoutliers.log")
            stdWriter.setFileTag("offoutliers", "log")
            stdWriter.setFileTag("offoutliers", "err")
            stdWriter.setFileTag("offoutliers", "out")
            cullOff.setStdWriter(stdWriter)


            #run it
            cullOff.offoutliers()

            refinedOffsets = cullOff.getRefinedOffsetField()
            numLeft = len(refinedOffsets._offsets)
            print('Number of offsets left after %2dth culling: %5d'%(i, numLeft))
            if numLeft < numCullOffsetsLimit:
                print('******************************************************************')
                print('WARNING: There are not enough offsets left, so we are forced to')
                print('         use offset without culling')
                print('******************************************************************')
                refinedOffsets = offsets
                break
                #raise Exception('Too few points left after culling: {} left'.format(numLeft))

        cullOffsetFile = 'cull.off'
        writeOffset(refinedOffsets, cullOffsetFile)



#####################################################################
#   STEP 3. form interferogram and amplitdue image
#####################################################################
    
    nrlks = inps.rlks
    nalks = inps.alks
    intWidth = int(masterWidth / nrlks)

    if inps.master2 == None:
        master2 = inps.master
    else:
        master2 = inps.master2
    if inps.slave2 == None:
        slave2 = inps.slave
    else:
        slave2 = inps.slave2

    #re-create master
    mSLC = isceobj.createSlcImage()
    mSLC.setFilename(master2)
    mSLC.setWidth(masterWidth)
    mSLC.setLength(masterLength)
    mSLC.setAccessMode('read')
    mSLC.createImage()

    #re-create slave
    sSLC = isceobj.createSlcImage()
    sSLC.setFilename(slave2)
    sSLC.setWidth(slaveWidth)
    sSLC.setLength(slaveLength)
    sSLC.setAccessMode('read')
    sSLC.createImage()

    #create interferogram
    interf = isceobj.createIntImage()
    interf.setFilename(inps.interferogram)
    interf.setWidth(intWidth)
    interf.setAccessMode('write')
    interf.createImage()

    #create amplitdue
    amplitude = isceobj.createAmpImage()
    amplitude.setFilename(inps.amplitude)
    amplitude.setWidth(intWidth)
    amplitude.setAccessMode('write')
    amplitude.createImage()

    #create a writer for resamp
    stdWriter = create_writer("log", "", True, filename="resamp.log")
    stdWriter.setFileTag("resamp", "log")
    stdWriter.setFileTag("resamp", "err")
    stdWriter.setFileTag("resamp", "out")

    #open slave frame to get SAR parameters
    #with open(inps.slave + '.pck', 'rb') as f:
        #slaveFrame = pickle.load(f)


    #set up resampling program now
    #The setting has been compared with resamp_roi's setting in ROI_pac item by item.
    #The two kinds of setting are exactly the same. The number of setting items are
    #exactly the same
    objResamp = stdproc.createResamp()

    objResamp.wireInputPort(name='offsets', object=refinedOffsets)
    #objResamp.wireInputPort(name='instrument', object=slaveFrame.getInstrument())
    objResamp.stdWriter = stdWriter
    
    objResamp.setNumberFitCoefficients(6)

    objResamp.setNumberRangeBin1(masterWidth)
    objResamp.setNumberRangeBin2(slaveWidth)    
    objResamp.setStartLine(1)
    objResamp.setNumberLines(masterLength)
    objResamp.setFirstLineOffset(1)

    #check if this is consistent with slaveFrame.dopCoeff! No problem
    objResamp.setDopplerCentroidCoefficients(slaveFrame.dopCoeff)
    objResamp.setRadarWavelength(slaveFrame.radarWavelegth)
    objResamp.setSlantRangePixelSpacing(0.5 * SPEED_OF_LIGHT / slaveFrame.rangeSamplingRate)

    objResamp.setNumberRangeLooks(nrlks)
    objResamp.setNumberAzimuthLooks(nalks)

    objResamp.setFlattenWithOffsetFitFlag(0)

    #run it
    objResamp.resamp(mSLC, sSLC, interf, amplitude) 
    
    #finialize images
    mSLC.finalizeImage()
    sSLC.finalizeImage()
    interf.finalizeImage()
    amplitude.finalizeImage()

    #trim amplitude
    cmd = "imageMath.py -e='a_0*(a_1>0);a_1*(a_0>0)' --a={} -o={} -s BIP -t float".format(
        inps.amplitude, 
        'tmp.amp'
        )
    runCmd(cmd)
    os.remove(inps.amplitude)
    os.remove('tmp.amp.xml')
    os.remove('tmp.amp.vrt')
    os.rename('tmp.amp', inps.amplitude)
######################################################################


#example:
#no external offset provided
#match_resamp_isce.py -m 150515_ori.slc -s 151016_ori.slc -mframe 150515_ori.slc.pck -sframe 151016_ori.slc.pck -m2 150515.slc -s2 151016.slc -int 150515-151016.int -amp 150515-151016.amp -nr 40 -na 100 -rfw 64 -afw 512 -rsm 32 -asm 32


#with external offset provided:
#./match_resamp_isce.py -m 150515.slc -s 151016.slc -offset ../cull.off -int 150515-151016.int -amp 150515-151016.amp -sframe ../151016.slc.pck





