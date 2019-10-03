#!/usr/bin/env python3

#Cunren Liang, 04-JUN-2015
#JPL/Caltech

#######################################################################################
#program for calcualting the offsets of subswaths for wide swath SAR using slc images
#starting subswath can be any of the subswaths of a wide swath SAR except the last one
#any number of consecutive subswaths larger than two can be specified for mosaicking
#
#update history:
# 10-apr-2017, CL. support burst-by-burst processing
#######################################################################################


import os
import sys
import argparse
import datetime
import pickle
import ntpath
import shutil
import numpy as np
import scipy.signal as ss
from xml.etree.ElementTree import ElementTree

import isce
import isceobj
import mroipac
from mroipac.ampcor.Ampcor import Ampcor
from isceobj.Orbit.Orbit import Orbit
from iscesys.StdOEL.StdOELPy import create_writer
from isceobj.Constants import SPEED_OF_LIGHT

from crlpac import runCmd
from crlpac import create_xml
from crlpac import writeOffset
from crlpac import meanOffset
from crlpac import cullOffset
from crlpac import getOffset
from crlpac import cal_coherence


def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="calculate offsets of subswaths of wide swath SAR\
                    \n\
                    ")
    parser.add_argument('-s', action='append', dest='slc', default=[],
            help = 'subswath SLCs or amplitude images IN ORDER')
    parser.add_argument('-f', action='append', dest='frame', default=[],
            help = 'subswath frame pickle files IN ORDER')
    parser.add_argument('-o', dest='output', type=str, required=True,
            help = 'output txt file containing the offsets')
    parser.add_argument('-rlks', dest='rlks', type=int, default=1,
            help = 'number of range looks to take when matching. default: 1')
    parser.add_argument('-alks', dest='alks', type=int, default=10,
            help = 'number of azimuth looks to take when matching. default: 10')

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

    ns = len(slcs)

    ############################
    pcks = []
    
    #dopCoeff = []
    #orbit = []
    
    width = []
    length = []
    startingRange = []
    farRange = []
    sensingStart = []
    sensingStop = []

    azimuthLineInterval = [] 
    rangePixelSize = []
    ############################

    for frame in frames:
        with open(frame, 'rb') as fid:
            pcks.append(pickle.load(fid))

    for pck in pcks:
        width.append(pck.getNumberOfSamples())
        length.append(pck.getNumberOfLines())
        startingRange.append(pck.startingRange)
        farRange.append(pck.getFarRange())
        sensingStart.append(pck.getSensingStart())
        sensingStop.append(pck.getSensingStop())
        azimuthLineInterval.append(1.0/pck.PRF)
        rangePixelSize.append(0.5 * SPEED_OF_LIGHT / pck.rangeSamplingRate)
    
    

    #get offset from gemoetrical parameters
    #measured by subswath 0 pixel size
    rangeOffset1 = []
    azimuthOffset1 = []
    rangeOffset1.append(0.0)
    azimuthOffset1.append(0.0)
    for i in range(1, ns):
        rangeOffset1.append(-(startingRange[i] - startingRange[i-1]) / rangePixelSize[0]) 
        azTimeDiff = (sensingStart[i] - sensingStart[i-1]).total_seconds()
        azimuthOffset1.append(-azTimeDiff / azimuthLineInterval[0])


    #################################################
    #prepare magnitude image for matching
    
    #number of azimuth looks to be taken: 5
    #this is approximately appropriate for ALOS2 full-aperture products
    #to suppress the spikes and achieve an appropriate oversampling rate
    #considering its azimuth resolution.
    #some discussions can be found in:
    #R. Bamler and M. Eineder, ScanSAR Processing Using Standard High Precision
    #SAR Algorithms, IEEE Trans. Geosci. Remote Sens., vol. 43, no. 1, Jan. 1996
    nalks = inps.alks
    nrlks = inps.rlks


    rangeScale = []
    azimuthScale = []
    rectWidth = []
    rectLength = []
    lookRectWidth = []
    lookRectLength = []
    lookRectMag = []
    for i, slc in enumerate(slcs):
        #STEP 1. get power
        base = os.path.splitext(ntpath.basename(slc))[0]
        powName = "{}_{}.pow".format(base, i)
        
        #works for both *.slc and *.amp (burst-by-burst processing)
        cmd = 'imageMath.py -e="abs(a)*abs(a)" --a={} -o {}'.format(slc, powName)
        runCmd(cmd)
        #resample to the sample sizes of subswath 1

        #STEP 2. rectify power image
        rangeScale.append(rangePixelSize[0] / rangePixelSize[i])
        azimuthScale.append(azimuthLineInterval[0] / azimuthLineInterval[i])
        rectPowName = "rect_{}".format(powName)
        if i == 0:
            rectWidth.append( width[i] )
            rectLength.append( length[i] )
            cmd = 'mv {} {}'.format(powName, rectPowName)
            runCmd(cmd)
        else:
            rectWidth.append( int(1.0 / rangeScale[i] * width[i]) )
            rectLength.append( int(1.0 / azimuthScale[i] * length[i]) )

            rectInput = '''Input Image File Name  (-) = {inFile}        ! dimension of file to be rectified
Output Image File Name (-) = {outFile}       ! dimension of output
Input Dimensions       (-) = {inWidth} {inLength} ! across, down
Output Dimensions      (-) = {outWidth} {outLength} ! across, down
Affine Matrix Row 1    (-) = {m11} {m12}      ! a b
Affine Matrix Row 2    (-) = {m21} {m22}      ! c d
Affine Offset Vector   (-) = {t1} {t2}        ! e f
Looks of the Offsets   (-) = {rlks} {alks}        ! lac, ldn
Looks of the Image File   (-) = 1 1        ! lac0, ldn0
File Type              (-) = {format}        ! [RMG, COMPLEX]
Interpolation Method   (-) = {method}        ! [NN, Bilinear, Sinc]
'''.format(
            inFile = powName,
            outFile = rectPowName,
            inWidth = width[i], inLength = length[i],
            outWidth = rectWidth[i], outLength = rectLength[i],
            m11 = rangeScale[i], m12 = 0.0000000000,
            m21 = 0.0000000000, m22 = azimuthScale[i],
            t1 = 0.0, t2 = 0.0,
            rlks = 1, alks = 1,
            format = 'REAL',
            method = 'Bilinear'
            )

            rectInputFile = 'rect_with_looks_{}.in'.format(i)
            with open(rectInputFile, 'w') as ff:
                ff.write(rectInput)
            cmd = '$INSAR_ZERODOP_BIN/rect_with_looks {}'.format(rectInputFile)
            runCmd(cmd)

        #create xml file
        create_xml(rectPowName, rectWidth[i], rectLength[i], 'float')

        #STEP 3. take looks
        lookRectPowName = '{}_{}rlks_{}alks.pow'.format(os.path.splitext(rectPowName)[0], nrlks, nalks)
        lookRectWidth.append( int(rectWidth[i] /  nrlks) )
        lookRectLength.append( int(rectLength[i] / nalks) ) 
        if inps.rlks == 1 and inps.alks == 1:
            os.symlink(rectPowName, lookRectPowName)
            create_xml(lookRectPowName, lookRectWidth[i], lookRectLength[i], 'float')
        else:
            cmd = 'looks.py -i {} -o {} -r {} -a {}'.format(rectPowName, lookRectPowName, nrlks, nalks)
            runCmd(cmd)

        #STEP 4. get magnitude image
        lookRectMagName = "{}.mag".format(os.path.splitext(lookRectPowName)[0])
        cmd = 'imageMath.py -e="sqrt(a)" --a={} -o {}'.format(lookRectPowName, lookRectMagName)
        runCmd(cmd)
        lookRectMag.append(lookRectMagName)


    #get offset from matching
    rangeOffset2 = []
    azimuthOffset2 = []
    rangeOffset2.append(0.0)
    azimuthOffset2.append(0.0)
    #do the matching between two consecutive subswaths:
    for i in range(1, ns):

        mMag = isceobj.createImage()
        mMag.setFilename(lookRectMag[i-1])
        mMag.dataType = 'FLOAT'
        mMag.setWidth(lookRectWidth[i-1])
        mMag.setLength(lookRectLength[i-1])
        mMag.setAccessMode('read')
        mMag.createImage()

        sMag = isceobj.createImage()
        sMag.setFilename(lookRectMag[i])
        sMag.dataType = 'FLOAT'
        sMag.setWidth(lookRectWidth[i])
        sMag.setLength(lookRectLength[i])
        sMag.setAccessMode('read')
        sMag.createImage()

        #still keep the name here: insarapp_slcs_ampcor
        ampcor = Ampcor(name='insarapp_slcs_ampcor')
        ampcor.configure()

        #DATA TYPE
        ampcor.setImageDataType1('real')
        ampcor.setImageDataType2('real')

        #INPUT/OUTPUT FILES
        ampcor.setMasterSlcImage(mMag)
        ampcor.setSlaveSlcImage(sMag)

        #MATCH REGION
        ########################################
        numRange = 20;
        numAzimuth = 100;
        firstSample = 1
        firstLine = 1
        rangeMeanOffset = int(rangeOffset1[i]/nrlks)
        azimuthMeanOffset = int(azimuthOffset1[i]/nalks)

        #it seems that we cannot use 0, don't know why
        if rangeMeanOffset == 0:
            rangeMeanOffset = 1
        if azimuthMeanOffset == 0:
            azimuthMeanOffset = 1

        if rangeMeanOffset < 0:
            firstSample = int(35 - rangeMeanOffset)
        if azimuthMeanOffset < 0:
            firstLine = int(35 - azimuthMeanOffset)
        ########################################
        ampcor.setFirstSampleAcross(firstSample)
        ampcor.setLastSampleAcross(lookRectWidth[i-1])
        ampcor.setNumberLocationAcross(numRange)
        ampcor.setFirstSampleDown(firstLine)
        ampcor.setLastSampleDown(lookRectLength[i-1])
        ampcor.setNumberLocationDown(numAzimuth)

        #MATCH PARAMETERS
        ########################################
        fftWidth = 32
        fftLength = 32
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
        ampcor.setAcrossGrossOffset(rangeMeanOffset)
        ampcor.setDownGrossOffset(azimuthMeanOffset)
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
        #offsetFile = 'offset_{}_{}.off'.format(i-1, i)
        #cullOffsetFile = 'cull_{}_{}.off'.format(i-1, i)
        #dumpFile = 'fitoff_{}_{}.out'.format(i-1, i)
        #atParam = getOffset(offsets, offsetFile, cullOffsetFile, dumpFile)
        #rangeOffset2.append(atParam[4])
        #azimuthOffset2.append(atParam[5])

        #record offset
        #ampcorOffsetFile = 'ampcor_{}_{}.off'.format(i-1, i)
        #writeOffset(offsets, ampcorOffsetFile)

        #cull offsets
        distances = (10,5,3)
        numCullOffsetsLimits = (100, 75, 50)
        refinedOffsets = cullOffset(offsets, distances, numCullOffsetsLimits)

        #record new offsets
        #cullOffsetFile = 'cull_{}_{}.off'.format(i-1, i)
        #writeOffset(refinedOffsets, cullOffsetFile)

        #get mean offset    
        #can use getFitPolynomials
        #test with mean!!!!!!!!!!
        #!!!!!!!!!!!!!!!!!!!!!!!!!!
        if refinedOffsets != None:
            rangeOffset2.append(meanOffset(refinedOffsets)[0]*nrlks)
            azimuthOffset2.append(meanOffset(refinedOffsets)[1]*nalks)
            #print("++++++++{}, {}, {}".format(i, rangeOffset2[i], azimuthOffset2[i]))
        else:
            print('******************************************************************')
            print('WARNING: There are not enough offsets left, so we are forced to')
            print('         use geometrical offset for subswath mosaicking')
            print('******************************************************************')
            rangeOffset2.append(rangeOffset1[i])
            azimuthOffset2.append(azimuthOffset1[i])


    offsetComp = "\n\ncomparision of offsets:\n\n"
    offsetComp += "offset type       i     geometrical           match      difference\n"
    offsetComp += "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
    for i, (offset1, offset2) in enumerate(zip(rangeOffset1, rangeOffset2)):
        offsetComp += "range offset     {:2d}   {:13.3f}   {:13.3f}   {:13.3f}\n".format(i, offset1, offset2, offset1 - offset2)
    for i, (offset1, offset2) in enumerate(zip(azimuthOffset1, azimuthOffset2)):
        offsetComp += "azimuth offset   {:2d}   {:13.3f}   {:13.3f}   {:13.3f}\n".format(i, offset1, offset2, offset1 - offset2)

    #write and report offsets
    with open(inps.output, 'w') as f:
        f.write(offsetComp)
    print("{}".format(offsetComp))


#######################################################################################

    #tidy up: delete less important files
    for i, slc in enumerate(slcs):
        #
        base = os.path.splitext(ntpath.basename(slc))[0]
        powName = "{}_{}.pow".format(base, i)
        rectPowName = "rect_{}".format(powName)
        lookRectPowName = '{}_{}rlks_{}alks.pow'.format(os.path.splitext(rectPowName)[0], nrlks, nalks)
        lookRectMagName = "{}.mag".format(os.path.splitext(lookRectPowName)[0])
        
        #for the first subswath, power image is renamed as rectified power image
        if i != 0:
            os.remove(powName)
        os.remove(rectPowName)
        os.remove(lookRectPowName)
        os.remove(lookRectMagName)

        os.remove(powName + '.xml')
        os.remove(rectPowName + '.xml')
        os.remove(lookRectPowName + '.xml')
        os.remove(lookRectMagName + '.xml')

        os.remove(powName + '.vrt')
        os.remove(rectPowName + '.vrt')
        os.remove(lookRectPowName + '.vrt')
        os.remove(lookRectMagName + '.vrt')

    for i in range(1, ns):
        rectInputFile = 'rect_with_looks_{}.in'.format(i)
        os.remove(rectInputFile)


