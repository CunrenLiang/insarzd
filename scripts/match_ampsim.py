#!/usr/bin/env python3

#Cunren Liang, 5-MAY-2015
#JPL/Caltech

import os
import sys
import shutil
import argparse
import pickle

import isce
import isceobj
from isceobj.Location.Offset import OffsetField,Offset
from mroipac.ampcor.Ampcor import Ampcor

from crlpac import getWidth
from crlpac import getLength
from crlpac import runCmd
from crlpac import cullOffset
from crlpac import getOffset


def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser( description='matching between radar and simulation')
    parser.add_argument('-m', '--amp', dest='amp', type=str, required=True,
            help = 'amplitude image')
    parser.add_argument('-s', '--sim', dest='sim', type=str, required=True,
            help = 'simulated radar image')
    parser.add_argument('-f', '--aff', dest='aff', type=str, required=True,
            help = '(output) affine transformation')
    parser.add_argument('-r','--rlks', dest='rlks', type=int, default=1,
            help = 'range looks')
    parser.add_argument('-a', '--alks', dest='alks', type=int, default=1,
            help = 'azimuth looks')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':
    
    inps = cmdLineParse()

    ampWidth = getWidth(inps.amp + '.xml')
    ampLength = getLength(inps.amp + '.xml')
    lookAmpWidth = int(ampWidth/inps.rlks)
    lookAmpLength = int(ampLength/inps.alks)

    simWidth = getWidth(inps.sim + '.xml')
    simLength = getLength(inps.sim + '.xml')
    lookSimWidth = int(simWidth/inps.rlks)
    lookSimLength = int(simLength/inps.alks)

    #prepare parameters for ampcor
    lookAmp = 'float_{}rlks_{}alks.amp'.format(inps.rlks, inps.alks)
    lookSim = 'float_{}rlks_{}alks.sim'.format(inps.rlks, inps.alks)
    offsetFile = 'ampsim.off'

    numAzimuth = 30 #number of matches in azimuth
    numRange = 30  #number of matches in range
    skipLines = int(lookAmpLength / numAzimuth)
    skipSamples = int(lookAmpWidth / numRange)
    fftWidth = 64 #matching window width
    fftLength = 64  #matching window length
    searchMaxWidth = 60
    searchMaxLength = 60
    s1 = 0
    s2 = 1

    #take looks and prepare data here!!!!!!!!!
    #do matching at matchdir
    matchdir = 'match_ampsim'
    os.mkdir(matchdir)
    os.chdir(matchdir)

    tmpAmp = 'tmp'
    cmd = '$INSAR_ZERODOP_SCR/look.py -i {} -o {} -r {} -a {}'.format('../' + inps.amp, tmpAmp, inps.rlks, inps.alks)
    runCmd(cmd)
    cmd = "imageMath.py -e='sqrt(a_0*a_0+a_1*a_1)' --a={} -o {} -t float".format(tmpAmp, lookAmp)
    runCmd(cmd)
    os.remove(tmpAmp)
    os.remove(tmpAmp + '.xml')
    os.remove(tmpAmp + '.vrt')

    if inps.rlks != 1 and inps.alks != 1:
        cmd = 'looks.py -i {} -o {} -r {} -a {}'.format('../' + inps.sim, lookSim, inps.rlks, inps.alks)
        runCmd(cmd)
    else:
        cmd = "cp {} {}".format('../' + inps.sim, lookSim)
        runCmd(cmd)
        cmd = "cp {} {}".format('../' + inps.sim + '.xml', lookSim + '.xml')
        runCmd(cmd)
        cmd = "cp {} {}".format('../' + inps.sim + '.vrt', lookSim + '.vrt')
        runCmd(cmd)

        img = isceobj.createImage()
        img.load(lookSim + '.xml')
        img.filename = lookSim
        img.extraFilename = lookSim + '.vrt'
        img.renderHdr()


#############################################################################################################
#   Line in match_ampsim_roipac.py in replaced by the following. Other part remains the same.

    #set amp image
    objAmp = isceobj.createImage()
    objAmp.setFilename(lookAmp)
    objAmp.setWidth(lookAmpWidth)
    objAmp.setLength(lookAmpLength)
    objAmp.dataType='FLOAT'
    objAmp.setAccessMode('read')
    objAmp.createImage()

    #set sim image
    objSim = isceobj.createImage()
    objSim.setFilename(lookSim)
    objSim.setWidth(lookSimWidth)
    objSim.setLength(lookSimLength)
    objSim.dataType='FLOAT'
    objSim.setAccessMode('read')
    objSim.createImage()

    objAmpcor = Ampcor(name='insarapp_intsim_ampcor')
    objAmpcor.configure()

    #DATA TYPE
    objAmpcor.setImageDataType1('real')
    objAmpcor.setImageDataType2('real')

    #INPUT/OUTPUT FILES
    objAmpcor.setMasterSlcImage(objAmp)
    objAmpcor.setSlaveSlcImage(objSim)

    #MATCH REGION
    ########################################
    #avoid the complain of Ampcor.py
    xMargin = 2 * searchMaxWidth + fftWidth
    yMargin = 2 * searchMaxLength + fftLength
    #make it smaller
    xMargin = fftWidth / 2 + 5
    yMargin = fftLength / 2 + 5

    firstSampleAcross = xMargin
    firstSampleDown = yMargin
    lastSampleAcross = lookAmpWidth - xMargin
    lastSampleDown = lookAmpLength - yMargin
    ########################################
    objAmpcor.setFirstSampleAcross(firstSampleAcross)
    objAmpcor.setLastSampleAcross(lastSampleAcross)
    objAmpcor.setNumberLocationAcross(numRange)

    objAmpcor.setFirstSampleDown(firstSampleDown)
    objAmpcor.setLastSampleDown(lastSampleDown)
    objAmpcor.setNumberLocationDown(numAzimuth)

    #MATCH PARAMETERS
    objAmpcor.setWindowSizeWidth(fftWidth)
    objAmpcor.setWindowSizeHeight(fftLength)
    objAmpcor.setSearchWindowSizeWidth(searchMaxWidth)
    objAmpcor.setSearchWindowSizeHeight(searchMaxLength)
    objAmpcor.setAcrossLooks(1)
    objAmpcor.setDownLooks(1)
    objAmpcor.setOversamplingFactor(64)
    objAmpcor.setZoomWindowSize(16)
    objAmpcor.setAcrossGrossOffset(0)
    objAmpcor.setDownGrossOffset(0)
    #1. The following not set
    #Matching Scale for Sample/Line Directions                       (-)    = 1. 1.
    #should add the following in Ampcor.py?
    #if not set, in this case, Ampcor.py'value is also 1. 1.
    #objAmpcor.setScaleFactorX(1.)
    #objAmpcor.setScaleFactorY(1.)



    #MATCH THRESHOLDS AND DEBUG DATA
    #2. The following not set
    #in roi_pac the value is set to 0 1
    #in isce the value is set to 0.001 1000.0
    #SNR and Covariance Thresholds                                   (-)    =  {s1} {s2}
    #should add the following in Ampcor?
    #THIS SHOULD BE THE ONLY THING THAT IS DIFFERENT FROM THAT OF ROI_PAC
    #objAmpcor.setThresholdSNR(0)
    #objAmpcor.setThresholdCov(1)
    objAmpcor.setDebugFlag(False)
    objAmpcor.setDisplayFlag(False)

    #in summary, only two things not set which are indicated by 'The following not set' above.
    
    #run ampcor
    objAmpcor.ampcor()

    #get offsets
    offsets = objAmpcor.getOffsetField()


    offsetsPlain = ''
    for offsetx in offsets:
        offsetsPlainx = "{}".format(offsetx)
        offsetsPlainx = offsetsPlainx.split()
        offsetsPlain = offsetsPlain + "{:8d} {:10.3f} {:8d} {:12.3f} {:11.5f} {:11.6f} {:11.6f} {:11.6f}\n".format(
            int(offsetsPlainx[0]),
            float(offsetsPlainx[1]),
            int(offsetsPlainx[2]),
            float(offsetsPlainx[3]),
            float(offsetsPlainx[4]),
            float(offsetsPlainx[5]),
            float(offsetsPlainx[6]),
            float(offsetsPlainx[7])
            )

    with open(offsetFile, 'w') as f:
        f.write(offsetsPlain)

#############################################################################################################

    cullOffsetFile = 'ampsim_cull.off'
    dumpFile = inps.aff
    #run fitoff here
    cmd = '$INSAR_ZERODOP_BIN/fitoff {} {} 1.5 .5 50 > {}'.format(offsetFile, cullOffsetFile, dumpFile)
    print("{}".format(cmd))
    runCmd(cmd)

    #check number of matching points left
    with open(cullOffsetFile, 'r') as ff:
        numCullOffsets = sum(1 for linex in ff)
        ff.close
    if numCullOffsets < 50:
        raise Exception('Too few points left after culling, {} left'.format(numCullOffsets))

    #cp it to the parent directory
    shutil.copy(dumpFile, '../' + dumpFile)
#############################################################################################################


#./match_ampsim.py -m 141018-141123.amp -s 141018-141123.sim -f ampsim_16rlks_16alks.aff -r 16 -a 16

