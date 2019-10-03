#!/usr/bin/env python3

#Cunren Liang, 04-JUN-2015
#JPL/Caltech

########################################################################################
#program for mosaicking wide swath SAR subswath interferograms and amplitude images
#starting subswath can be any of the wide swath SAR subswaths except the last one
#any number of consecutive subswaths larger than two can be specified for mosaicking

#convention: paramters of frame pickle file are single-look parameters
#            all the range/azimuth number of looks of the subswath interferograms/amplitudes are the same

#update history:
# 12-APR-2016, CL. This version compensates for phase difference for subsequent
#                  interferograms to make them consistent with subswath 1 interferogram
#
# 10-apr-2017, CL. support burst-by-burst processing
# 24-apr-2017, CL. support mosaicking multiple-look interferograms and amplitudes 
########################################################################################


import os
import re
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
            epilog="mosaic subswaths of wide swath SAR\
                    \n\
                    ")
    parser.add_argument('-f', action='append', dest='frame', default=[],
            help = 'subswath frame pickle files of master IN ORDER')
    parser.add_argument('-i', action='append', dest='inf', default=[],
            help = 'subswath interferograms IN ORDER')
    parser.add_argument('-a', action='append', dest='amp', default=[],
            help = 'subswath amplitudes IN ORDER')
    parser.add_argument('-offset', dest='offset', type=str, required=True,
            help = 'subswath offset file')
    parser.add_argument('-io', dest='infoutput', type=str, required=True,
            help = 'output interferogram')
    parser.add_argument('-ao', dest='ampoutput', type=str, required=True,
            help = 'output amplitude image')
    parser.add_argument('-fo', dest='frameoutput', type=str, required=True,
            help = 'output frame pickle file for mosaicked image. also generate *.par.xml file based on this')
    parser.add_argument('-imethod', dest='infmethod', type=int, default=3,
            help = 'interferogram mosaicking method: 0: mosaic at the center of overlap area. 1: mosaic at the right egde of overlap area. 2: mosaic at the left edge of overlap area. 3: add overlap area (default). 4: Liang weighting method')
    parser.add_argument('-amethod', dest='ampmethod', type=int, default=0,
            help = 'amplitude mosaicking method: 0: mosaic at the center of overlap area (default). 1: mosaic at the right egde of overlap area. 2: mosaic at the left edge of overlap area. 3: add overlap area. 4: Liang weighting method')
    parser.add_argument('-rlks', dest='rlks', type=int, default=1,
            help = 'number of range looks of the interferograms and amplitudes. default: 1')
    parser.add_argument('-alks', dest='alks', type=int, default=1,
            help = 'number of azimuth looks of the interferograms and amplitudes. default: 1')
    parser.add_argument('-phc', dest='phc', type=int, default=0,
            help = 'whether do phase correction for subswath interferograms (except the first subswath). 0: No (default). 1: Yes.')
    parser.add_argument('-drlks', dest='drlks', type=int, default=1,
            help = 'number of range looks of the differential interferogram in the overlap area. default: 1')
    parser.add_argument('-dalks', dest='dalks', type=int, default=1,
            help = 'number of azimuth looks of the differential interferogram in the overlap area. default: 1')
    parser.add_argument('-filt', dest='filt', type=int, default=0,
            help = 'whether do filtering for interferograms used for subswath phase difference calculation. 0: No (default). 1: Yes.')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':

    BIN_DIR='$INSAR_ZERODOP_BIN'

    inps = cmdLineParse()

    frames = inps.frame
    infs = inps.inf
    amps = inps.amp
    rlks = inps.rlks
    alks = inps.alks

    if not (len(frames) == len(infs) == len(amps)):
        raise Exception('the numbers of files are not equal!')
    if len(frames) < 2:
        raise Exception('at least 2 subswaths should be specified!')

    ns = len(frames)

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
        with open(frame, 'rb') as f:
            pcks.append(pickle.load(f))

    for pck in pcks:
        width.append(pck.getNumberOfSamples())
        length.append(pck.getNumberOfLines())
        startingRange.append(pck.startingRange)
        farRange.append(pck.getFarRange())
        sensingStart.append(pck.getSensingStart())
        sensingStop.append(pck.getSensingStop())
        azimuthLineInterval.append(1.0/pck.PRF)
        rangePixelSize.append(0.5 * SPEED_OF_LIGHT / pck.rangeSamplingRate)

    rangeScale = []
    azimuthScale = []
    rectWidth = []
    rectLength = []
    for i in range(ns):
        rangeScale.append(rangePixelSize[0] / rangePixelSize[i])
        azimuthScale.append(azimuthLineInterval[0] / azimuthLineInterval[i])
        if i == 0:
            rectWidth.append( int(width[i] / rlks) )
            rectLength.append( int(length[i] / alks) )
        else:
            rectWidth.append( int(1.0 / rangeScale[i] * int(width[i] / rlks)) )
            rectLength.append( int(1.0 / azimuthScale[i] * int(length[i] / alks)) )


    #read offset
    #################################
    #offset from gemoetrical parameters
    rangeOffset1 = []
    azimuthOffset1 = []
    #offset from matching
    rangeOffset2 = []
    azimuthOffset2 = []

    with open(inps.offset, 'r') as f:
        lines = f.readlines()
    for linex in lines:
        if 'range offset' in linex:
            rangeOffset1.append(float(re.split('\s+', linex)[3]))
            rangeOffset2.append(float(re.split('\s+', linex)[4]))
        if 'azimuth offset' in linex:
            azimuthOffset1.append(float(re.split('\s+', linex)[3]))
            azimuthOffset2.append(float(re.split('\s+', linex)[4]))
    #################################

    #convert original offset to offset for images with looks
    rangeOffset1 = list(np.asarray(rangeOffset1)/rlks)
    azimuthOffset1 = list(np.asarray(azimuthOffset1)/alks)
    rangeOffset2 = list(np.asarray(rangeOffset2)/rlks)
    azimuthOffset2 = list(np.asarray(azimuthOffset2)/alks)

    #get offset relative to the first frame
    rangeOffset = []
    azimuthOffset = []
    rangeOffset.append(0.0)
    azimuthOffset.append(0.0)
    for i in range(1, ns):
        rangeOffset.append(0.0)
        azimuthOffset.append(0.0)
        for j in range(1, i+1):
            rangeOffset[i] += rangeOffset2[j]
            #here we still use azimuth offset calculated from geometrical parameters
            #as it should be accurate
            azimuthOffset[i] += azimuthOffset1[j]


    #resample each subswath
    rinfs = []
    ramps = []
    for i, (inf, amp) in enumerate(zip(infs, amps)):
        infBase = os.path.splitext(ntpath.basename(inf))[0]
        ampBase = os.path.splitext(ntpath.basename(amp))[0]
        rinfs.append("{}_{}.int".format(infBase, i))
        ramps.append("{}_{}.amp".format(ampBase, i))

        if i == 0:
            #shutil.copyfile(inf, rinfs[i])
            #shutil.copyfile(amp, ramps[i])
            os.symlink(inf, rinfs[i])
            os.symlink(amp, ramps[i])
            continue

        rangeOffsetFrac = rangeOffset[i] - int(rangeOffset[i])
        azimuthOffsetFrac = azimuthOffset[i] - int(azimuthOffset[i])


        rectInputInf = '''Input Image File Name  (-) = {inFile}        ! dimension of file to be rectified
Output Image File Name (-) = {outFile}       ! dimension of output
Input Dimensions       (-) = {inWidth} {inLength} ! across, down
Output Dimensions      (-) = {outWidth} {outLength} ! across, down
Affine Matrix Row 1    (-) = {m11} {m12}      ! a b
Affine Matrix Row 2    (-) = {m21} {m22}      ! c d
Affine Offset Vector   (-) = {t1} {t2}        ! e f
Looks of the Offsets   (-) = 1 1        ! lac, ldn
Looks of the Image File   (-) = 1 1        ! lac0, ldn0
File Type              (-) = {format}        ! [RMG, COMPLEX]
Interpolation Method   (-) = {method}        ! [NN, Bilinear, Sinc]
'''.format(
        inFile = inf,
        outFile = rinfs[i],
        inWidth = int(width[i] / rlks), inLength = int(length[i] / alks),
        outWidth = rectWidth[i], outLength = rectLength[i],
        m11 = rangeScale[i], m12 = 0.0000000000,
        m21 = 0.0000000000, m22 = azimuthScale[i],
        t1 = rangeOffsetFrac * rangeScale[i], t2 = azimuthOffsetFrac * azimuthScale[i],
        format = 'COMPLEX',
        method = 'Sinc'
        )

        rectInputAmp = '''Input Image File Name  (-) = {inFile}        ! dimension of file to be rectified
Output Image File Name (-) = {outFile}       ! dimension of output
Input Dimensions       (-) = {inWidth} {inLength} ! across, down
Output Dimensions      (-) = {outWidth} {outLength} ! across, down
Affine Matrix Row 1    (-) = {m11} {m12}      ! a b
Affine Matrix Row 2    (-) = {m21} {m22}      ! c d
Affine Offset Vector   (-) = {t1} {t2}        ! e f
Looks of the Offsets   (-) = 1 1        ! lac, ldn
Looks of the Image File   (-) = 1 1        ! lac0, ldn0
File Type              (-) = {format}        ! [RMG, COMPLEX]
Interpolation Method   (-) = {method}        ! [NN, Bilinear, Sinc]
'''.format(
        inFile = amp,
        outFile = ramps[i],
        inWidth = int(width[i] / rlks), inLength = int(length[i] / alks),
        outWidth = rectWidth[i], outLength = rectLength[i],
        m11 = rangeScale[i], m12 = 0.0000000000,
        m21 = 0.0000000000, m22 = azimuthScale[i],
        t1 = rangeOffsetFrac * rangeScale[i], t2 = azimuthOffsetFrac * azimuthScale[i],
        format = 'COMPLEX',
        method = 'Sinc'
        )

        rectInputFileInf = 'rect_with_looks_inf_{}.in'.format(i)
        with open(rectInputFileInf, 'w') as f:
            f.write(rectInputInf)
        cmd = '$INSAR_ZERODOP_BIN/rect_with_looks {}'.format(rectInputFileInf)
        runCmd(cmd)

        rectInputFileAmp = 'rect_with_looks_amp_{}.in'.format(i)
        with open(rectInputFileAmp, 'w') as f:
            f.write(rectInputAmp)
        cmd = '$INSAR_ZERODOP_BIN/rect_with_looks {}'.format(rectInputFileAmp)
        runCmd(cmd)


    #determine output width and length
    #actually no need to calculate in range direction
    xs = []
    xe = []
    ys = []
    ye = []
    for i in range(ns):
        if i == 0:
            xs.append(0)
            xe.append(rectWidth[i] - 1)
            ys.append(0)
            ye.append(rectLength[i] - 1)
        else:
            xs.append(0 - int(rangeOffset[i]))
            xe.append(rectWidth[i] - 1 - int(rangeOffset[i]))
            ys.append(0 - int(azimuthOffset[i]))
            ye.append(rectLength[i] - 1 - int(azimuthOffset[i]))

    (xmin, xminIndex) = min((v,i) for i,v in enumerate(xs))
    (xmax, xmaxIndex) = max((v,i) for i,v in enumerate(xe))
    (ymin, yminIndex) = min((v,i) for i,v in enumerate(ys))
    (ymax, ymaxIndex) = max((v,i) for i,v in enumerate(ye))

    outWidth = xmax - xmin + 1
    outLength = ymax - ymin + 1


    #prepare offset for mosaicing
    rangeOffset3 = []
    azimuthOffset3 = []
    for i in range(ns):
        azimuthOffset3.append(int(azimuthOffset[i]) - int(azimuthOffset[yminIndex]))
        if i != 0:
            rangeOffset3.append(int(rangeOffset[i]) - int(rangeOffset[i-1]))
        else:
            rangeOffset3.append(0)


    #mosaic interferograms
    delta = int(30 / rlks)
    diffflag = 0
    oflag = inps.infmethod
    infOutput = inps.infoutput
    cmd = "{}/mosaicsubswath {} {} {} {} {} {}".format(BIN_DIR, infOutput, outWidth, outLength, delta, diffflag, ns)
    for i in range(ns):
        cmdi = " {} {} {} {} {}".format(rinfs[i], rectWidth[i], rangeOffset3[i], azimuthOffset3[i], oflag) 
        cmd += cmdi
    runCmd(cmd)

    #re-run mosaicsubswath using phase-compensated interferograms later
    #run this after compensating interferograms
    diffflag = 1 #no need to output difference for re-running
    cmdbak = "{}/mosaicsubswath {} {} {} {} {} {}".format(BIN_DIR, infOutput, outWidth, outLength, delta, diffflag, ns)
    for i in range(ns):
        cmdi = " {} {} {} {} {}".format(rinfs[i], rectWidth[i], rangeOffset3[i], azimuthOffset3[i], oflag) 
        cmdbak += cmdi

    #mosaic amplitudes
    diffflag = 1  #do not output overlap area difference for amplitude image
    oflag = inps.ampmethod
    ampOutput = inps.ampoutput
    cmd = "{}/mosaicsubswath {} {} {} {} {} {}".format(BIN_DIR, ampOutput, outWidth, outLength, delta, diffflag, ns)
    for i in range(ns):
        cmdi = " {} {} {} {} {}".format(ramps[i], rectWidth[i], rangeOffset3[i], azimuthOffset3[i], oflag) 
        cmd += cmdi
    runCmd(cmd)

    #update pickle file
    #should update doppler for the whole scene, but it won't be used in the following interferometric
    #processing, so I don't do this here
    #no need to update orbit, as we are mosaicking subswaths
    #assume there are outWidth*rlks samples, and outLength*alks lines
    outStartingRange = startingRange[0]
    outFarRange = outStartingRange + (outWidth*rlks - 1) * rangePixelSize[0]
    outSensingStart = sensingStart[0] + datetime.timedelta(seconds=azimuthOffset3[0]*alks * azimuthLineInterval[0])
    outSensingStop = outSensingStart + datetime.timedelta(seconds=(outLength*alks - 1) * azimuthLineInterval[0])

    outpck = pcks[0]
    outpck.setNumberOfSamples(outWidth*rlks)
    outpck.setNumberOfLines(outLength*alks)
    outpck.startingRange = outStartingRange
    outpck.setFarRange(outFarRange)
    outpck.setSensingStart(outSensingStart)
    outpck.setSensingStop(outSensingStop)
    #sensing mid
    outpck.setSensingMid(outpck.getSensingStart() + datetime.timedelta(seconds=(outpck.getSensingStop()-outpck.getSensingStart()).total_seconds()/2.0))

    with open(inps.frameoutput, 'wb') as f:
        pickle.dump(outpck, f)

    #create XML file for mosaicked interferogram
    create_xml(infOutput, outWidth, outLength, 'int')
    #create XML file for mosaicked amplitude
    create_xml(ampOutput, outWidth, outLength, 'amp')


#######################################################################################
#adjust phase of subsequent subswath interferograms here
#WARNING: 1. NOW ONLY SUPPORT THE CASE THAT ALL SUBSWATH RANGE SAMPLING INTERVALS ARE THE
#            SAME.
#         2. phase difference calculation should also be able to work for ScanSAR-stripmap
#            interferometry, although the overlap area may not be the center of 
#            stripmap swath.


    #parameters for this part
    nrlks2 = inps.drlks
    nalks2 = inps.dalks
    corth = 0.85

    if inps.filt==1:
        corth = 0.90

    diffMean = []
    diffMean.append(0.0)
    for i in range(ns - 1):
        diffFile = '{}-{}.int'.format(i, i+1)
        diffFile1 = '{}-{}_{}.int'.format(i, i+1, i)
        diffFile2 = '{}-{}_{}.int'.format(i, i+1, i+1)
        diffLookFile = '{}-{}_{}rlks_{}alks.int'.format(i, i+1, nrlks2, nalks2)
        diffLookFile1 = '{}-{}_{}_{}rlks_{}alks.int'.format(i, i+1, i, nrlks2, nalks2)
        diffLookFile2 = '{}-{}_{}_{}rlks_{}alks.int'.format(i, i+1, i+1, nrlks2, nalks2)

        diffWidth = rectWidth[i] + rangeOffset3[i+1] - 2 * delta
        diffLength = os.stat(diffFile).st_size / 8 / diffWidth
        #to make it integer. crl sep-10-2017
        #diffLength /= 2
        diffLength = int(diffLength/2)

        #STEP 1. extract subswath interferograms of overlap area
        diff1 = (np.fromfile(diffFile, dtype=np.complex64).reshape(diffLength*2, diffWidth))[0:diffLength*2:2, :]
        diff2 = (np.fromfile(diffFile, dtype=np.complex64).reshape(diffLength*2, diffWidth))[1:diffLength*2:2, :]
        diff1.astype(np.complex64).tofile(diffFile1)
        diff2.astype(np.complex64).tofile(diffFile2)
        create_xml(diffFile1, diffWidth, diffLength, 'int')
        create_xml(diffFile2, diffWidth, diffLength, 'int')

        #STEP 2. take looks
        if nrlks2 == 1 and nalks2 == 1:
            runCmd('cp {} {}'.format(diffFile1, diffLookFile1))
            runCmd('cp {} {}'.format(diffFile1+'.xml', diffLookFile1+'.xml'))
            runCmd('cp {} {}'.format(diffFile1+'.vrt', diffLookFile1+'.vrt'))
            runCmd('cp {} {}'.format(diffFile2, diffLookFile2))
            runCmd('cp {} {}'.format(diffFile2+'.xml', diffLookFile2+'.xml'))
            runCmd('cp {} {}'.format(diffFile2+'.vrt', diffLookFile2+'.vrt'))
        else:
            cmd = '$INSAR_ZERODOP_SCR/look.py -i {} -o {} -r {} -a {}'.format(diffFile1, diffLookFile1, nrlks2, nalks2)
            runCmd(cmd)
            cmd = '$INSAR_ZERODOP_SCR/look.py -i {} -o {} -r {} -a {}'.format(diffFile2, diffLookFile2, nrlks2, nalks2)
            runCmd(cmd)

        if inps.filt==1:
            tmpfile1='tmp1'
            tmpfile2='tmp2'
            cmd = "imageMath.py -e='a/(abs(a)+(abs(a)==0))' --a={} -o {} -t cfloat -s BIP".format(diffLookFile1, tmpfile1)
            runCmd(cmd)
            cmd = "imageMath.py -e='a/(abs(a)+(abs(a)==0))' --a={} -o {} -t cfloat -s BIP".format(diffLookFile2, tmpfile2)
            runCmd(cmd)
            cmd = "psfilt1 {} {} {} 3.0 64 1".format(tmpfile1, diffLookFile1, int(diffWidth/nrlks2))
            runCmd(cmd)
            cmd = "psfilt1 {} {} {} 3.0 64 1".format(tmpfile2, diffLookFile2, int(diffWidth/nrlks2))
            runCmd(cmd)
            os.remove(tmpfile1)
            os.remove(tmpfile2)
            os.remove(tmpfile1+'.xml')
            os.remove(tmpfile2+'.xml')
            os.remove(tmpfile1+'.vrt')
            os.remove(tmpfile2+'.vrt')

        #STEP 3. get difference
        diffLook1 = np.fromfile(diffLookFile1, dtype=np.complex64).reshape(int(diffLength/nalks2), int(diffWidth/nrlks2))
        diffLook2 = np.fromfile(diffLookFile2, dtype=np.complex64).reshape(int(diffLength/nalks2), int(diffWidth/nrlks2))
        diffLook = diffLook1 * np.conj(diffLook2)
        diffLook.astype(np.complex64).tofile(diffLookFile)
        create_xml(diffLookFile, int(diffWidth/nrlks2), int(diffLength/nalks2), 'int')

        #STEP 4. calculate mean angle
        cor = cal_coherence(diffLook)
        index = np.nonzero(np.logical_and(cor>corth, np.absolute(diffLook)!=0))
        angle = np.mean(np.angle(diffLook[index]), dtype=np.float64)
        diffMean.append(angle)
        print('phase offset: subswath{} - subswath{}: {}'.format(i, i+1, angle))

        #clear files
        os.remove(diffFile)
        os.remove(diffFile1)
        os.remove(diffFile1 + '.xml')
        os.remove(diffFile1 + '.vrt')
        os.remove(diffFile2)
        os.remove(diffFile2 + '.xml')
        os.remove(diffFile2 + '.vrt')
        os.remove(diffLookFile1)
        os.remove(diffLookFile1 + '.xml')
        os.remove(diffLookFile1 + '.vrt')
        os.remove(diffLookFile2)
        os.remove(diffLookFile2 + '.xml')
        os.remove(diffLookFile2 + '.vrt')

    if inps.phc == 1:
        #phase compensation for subswath interferograms (except the first subswath)
        cj = np.complex64(1j)
        for i in range(1, ns):
            print('correcting phase of subswath {}'.format(i+1))
            rinfc = np.fromfile(rinfs[i], dtype=np.complex64).reshape(rectLength[i], rectWidth[i])
            pha_offset = 1.0
            for j in range(1, i+1):
                pha_offset *= np.exp(cj * diffMean[j])
            rinfc *= pha_offset
            #overwrite original interferograms
            rinfc.astype(np.complex64).tofile(rinfs[i])

        #run again the command for mosaicking interferograms
        runCmd(cmdbak)


#######################################################################################
    ##tidy up: move important files to a backup dir
    procdir = 'subswath_mosaic_int_amp'
    os.mkdir(procdir)
    for i in range(1, ns):
        diffFile = '{}-{}_{}rlks_{}alks.int'.format(i-1, i, nrlks2, nalks2)
        shutil.move(diffFile, procdir)
        shutil.move(diffFile + '.xml', procdir)
        shutil.move(diffFile + '.vrt', procdir)
    
    #tidy up: delete less important files
    for i in range(1, ns):
        rectInputFileInf = 'rect_with_looks_inf_{}.in'.format(i)
        rectInputFileAmp = 'rect_with_looks_amp_{}.in'.format(i)
        os.remove(rectInputFileInf)
        os.remove(rectInputFileAmp)

    for rinf, ramp in zip(rinfs, ramps):
        os.remove(rinf)
        os.remove(ramp)

