#!/usr/bin/env python3

#Cunren Liang, 22-MAY-2015
#JPL/Caltech

############################################################################
#program for mosaicking slcs/interferograms/amplitude images of consecutive frames
#
#convention: paramters of frame pickle file are single-look parameters
#            all the range/azimuth number of looks of the frame images are the same
#
#10-apr-2017, crl, support burst-by-burst processing
#24-apr-2017, CL. support mosaicking multiple-look images 
############################################################################


import os
import re
import sys
import argparse
import datetime
import pickle
import ntpath
import shutil
import numpy as np
from xml.etree.ElementTree import ElementTree

import isce
import isceobj
import mroipac
from mroipac.ampcor.Ampcor import Ampcor
from isceobj.Orbit.Orbit import Orbit
from iscesys.StdOEL.StdOELPy import create_writer
from isceobj.Constants import SPEED_OF_LIGHT

from crlpac import create_xml
from crlpac import getWidth
from crlpac import getLength
from crlpac import runCmd
from crlpac import writeOffset
from crlpac import meanOffset
from crlpac import cullOffset
from crlpac import getOffset


def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="program for mosaicking slcs/interferograms/amplitude images of consecutive frames\
                    \n\
                    ")
    parser.add_argument('-f', action='append', dest='frame', default=[],
            help = 'frame pickle files IN THE ORDER OF acquisition time')
    parser.add_argument('-i', action='append', dest='inf', default=[],
            help = 'frame slcs/interferograms/amplitude images IN THE ORDER OF acquisition time')
    parser.add_argument('-offset', dest='offset', type=str, required=True,
            help = 'frame offset file')
    parser.add_argument('-io', dest='infoutput', type=str, required=True,
            help = 'mosaicked slcs/interferograms/amplitude image')
    parser.add_argument('-fo', dest='frameoutput', type=str, required=True,
            help = 'output frame pickle file for mosaicked image. also generate *.par.xml file based on this')
    parser.add_argument('-rlks', dest='rlks', type=int, default=1,
            help = 'number of range looks of the slcs/interferograms/amplitude images. default: 1')
    parser.add_argument('-alks', dest='alks', type=int, default=1,
            help = 'number of azimuth looks of the slcs/interferograms/amplitude images. default: 1')
    parser.add_argument('-phc', dest='phc', type=int, default=0,
            help = 'whether do phase correction to frame file (except the first frame). 0: No (default). 1: Yes.')

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
    rlks = inps.rlks
    alks = inps.alks

    if not (len(frames) == len(infs)):
        raise Exception('the numbers of files are not equal!')
    if len(frames) < 2:
        raise Exception('at least 2 frames should be specified!')

    nframe = len(frames)

    if os.path.splitext(ntpath.basename(infs[0]))[1] == '.int':
        file_type = 'int'
    elif os.path.splitext(ntpath.basename(infs[0]))[1] == '.amp':
        file_type = 'amp'
    elif os.path.splitext(ntpath.basename(infs[0]))[1] == '.slc':
        file_type = 'slc'
    else:
        raise Exception('file type not supported yet!')

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
        with open(frame, 'rb') as f:
            pcks.append(pickle.load(f))

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


    #read offset
    #################################
    #get offset from gemoetrical parameters
    rangeOffset1 = []
    azimuthOffset1 = []
    #get offset from matching
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


    rangeOffset = []
    azimuthOffset = []
    rangeOffset.append(0.0)
    azimuthOffset.append(0.0)
    #get offset relative to the first frame
    for i in range(1, nframe):
        rangeOffset.append(0.0)
        azimuthOffset.append(0.0)
        for j in range(1, i+1):
            rangeOffset[i] += rangeOffset2[j]
            azimuthOffset[i] += azimuthOffset2[j]


    #resample each frame
    rinfs = []
    crinfs = []
    for i, inf in enumerate(infs):
        infBase = os.path.splitext(ntpath.basename(inf))[0]
        infExt = os.path.splitext(ntpath.basename(inf))[1]
        rinfs.append("{}_{}{}".format(infBase, i, infExt))
        crinfs.append("{}_{}_corrected{}".format(infBase, i, infExt))

        if i == 0:
            os.symlink(inf, rinfs[i])
            #create xml file for rectified interferogram
            create_xml(rinfs[i], int(width[i]/rlks), int(length[i]/alks), file_type)
            if inps.phc == 1:
                os.symlink(inf, crinfs[i])
                #create xml file for corrected interferogram
                create_xml(crinfs[i], int(width[i]/rlks), int(length[i]/alks), file_type)
            continue

        rangeOffsetFrac = rangeOffset[i] - int(rangeOffset[i])
        azimuthOffsetFrac = azimuthOffset[i] - int(azimuthOffset[i])

        if file_type == 'int' or file_type == 'amp':
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
            inWidth = int(width[i]/rlks), inLength = int(length[i]/alks),
            outWidth = int(width[i]/rlks), outLength = int(length[i]/alks),
            m11 = 1.0000000000, m12 = 0.0000000000,
            m21 = 0.0000000000, m22 = 1.0000000000,
            t1 = rangeOffsetFrac, t2 = azimuthOffsetFrac,
            format = 'COMPLEX',
            method = 'Sinc'
            )
            rectInputFileInf = 'rect_with_looks_inf_{}.in'.format(i)
            with open(rectInputFileInf, 'w') as f:
                f.write(rectInputInf)
            cmd = '{}/rect_with_looks {}'.format(BIN_DIR, rectInputFileInf)
        elif file_type == 'slc':
            if rlks != 1 or alks != 1:
                print('***WARNING: you are resampling slc, but it is not single look!\n\n')
            #change the location of resamp_mosaic here!!!!!!!!!!!!
            cmd = "{}/resamp {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}".format(BIN_DIR,
                inf, 
                rinfs[i], 
                'fake', 
                'fake', 
                int(width[i]/rlks),
                int(length[i]/alks),
                int(width[i]/rlks),
                int(length[i]/alks),
                prf[i],
                dopCoeff[i][0], dopCoeff[i][1], dopCoeff[i][2], dopCoeff[i][3],
                rangeOffsetFrac, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                azimuthOffsetFrac
                )
        runCmd(cmd)
        #create xml file for rectified interferogram
        create_xml(rinfs[i], int(width[i]/rlks), int(length[i]/alks), file_type)

        if inps.phc == 1:
            #correct phase
            cmd = '{BIN_DIR}/correctphase {input0} {nrgin0} {input1} {nrgin1} {nrgoff} {nazoff} {output}'.format(BIN_DIR=BIN_DIR,
                input0 = rinfs[i-1],
                nrgin0 = int(width[i-1]/rlks),
                input1 = rinfs[i],
                nrgin1 = int(width[i]/rlks),
                nrgoff = int(rangeOffset2[i]),
                nazoff = int(azimuthOffset2[i]),
                output = crinfs[i]
                )
            runCmd(cmd)
            #create xml file for corrected interferogram
            create_xml(crinfs[i], int(width[i]/rlks), int(length[i]/alks), file_type)


    xs = []
    xe = []
    ys = []
    ye = []
    #determine output width and length
    #actually no need to calculate in azimuth direction
    for i in range(nframe):
        if i == 0:
            xs.append(0)
            xe.append(int(width[i]/rlks) - 1)
            ys.append(0)
            ye.append(int(length[i]/alks) - 1)
        else:
            xs.append(0 - int(rangeOffset[i]))
            xe.append(int(width[i]/rlks) - 1 - int(rangeOffset[i]))
            ys.append(0 - int(azimuthOffset[i]))
            ye.append(int(length[i]/alks) - 1 - int(azimuthOffset[i]))

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

    #mosaic interferogram later, so that the interferogram difference can override amplitude difference file
    infOutput = inps.infoutput
    cmd = "{}/mosaicframe {} {} {} {}".format(BIN_DIR, infOutput, outWidth, outLength, nframe)
    for i in range(nframe):
        if inps.phc == 1:
            input_file = crinfs[i]
        else:
            input_file = rinfs[i]
        #mosaic at the center of the overlap area
        cmdi = " {} {} {} {} {}".format(input_file, int(width[i]/rlks), rangeOffset3[i], azimuthOffset3[i], 0)
        cmd += cmdi
    runCmd(cmd)
    #create xml file for output image
    create_xml(infOutput, outWidth, outLength, file_type)

    #create/update auxiliary files
    ############################################################################
    #1. update pickle file
    outDopCoeff = [0, 0, 0, 0]
    for dopCoeffx in dopCoeff:
        for i in range(4):
            outDopCoeff[i] += dopCoeffx[i] / nframe

    #assume there are outWidth*rlks samples, and outLength*alks lines
    #outStartingRange = startingRange[0] - rangeOffset3[0] * rangePixelSize[0]
    #fixed the uppler line's bug, CL, 03-JUN-2015
    outStartingRange = startingRange[0] + rangeOffset3[0]*rlks * rangePixelSize[0]
    outFarRange = outStartingRange + (outWidth*rlks - 1) * rangePixelSize[0]
    outSensingStart = sensingStart[0]
    outSensingStop = outSensingStart + datetime.timedelta(seconds=(outLength*alks - 1) / prf[0])

    outpck = pcks[0]
    outpck.dopCoeff = outDopCoeff
    outpck.setNumberOfSamples(outWidth*rlks)
    outpck.setNumberOfLines(outLength*alks)
    outpck.startingRange = outStartingRange
    outpck.setFarRange(outFarRange)
    outpck.setSensingStart(outSensingStart)
    outpck.setSensingStop(outSensingStop)
    ##still needs to determin which orbit to use:
    #outpck.setOrbit(orbit[nframe - 1])
    #sensing mid
    outpck.setSensingMid(outpck.getSensingStart() + datetime.timedelta(seconds=(outpck.getSensingStop()-outpck.getSensingStart()).total_seconds()/2.0))


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

    ############################################################################
    #with open('./alos2/new/150221-150502/f560/150221.slc.pck', 'rb') as fid:
    #    tmppck = pickle.load(fid)
    #    outpck.setOrbit(tmppck.getOrbit())
    #with open('./alos2/new/150221-150502/f560/150502.slc.pck', 'rb') as fid:
    #    tmppck = pickle.load(fid)
    #    outpck.setOrbit(tmppck.getOrbit())
    ############################################################################

    with open(inps.frameoutput, 'wb') as f:
        pickle.dump(outpck, f)


    # tidy up
    ############################################################################
    #1. move important files to a backup dir
    if file_type == 'slc':
        output = os.path.splitext(ntpath.basename(frames[0]))[0]
        procdir = 'frame_mosaic_' + os.path.splitext(output)[0]
    else:
        procdir = 'frame_mosaic_' + file_type
    os.mkdir(procdir)

    for i in range(nframe - 1):
        diffFile = '{}-{}.int'.format(i, i+1)
        diffLength = os.stat(diffFile).st_size / 8 / outWidth
        create_xml(diffFile, outWidth, diffLength, 'int')
        shutil.move(diffFile, procdir)
        shutil.move(diffFile + '.xml', procdir)
        shutil.move(diffFile + '.vrt', procdir)


    #2. delete less important files
    for rinf, crinf in zip(rinfs, crinfs):
        os.remove(rinf)
        os.remove(rinf + '.xml')
        os.remove(rinf + '.vrt')
        if inps.phc == 1:
            os.remove(crinf)
            os.remove(crinf + '.xml')
            os.remove(crinf + '.vrt')

    if file_type == 'int' or file_type == 'amp':
        for i in range(1, nframe):
            rectInputFileInf = 'rect_with_looks_inf_{}.in'.format(i)
            os.remove(rectInputFileInf)


    tobedone = '\n\n'
    tobedone += 'Mosaicking frame completed. Please copy the following\n'
    tobedone += 'configure files from frame 1 directory to mosaicked frame directory.\n'
    tobedone += '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'
    tobedone += 'master.xml\n'
    tobedone += 'slave.xml\n'
    tobedone += 'insarzd.xml\n\n'
    tobedone += 'And mosaic pickle file for slave using frame_frame.py\n'
    tobedone += 're-calculate baseline using $INSAR_ZERODOP_SCR/calBaseline.py\n\n'

    print(tobedone)



