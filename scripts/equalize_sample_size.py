#!/usr/bin/env python3

#program for preparing stripmap data for doing ScanSAR-stripmap interferometry
#Cunren Liang, JPL/Caltech

import os
import sys
import shutil
import pickle
import datetime
import argparse
import math
import numpy as np

import isce
import isceobj
from imageMath import IML
from isceobj.Constants import SPEED_OF_LIGHT

from crlpac import runCmd
from crlpac import create_xml
from crlpac import overlapFrequency


def cmdLineParse():
    '''
    Command line parser.
    '''
    parser = argparse.ArgumentParser( description='equalize sample size')
    parser.add_argument('-slc', dest='slc', type=str, required=True,
            help = 'slc file')
    parser.add_argument('-frame', dest='frame', type=str, required=True,
            help = 'frame pickle file')
    parser.add_argument('-frame2', dest='frame2', type=str, required=True,
            help = 'the other frame pickle file')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':

    inps = cmdLineParse()

    #frame
    with open(inps.frame, 'rb') as f:
        frame = pickle.load(f)
    #the other frame
    with open(inps.frame2, 'rb') as f:
        frame2 = pickle.load(f)


    if math.fabs(frame.instrument.rangeSamplingRate - frame2.instrument.rangeSamplingRate) < 1.0 and math.fabs(frame.instrument.PRF - frame2.instrument.PRF) < 1.0:
        print('no need to resample {}. nothing is done by this program.'.format(inps.slc))
    else:
        os.rename(inps.slc, 'tmp.slc')
        #os.remove(inps.slc + '.par.xml')
        os.remove(inps.slc + '.pck')
        os.remove(inps.slc + '.xml')
        os.remove(inps.slc + '.vrt')

        outWidth = round(frame.getNumberOfSamples() * (1.0/frame.instrument.rangeSamplingRate) / (1.0/frame2.instrument.rangeSamplingRate))
        outLength = round(frame.getNumberOfLines() * (1.0/frame.instrument.PRF) / (1.0/frame2.instrument.PRF))
        cmd = '$INSAR_ZERODOP_BIN/resamp {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}'.format(
            'tmp.slc',
            inps.slc,
            'fake', 
            'fake', 
            outWidth,
            outLength,
            frame.getNumberOfSamples(),
            frame.getNumberOfLines(),
            frame.instrument.PRF,
            frame.dopCoeff[0], frame.dopCoeff[1], frame.dopCoeff[2], frame.dopCoeff[3],
            0.0, (1.0/frame2.instrument.rangeSamplingRate) / (1.0/frame.instrument.rangeSamplingRate) - 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, (1.0/frame2.instrument.PRF) / (1.0/frame.instrument.PRF) - 1.0
            )
        runCmd(cmd)


        #index = np.array(range(frame.getNumberOfSamples()))
        index2  = np.array(range(outWidth))
        index = np.array(range(outWidth)) * (1.0/frame2.instrument.rangeSamplingRate) / (1.0/frame.instrument.rangeSamplingRate)

        dop = frame.dopCoeff[0] + frame.dopCoeff[1] * index + frame.dopCoeff[2] * index**2 + frame.dopCoeff[3] * index**3
        p = np.polyfit(index2, dop, 3) 
        dopCoeff=list(np.flipud(p))

        ka = frame.fmrateCoeff[0] + frame.fmrateCoeff[1] * index**1 + frame.fmrateCoeff[2] * index**2
        p = np.polyfit(index2, ka, 2) 
        fmrateCoeff=list(np.flipud(p))

        #set output frame
        outframe = frame
        outframe.setNumberOfSamples(outWidth)
        outframe.setNumberOfLines(outLength)
        outframe.setFarRange(frame.getStartingRange() + (outWidth - 1.0) *0.5 * SPEED_OF_LIGHT / frame2.rangeSamplingRate)
        outframe.setSensingStop(frame.getSensingStart() + datetime.timedelta(seconds=(outLength - 1.0) / frame2.instrument.PRF))
        outframe.setSensingMid(frame.getSensingStart() + datetime.timedelta(seconds=(outLength - 1.0) / frame2.instrument.PRF/2.0))

        outframe.dopCoeff = dopCoeff
        outframe.fmrateCoeff = fmrateCoeff

        outframe.getInstrument().setPulseRepetitionFrequency(frame2.instrument.PRF)
        outframe.getInstrument().setRangeSamplingRate(frame2.instrument.rangeSamplingRate)
        outframe.getInstrument().setRangePixelSize(frame2.getInstrument().getRangePixelSize())


        #outframe.getInstrument().setRadarWavelength(SPEED_OF_LIGHT/centerfreq)
        #outframe.getInstrument().setPulseLength(math.fabs(overlapbandwidth/frame.instrument.chirpSlope))


        #remove orignal slc, leaving only filtered slc
        os.remove('tmp.slc')

        #output new frame pickle file
        with open(inps.frame, 'wb') as f:
            pickle.dump(outframe, f)

        create_xml(inps.slc, outWidth, outLength, 'slc')

















