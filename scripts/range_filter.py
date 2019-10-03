#!/usr/bin/env python3

#Cunren Liang, JPL/Caltech
#2016

#program for range filtering


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
    parser = argparse.ArgumentParser( description='range filtering of an InSAR pair')
    parser.add_argument('-mslc', dest='mslc', type=str, required=True,
            help = 'master slc file')
    parser.add_argument('-mframe', dest='mframe', type=str, required=True,
            help = 'master frame pickle file')
    parser.add_argument('-sslc', dest='sslc', type=str, required=True,
            help = 'slave slc file')
    parser.add_argument('-sframe', dest='sframe', type=str, required=True,
            help = 'slave frame pickle file')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':

    inps = cmdLineParse()

    with open(inps.mframe, 'rb') as f:
        mframe = pickle.load(f)
    with open(inps.sframe, 'rb') as f:
        sframe = pickle.load(f)

    centerfreq1 = SPEED_OF_LIGHT / mframe.radarWavelegth
    bandwidth1 = math.fabs(mframe.instrument.pulseLength * mframe.instrument.chirpSlope)
    centerfreq2 = SPEED_OF_LIGHT / sframe.radarWavelegth
    bandwidth2 = math.fabs(sframe.instrument.pulseLength * sframe.instrument.chirpSlope)
    overlapfreq = overlapFrequency(centerfreq1, bandwidth1, centerfreq2, bandwidth2)

    if overlapfreq == None:
        raise Exception('there is no overlap bandwidth in range. nothing is done by this program.\n')
    overlapbandwidth = overlapfreq[1] - overlapfreq[0]
    if overlapbandwidth < 3e6:
        print('overlap bandwidth: {}, percentage: {}%'.format(overlapbandwidth, 100.0*overlapbandwidth/bandwidth1))
        raise Exception('there is not enough overlap bandwidth. nothing is done by this program.\n')
    centerfreq = (overlapfreq[1] + overlapfreq[0]) / 2.0

    #if math.fabs(centerfreq1 - centerfreq2) > 1.0:
    #    raise Exception('not good to do interferometry with different wavelengths for ALOS-2.\n')


    #1. filtering master
    if math.fabs(centerfreq1 - centerfreq) < 1.0 and (bandwidth1 - 1.0) < overlapbandwidth:
        print('no need to range filter {}. nothing is done by this program.'.format(inps.mslc))
    else:
        os.rename(inps.mslc, 'tmp.slc')
        #os.remove(inps.mslc + '.par.xml')
        os.remove(inps.mslc + '.pck')
        os.remove(inps.mslc + '.xml')
        os.remove(inps.mslc + '.vrt')

        cmd = '$INSAR_ZERODOP_BIN/rg_filter {inputf} {outputf} {nrg} {bw} {bc} {nfilter} {nfft} {beta}'.format(
            inputf = 'tmp.slc',
            outputf = inps.mslc,
            nrg = mframe.getNumberOfSamples(),
            bw = overlapbandwidth / mframe.instrument.rangeSamplingRate,
            bc = (centerfreq - centerfreq1) / mframe.instrument.rangeSamplingRate,
            nfilter = 129,
            nfft = 512,
            beta = 0.1
            )
        runCmd(cmd)

        #remove orignal slc, leaving only filtered slc
        os.remove('tmp.slc')

        #output new frame pickle file
        outframe = mframe
        outframe.getInstrument().setRadarWavelength(SPEED_OF_LIGHT/centerfreq)
        outframe.getInstrument().setPulseLength(math.fabs(overlapbandwidth/mframe.instrument.chirpSlope))
        with open(inps.mframe, 'wb') as f:
            pickle.dump(outframe, f)

        create_xml(inps.mslc, mframe.getNumberOfSamples(), mframe.getNumberOfLines(), 'slc')


    #2. filtering slave
    if math.fabs(centerfreq2 - centerfreq) < 1.0 and (bandwidth2 - 1.0) < overlapbandwidth:
        print('no need to range filter {}. nothing is done by this program.'.format(inps.sslc))
    else:
        os.rename(inps.sslc, 'tmp.slc')
        #os.remove(inps.sslc + '.par.xml')
        os.remove(inps.sslc + '.pck')
        os.remove(inps.sslc + '.xml')
        os.remove(inps.sslc + '.vrt')

        cmd = '$INSAR_ZERODOP_BIN/rg_filter {inputf} {outputf} {nrg} {bw} {bc} {nfilter} {nfft} {beta}'.format(
            inputf = 'tmp.slc',
            outputf = inps.sslc,
            nrg = sframe.getNumberOfSamples(),
            bw = overlapbandwidth / sframe.instrument.rangeSamplingRate,
            bc = (centerfreq - centerfreq2) / sframe.instrument.rangeSamplingRate,
            nfilter = 129,
            nfft = 512,
            beta = 0.1
            )
        runCmd(cmd)

        #remove orignal slc, leaving only filtered slc
        os.remove('tmp.slc')

        #output new frame pickle file
        outframe = sframe
        outframe.getInstrument().setRadarWavelength(SPEED_OF_LIGHT/centerfreq)
        outframe.getInstrument().setPulseLength(math.fabs(overlapbandwidth/sframe.instrument.chirpSlope))
        with open(inps.sframe, 'wb') as f:
            pickle.dump(outframe, f)

        create_xml(inps.sslc, sframe.getNumberOfSamples(), sframe.getNumberOfLines(), 'slc')

