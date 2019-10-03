#!/usr/bin/env python3

#Cunren Liang, JPL/Caltech

import os
import re
import sys
import glob
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


def cmdLineParse():
    '''
    Command line parser.
    '''
    parser = argparse.ArgumentParser( description='range filter image')
    parser.add_argument('-islc', dest='islc', type=str, required=True,
            help = 'input slc')
    parser.add_argument('-iframe', dest='iframe', type=str, required=True,
            help = 'input frame pickle file')
    parser.add_argument('-oslc', dest='oslc', type=str, required=True,
            help = 'output slc')
    parser.add_argument('-oframe', dest='oframe', type=str, default='NotSpecified',
            help = 'output frame pickle file')
    parser.add_argument('-lu', dest='lu', type=int, default=0,
            help = 'lower or upper band? 0: lower (Default). 1: upper')
    parser.add_argument('-offset', dest='offset', type=str, default='NotSpecified',
            help = 'offset file containing the subswath offset')
    parser.add_argument('-sn', dest='sn', type=int, default=1,
            help = 'subswath (starting from 0). only works when offset file provided. default: 0')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':

    inps = cmdLineParse()

    with open(inps.iframe, 'rb') as f:
        iframe = pickle.load(f)

    rgbandwidth = abs(iframe.instrument.pulseLength * iframe.instrument.chirpSlope)
    if inps.lu == 0:
        bc0 = -rgbandwidth / 3.0 / iframe.instrument.rangeSamplingRate
        fx = SPEED_OF_LIGHT / iframe.radarWavelegth - rgbandwidth / 3.0
    else:
        bc0 = rgbandwidth / 3.0 / iframe.instrument.rangeSamplingRate
        fx = SPEED_OF_LIGHT / iframe.radarWavelegth + rgbandwidth / 3.0

#read offset from subswath offset file
###########################################################################################
    if inps.offset != 'NotSpecified':

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
        ns = len(rangeOffset1)

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

        range_offset = -1.0 * rangeOffset[inps.sn]
    else:
        range_offset = 0.0


    cmd = '$INSAR_ZERODOP_BIN/rg_filter {inputf} {outputf} {nrg} {bw} {bc} {nfilter} {nfft} {beta} {zero_cf} {offset}> /dev/null'.format(
    inputf = inps.islc,
    outputf = inps.oslc,
    nrg = iframe.getNumberOfSamples(),
    bw = rgbandwidth / 3.0 / iframe.instrument.rangeSamplingRate,
    bc = bc0,
    nfilter = 129,
    nfft = 512,
    beta = 0.1,
    zero_cf = 0,
    offset = range_offset
    )
    runCmd(cmd)
    create_xml(inps.oslc, iframe.getNumberOfSamples(), iframe.getNumberOfLines(), 'slc')


    iframe.instrument.setPulseLength(iframe.instrument.pulseLength/3.0)
    iframe.getInstrument().setRadarWavelength(float(SPEED_OF_LIGHT/fx))
    if inps.oframe != 'NotSpecified':
        with open(inps.oframe, 'wb') as f:
            pickle.dump(iframe, f)




