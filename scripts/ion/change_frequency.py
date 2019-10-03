#!/usr/bin/env python3

#program for changing range frequency
#Cunren Liang, 2016
#JPL/Caltech


import os
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


def cmdLineParse():
    '''
    Command line parser.
    '''
    parser = argparse.ArgumentParser( description='range filter bursts')
    parser.add_argument('-iframe', dest='iframe', type=str, required=True,
            help = 'input frame pickle file')
    parser.add_argument('-oframe', dest='oframe', type=str, default='NotSpecified',
            help = 'output frame pickle file')
    parser.add_argument('-lu', dest='lu', type=int, default=0,
            help = 'lower or upper band? 0: lower (Default). 1: upper')

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

    iframe.instrument.setPulseLength(iframe.instrument.pulseLength/3.0)
    iframe.getInstrument().setRadarWavelength(float(SPEED_OF_LIGHT/fx))
    if inps.oframe != 'NotSpecified':
        with open(inps.oframe, 'wb') as f:
            pickle.dump(iframe, f)



            


