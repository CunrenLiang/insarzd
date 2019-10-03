#!/usr/bin/env python3

#Cunren Liang, JPL/Caltech, 2017

import os
import sys
import argparse
import numpy as np

from crlpac import getWidth
from crlpac import getLength
from crlpac import create_xml


def cmdLineParse():
    '''
    Command line parser.
    '''
    parser = argparse.ArgumentParser( description='scale amplitude and interferogram')
    parser.add_argument('-inf', dest='inf', type=str, required=True,
            help = 'interferogram file')
    parser.add_argument('-amp', dest='amp', type=str, required=True,
            help = 'ampitude file')
    parser.add_argument('-ratio', dest='ratio', type=float, default=1.0e5,
            help='scaling ratio. default: 1.0e5')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':

    inps = cmdLineParse()

    #ratio = 1.0e10
    ratio = inps.ratio
    print('master and slave scaling ratio: {}'.format(ratio))

    width = getWidth(inps.inf + '.xml')
    length = getLength(inps.inf + '.xml')

    #read master burst
    inf = np.fromfile(inps.inf, dtype=np.complex64).reshape(length, width)
    amp = np.fromfile(inps.amp, dtype=np.complex64).reshape(length, width)

    flag = (inf!=0)*(amp.real!=0)*(amp.imag!=0)
    nvalid = np.sum(flag, dtype=np.float64)

    mpwr1 =  np.sqrt(np.sum(amp.real * amp.real * flag, dtype=np.float64) / nvalid)
    mpwr2 =  np.sqrt(np.sum(amp.imag * amp.imag * flag, dtype=np.float64) / nvalid)

    amp.real = amp.real / ratio
    amp.imag = amp.imag / ratio * mpwr1 / mpwr2
    inf = inf / ratio / ratio * mpwr1 / mpwr2

    amp.astype(np.complex64).tofile(inps.amp)
    inf.astype(np.complex64).tofile(inps.inf)

#according to equation (2)
#Howard A. Zebker and Katherine Chen, Accurate Estimation of Correlation in InSAR Observations
#IEEE GEOSCIENCE AND REMOTE SENSING LETTERS, VOL. 2, NO. 2, APRIL 2005
#
#the operation of the program does not affect coherence estimation



