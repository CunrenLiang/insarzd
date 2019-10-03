#!/usr/bin/env python3

#program for fixing relative phase unwrapping error
#Cunren Liang, 2016
#JPL/Caltech


import os
import sys
import glob
import shutil
import pickle
import datetime
import argparse
import numpy as np
import numpy.matlib
import scipy.signal as ss
from xml.etree.ElementTree import ElementTree

import isce
import isceobj
from imageMath import IML

from crlpac import getWidth
from crlpac import getLength
from crlpac import create_xml


def cmdLineParse():
    '''
    Command line parser.
    '''
    parser = argparse.ArgumentParser( description='fix relative unwrapping error of upper band interferogram')

    parser.add_argument('-unwl', dest='unwl', type=str, required=True,
            help = 'unwrapped interferogram of lower band (without mask)')
    parser.add_argument('-ccl', dest='ccl', type=str, required=True,
            help = 'connected components of unwrapped interferogram of lower band')
    parser.add_argument('-unwu', dest='unwu', type=str, required=True,
            help = 'unwrapped interferogram of upper band (without mask)')
    parser.add_argument('-ccu', dest='ccu', type=str, required=True,
            help = 'connected components of unwrapped interferogram of upper band')
    parser.add_argument('-cor', dest='cor', type=str, required=True,
            help = 'coherence of lower band')
    parser.add_argument('-unwc', dest='unwc', type=str, required=True,
            help = 'output corrected unwrapped interferogram of upper band')
    parser.add_argument('-cot', dest='cot', type=float, default=0.1,
            help = 'coherence threshhold for phase difference calculation. Default: 0.1')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':

    inps = cmdLineParse()

    width = getWidth(inps.unwl + '.xml')
    length = getLength(inps.unwl + '.xml')

    unwl = (np.fromfile(inps.unwl, dtype=np.float32).reshape(length*2, width))[1:length*2:2, :]
    unwu = (np.fromfile(inps.unwu, dtype=np.float32).reshape(length*2, width))[1:length*2:2, :]
    cor   = (np.fromfile(inps.cor,  dtype=np.float32).reshape(length*2, width))[1:length*2:2, :]
    
    ccl = np.fromfile(inps.ccl, dtype=np.int8).reshape(length, width)
    ccu = np.fromfile(inps.ccu, dtype=np.int8).reshape(length, width)

    flag = (unwl!=0)*(unwu!=0)*(cor>=inps.cot)*(ccl>0)*(ccu>0)
    index = np.nonzero(flag!=0)

    #get mean value
    mv = np.mean((unwl - unwu)[index], dtype=np.float64)
    print('mean value of phase difference: {}'.format(mv))


    flag2 = (unwl!=0)*(unwu!=0)
    index2 = np.nonzero(flag2)

    #phase for adjustment
    unwd = ((unwl - unwu)[index2] - mv) / (2.0*np.pi)
    unw_adj = np.around(unwd) * (2.0*np.pi)

    #ajust phase of upper band
    unwu[index2] += unw_adj

    unw_diff = (unwl - unwu)
    print('after adjustment:')
    print('max phase difference: {}'.format(np.amax(unw_diff)))
    print('min phase difference: {}'.format(np.amin(unw_diff)))




    #get two band rmg file of ajusted upper band
    unwc = np.zeros((length*2, width))
    unwc[0:length*2:2, :] = (np.fromfile(inps.unwu, dtype=np.float32).reshape(length*2, width))[0:length*2:2, :]
    unwc[1:length*2:2, :] = unwu

    #output
    unwc.astype(np.float32).tofile(inps.unwc)

    #create xml
    create_xml(inps.unwc, width, length, 'unw')













