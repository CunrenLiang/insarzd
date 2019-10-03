#!/usr/bin/env python3

#program for removing non-overlap spectrum
#Cunren Liang, JPL/Caltech

import os
import sys
import shutil
import pickle
import datetime
import argparse
import numpy as np

import isce
import isceobj
from imageMath import IML

from crlpac import runCmd


def cmdLineParse():
    '''
    Command line parser.
    '''
    parser = argparse.ArgumentParser( description='remove unsynchronized signal')
    parser.add_argument('-mslc', dest='mslc', type=str, required=True,
            help = 'master slc file')
    parser.add_argument('-mframe', dest='mframe', type=str, required=True,
            help = 'master frame pickle file')
    parser.add_argument('-sslc', dest='sslc', type=str, required=True,
            help = 'slave slc file')
    parser.add_argument('-sframe', dest='sframe', type=str, required=True,
            help = 'slave frame pickle file')
    parser.add_argument('-syncf', dest='syncf', type=str, required=True,
            help = 'burst synchronization file produced by burst_sync.py')
    parser.add_argument('-syncf2', dest='syncf2', type=str, required=True,
            help = 'burst synchronization file produced by burst_sync.py for determing if remove unsynchronized signal')
    parser.add_argument('-synthr', dest='synthr', type=float, default=85.0,
            help = 'burst synchronization threshhold (percent) for removing unsynchronized signal. Default: 85.0')
    parser.add_argument('-mfilt', dest='mfilt', type=int, default=1,
            help = 'filter master? 0: No. 1: Yes (default)')
    parser.add_argument('-sfilt', dest='sfilt', type=int, default=1,
            help = 'filter slave? 0: No. 1: Yes (default)')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':

    inps = cmdLineParse()


    #master frame
    with open(inps.mframe, 'rb') as f:
        mframe = pickle.load(f)

    #slave frame
    with open(inps.sframe, 'rb') as f:
        sframe = pickle.load(f)

################################################################################################
    if not (hasattr(mframe, 'nbraw') and hasattr(mframe, 'ncraw') and hasattr(sframe, 'nbraw') and hasattr(sframe, 'ncraw')):
        print('nothing is done by this program.')
        sys.exit(0)
################################################################################################



################################################################################################
#check if need to do mbf
    with open(inps.syncf2) as f:
        lines = f.readlines()
    for linex in lines:
        if 'burst synchronization:' in linex:
            sync_per = float(linex.split(':')[1].strip().strip('%'))
    if sync_per >= inps.synthr:
        print('no need to remove unsynchronized signal. nothing is done by this program.')
        sys.exit(0) 
################################################################################################










    with open(inps.syncf) as f:
        lines = f.readlines()

    for linex in lines:
        if 'image pair range offset' in linex:
            rgoffset = float(linex.split(':')[1].strip())
        if 'image pair azimuth offset' in linex:
            azoffset = float(linex.split(':')[1].strip())
        if 'number of unsynchronized lines in a burst' in linex:
            unsynLines = float(linex.split(':')[1].strip())


    if inps.mfilt == 1:
        os.rename(inps.mslc, 'tmp.slc')

        #calculate slave doppler in master's coordinates
        sdopCoeff = sframe.dopCoeff
        index = np.array(range(mframe.getNumberOfSamples())) + rgoffset
        sdop_m = sdopCoeff[0] + sdopCoeff[1] * index + sdopCoeff[2] * index**2 + sdopCoeff[3] * index**3
        p = np.polyfit(index - rgoffset, sdop_m, 3)
        sdopCoeff_m=list(np.flipud(p))

        #filter image
        cmd = '$INSAR_ZERODOP_BIN/mbf {inputf} {outputf} {nrg} {prf} {prf_frac} {nb} {nbg} {nboff} {bsl} {kacoeff0} {kacoeff1} {kacoeff2} {dopcoeff1_0} {dopcoeff1_1} {dopcoeff1_2} {dopcoeff1_3} {dopcoeff2_0} {dopcoeff2_1} {dopcoeff2_2} {dopcoeff2_3}'.format(
            inputf = 'tmp.slc',
            outputf = inps.mslc,
            nrg = mframe.getNumberOfSamples(),
            prf = mframe.PRF,
            prf_frac = 1.0,
            nb = mframe.nbraw,
            nbg = mframe.ncraw - mframe.nbraw,
            nboff = unsynLines,
            bsl = mframe.burstStartLines[0],
            kacoeff0 = mframe.fmrateCoeff[0],
            kacoeff1 = mframe.fmrateCoeff[1],
            kacoeff2 = mframe.fmrateCoeff[2],
            dopcoeff1_0 = mframe.dopCoeff[0],
            dopcoeff1_1 = mframe.dopCoeff[1],
            dopcoeff1_2 = mframe.dopCoeff[2],
            dopcoeff1_3 = mframe.dopCoeff[3],
            dopcoeff2_0 = sdopCoeff_m[0],
            dopcoeff2_1 = sdopCoeff_m[1],
            dopcoeff2_2 = sdopCoeff_m[2],
            dopcoeff2_3 = sdopCoeff_m[3]
            )
        runCmd(cmd)

        os.remove('tmp.slc')

    if inps.sfilt == 1:
        os.rename(inps.sslc, 'tmp.slc')

        #calculate master doppler in slave's coordinates
        mdopCoeff = mframe.dopCoeff
        index = np.array(range(sframe.getNumberOfSamples())) - rgoffset
        mdop_s = mdopCoeff[0] + mdopCoeff[1] * index + mdopCoeff[2] * index**2 + mdopCoeff[3] * index**3
        p = np.polyfit(index + rgoffset, mdop_s, 3)
        mdopCoeff_s=list(np.flipud(p))

        #filter image
        cmd = '$INSAR_ZERODOP_BIN/mbf {inputf} {outputf} {nrg} {prf} {prf_frac} {nb} {nbg} {nboff} {bsl} {kacoeff0} {kacoeff1} {kacoeff2} {dopcoeff1_0} {dopcoeff1_1} {dopcoeff1_2} {dopcoeff1_3} {dopcoeff2_0} {dopcoeff2_1} {dopcoeff2_2} {dopcoeff2_3}'.format(
            inputf = 'tmp.slc',
            outputf = inps.sslc,
            nrg = sframe.getNumberOfSamples(),
            prf = sframe.PRF,
            prf_frac = 1.0,
            nb = sframe.nbraw,
            nbg = sframe.ncraw - sframe.nbraw,
            nboff = -unsynLines,
            bsl = sframe.burstStartLines[0],
            kacoeff0 = sframe.fmrateCoeff[0],
            kacoeff1 = sframe.fmrateCoeff[1],
            kacoeff2 = sframe.fmrateCoeff[2],
            dopcoeff1_0 = sframe.dopCoeff[0],
            dopcoeff1_1 = sframe.dopCoeff[1],
            dopcoeff1_2 = sframe.dopCoeff[2],
            dopcoeff1_3 = sframe.dopCoeff[3],
            dopcoeff2_0 = mdopCoeff_s[0],
            dopcoeff2_1 = mdopCoeff_s[1],
            dopcoeff2_2 = mdopCoeff_s[2],
            dopcoeff2_3 = mdopCoeff_s[3]
            )
        runCmd(cmd)

        os.remove('tmp.slc')





















