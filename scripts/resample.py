#!/usr/bin/env python3

import os
import sys
import shutil
import argparse
import pickle

import isce
import isceobj
from imageMath import IML

from crlpac import getWidth
from crlpac import getLength
from crlpac import runCmd


def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser( description='resample slave image')
    parser.add_argument('-m', '--master', dest='master', type=str, required=True,
            help = 'master SLC')
    parser.add_argument('-s', '--slave', dest='slave', type=str, required=True,
            help = 'slave SLC')
    parser.add_argument('-r', '--rgoff', dest='rgoff', type=str, required=True,
            help = 'range offset file')
    parser.add_argument('-a', '--azoff', dest='azoff', type=str, required=True,
            help = 'azimuth offset file')
    parser.add_argument('-o', '--rslave', dest='rslave', type=str, required=True,
            help = '(output) resampled slave image')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':
    
    inps = cmdLineParse()

    #get information
    masterWidth = getWidth(inps.master + '.xml')
    masterLength = getLength(inps.master + '.xml')
    slaveWidth = getWidth(inps.slave + '.xml')
    slaveLength = getLength(inps.slave + '.xml')

    with open(inps.slave + '.pck', 'rb') as fid:
        slaveFrame = pickle.load(fid)
    prf = slaveFrame.PRF

    #can get doppler centroid frequency information from frame now
    dop0 = slaveFrame.dopCoeff[0];
    dop1 = slaveFrame.dopCoeff[1];
    dop2 = slaveFrame.dopCoeff[2];
    dop3 = slaveFrame.dopCoeff[3];

    #get xml file for resampled SLC
    shutil.copy(inps.master + '.xml', inps.rslave + '.xml')
    #fix file name
    img, dataName, metaName = IML.loadImage(inps.rslave + '.xml')
    img.filename = inps.rslave
    #img.setAccessMode('READ')
    img.renderHdr()



    #run resamp
    cmd = "$INSAR_ZERODOP_BIN/resamp {} {} {} {} {} {} {} {} {} {} {} {} {}".format(inps.slave, inps.rslave, inps.rgoff, inps.azoff, masterWidth, masterLength, slaveWidth, slaveLength, prf, dop0, dop1, dop2, dop3)
    #print("{}".format(cmd))
    runCmd(cmd)

#./resample.py -m 20130927.slc -s 20141211.slc -r 20130927-20141211_rg.off -a 20130927-20141211_az.off -o 20141211_resamp.slc



