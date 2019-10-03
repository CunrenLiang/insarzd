#!/usr/bin/env python3

import os
import sys
import argparse
import pickle

import isce
import isceobj
from isceobj.Constants import SPEED_OF_LIGHT

from crlpac import getWidth
from crlpac import getLength
from crlpac import runCmd
from crlpac import create_xml


def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser( description='interferometry')
    parser.add_argument('-m', '--mframe', dest='mframe', type=str, required=True,
            help = 'master frame (picke file)')
    parser.add_argument('-i', '--intf', dest='intf', type=str, required=True,
            help = 'interferogram')
    parser.add_argument('-r','--rgoff', dest='rgoff', type=str, required=True,
            help = 'range offset file')
    parser.add_argument('-d', '--dintf', dest='dintf', type=str, required=True,
            help = '(output) differential inteferogram')
    parser.add_argument('-rlks', dest='rlks', type=int, default=1,
            help = 'number of range looks')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':
    
    inps = cmdLineParse()

    with open(inps.mframe, 'rb') as fid:
        masterFrame = pickle.load(fid)

    slantRangePixelSpacing = inps.rlks * 0.5 * SPEED_OF_LIGHT / masterFrame.rangeSamplingRate
    radarWavelength = masterFrame.radarWavelegth

    intfWidth = getWidth(inps.intf + '.xml')
    intfLength = getLength(inps.intf + '.xml')

    #run it
    cmd = "$INSAR_ZERODOP_BIN/flat {} {} {} {} {}".format(inps.intf, inps.rgoff, inps.dintf, slantRangePixelSpacing, radarWavelength)
    #print("{}".format(cmd))
    runCmd(cmd)

    #get xml file for interferogram
    create_xml(inps.dintf, intfWidth, intfLength, 'int')

    #./flat.py -m 20130927.slc.pck -i 20130927-20141211.int -r 20130927-20141211_rg.off -d diff_20130927-20141211.int





