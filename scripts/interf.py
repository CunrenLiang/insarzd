#!/usr/bin/env python3

import os
import sys
import shutil
import argparse

import isce
import isceobj

from crlpac import getWidth
from crlpac import getLength
from crlpac import runCmd
from crlpac import create_xml


def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser( description='interferometry')
    parser.add_argument('-m', '--master', dest='master', type=str, required=True,
            help = 'master SLC')
    parser.add_argument('-s', '--slave', dest='slave', type=str, required=True,
            help = 'resampled slave SLC')
    parser.add_argument('-i', '--intf', dest='intf', type=str, required=True,
            help = '(output) interferogram')
    parser.add_argument('-a', '--amp', dest='amp', type=str, required=True,
            help = '(output) amplitudes of master and slave SLCs')

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

    #run interf
    cmd = "$INSAR_ZERODOP_BIN/interf {} {} {} {} {}".format(inps.master, inps.slave, inps.intf, inps.amp, masterWidth)
    #print("{}".format(cmd))
    runCmd(cmd)

    #get xml file for interferogram
    create_xml(inps.intf, masterWidth, masterLength, 'int')
    #get xml file for amplitude image
    create_xml(inps.amp, masterWidth, masterLength, 'amp')


    #./interf.py -m 20130927.slc -s 20141211.slc -i 20130927-20141211.int -a 20130927-20141211.amp




