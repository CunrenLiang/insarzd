#!/usr/bin/env python3

#Cunren Liang, JPL/Caltech

import os
import sys
import argparse

import isce
import isceobj

from crlpac import getWidth
from crlpac import getLength
from crlpac import runCmd
from crlpac import create_xml


def ampLooks(inps):


    inWidth = getWidth(inps.input + '.xml')
    inLength = getLength(inps.input + '.xml')
    outWidth = int(inWidth/inps.rlks)
    outLength = int(inLength/inps.alks)

    #run it
    #cmd = 'echo -e "{}\n{}\n{} {}\n{} {}\n" | $INSAR_ZERODOP_BIN/rilooks'.format(inps.input, inps.output, inWidth, inLength, inps.rlks, inps.alks)
    #it seems that echo does not require -e in this situation, strange
    cmd = '$INSAR_ZERODOP_BIN/look {} {} {} {} {} 4 1 0'.format(inps.input, inps.output, inWidth, inps.rlks, inps.alks)
    runCmd(cmd)

    #get xml file for amplitude image
    create_xml(inps.output, outWidth, outLength, 'amp')


def intLooks(inps):

    inWidth = getWidth(inps.input + '.xml')
    inLength = getLength(inps.input + '.xml')
    outWidth = int(inWidth/inps.rlks)
    outLength = int(inLength/inps.alks)

    #run program here
    cmd = '$INSAR_ZERODOP_BIN/look {} {} {} {} {} 4 0 0'.format(inps.input, inps.output, inWidth, inps.rlks, inps.alks)
    runCmd(cmd)
    
    #get xml file for interferogram
    create_xml(inps.output, outWidth, outLength, 'int')


def mskLooks(inps):

    inWidth = getWidth(inps.input + '.xml')
    inLength = getLength(inps.input + '.xml')
    outWidth = int(inWidth/inps.rlks)
    outLength = int(inLength/inps.alks)

    #look_msk infile outfile nrg nrlks nalks
    #run program here
    cmd = '$INSAR_ZERODOP_BIN/look_msk {} {} {} {} {}'.format(inps.input, inps.output, inWidth, inps.rlks, inps.alks)
    runCmd(cmd)
    
    #get xml file for interferogram
    image = isceobj.createImage()
    accessMode = 'read'
    dataType = 'BYTE'
    bands = 1
    scheme = 'BIL'
    width = outWidth
    image.initImage(inps.output, accessMode, width, dataType, bands=bands, scheme=scheme)
    descr = 'Radar shadow-layover mask. 1 - Radar Shadow. 2 - Radar Layover. 3 - Both.'
    image.addDescription(descr)
    image.renderHdr()
    #image.finalizeImage()


def hgtLooks(inps):

    inWidth = getWidth(inps.input + '.xml')
    inLength = getLength(inps.input + '.xml')
    outWidth = int(inWidth/inps.rlks)
    outLength = int(inLength/inps.alks)

    #look_msk infile outfile nrg nrlks nalks
    #run program here
    cmd = '$INSAR_ZERODOP_BIN/look {} {} {} {} {} 3 0 1'.format(inps.input, inps.output, inWidth, inps.rlks, inps.alks)
    runCmd(cmd)
    

    #get xml
    image = isceobj.createImage()
    accessMode = 'read'
    dataType = 'DOUBLE'
    width = outWidth
    image.initImage(inps.output,accessMode,width,dataType)

    image.addDescription('Pixel-by-pixel height in meters.')
    image.renderHdr()
    #image.finalizeImage()
    #image.createImage()



def cmdLineParse():
    '''
    Command line parser.
    '''
    parser = argparse.ArgumentParser( description='take looks')
    parser.add_argument('-i', '--input', dest='input', type=str, required=True,
            help = 'input file')
    parser.add_argument('-o', '--output', dest='output', type=str, required=True,
            help = 'output file')
    parser.add_argument('-r','--rlks', dest='rlks', type=int, default=4,
            help = 'number of range looks')
    parser.add_argument('-a','--alks', dest='alks', type=int, default=4,
            help = 'number of azimuth looks')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':

    inps = cmdLineParse()

    if inps.input.endswith('.amp'):
        ampLooks(inps)
    elif inps.input.endswith('.int'):
        intLooks(inps)
    elif inps.input.endswith('.msk'):
        mskLooks(inps)
    elif inps.input.endswith('.hgt') or inps.input.endswith('.lat') or inps.input.endswith('.lon'):
        hgtLooks(inps)
    else:
        raise Exception('file type not supported yet')
        #sys.exit()

#look.py -i diff_20130927-20141211.int -o diff_20130927-20141211_16rlks_16alks.int -r 16 -a 16
#look.py -i 20130927-20141211.amp -o 20130927-20141211_16rlks_16alks.amp -r 16 -a 16    
#look.py -i msk.msk -o 3_6.msk -r 3 -a 6



