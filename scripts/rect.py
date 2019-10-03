#!/usr/bin/env python3

import os
import sys
import shutil
import argparse

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

    parser = argparse.ArgumentParser( description='adjust an image by affine transformation and interpolation')
    parser.add_argument('-i', '--inimage', dest='inimage', type=str, required=True,
            help = 'input image')
    parser.add_argument('-o', '--outimage', dest='outimage', type=str, required=True,
            help = 'output image')
    parser.add_argument('-t','--aff', dest='aff', type=str, required=True,
            help = 'affine transformation file')
    parser.add_argument('-f', '--format', dest='format', type=str, required=True,
            help = 'image file format: CPX, RMG, BYTE, REAL, DOUBLE')
    parser.add_argument('-r','--rlks', dest='rlks', type=int, default=1,
            help = 'range looks of the affine transformation')
    parser.add_argument('-a', '--alks', dest='alks', type=int, default=1,
            help = 'azimuth looks of the affine transformation')
    parser.add_argument('-r0','--rlks0', dest='rlks0', type=int, default=1,
            help = 'range looks of the image')
    parser.add_argument('-a0', '--alks0', dest='alks0', type=int, default=1,
            help = 'azimuth looks of the image')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':
    
    inps = cmdLineParse()

    inWidth = getWidth(inps.inimage + '.xml')
    inLength = getLength(inps.inimage + '.xml')
    outWidth = inWidth
    outLength = inLength

    with open(inps.aff) as f:
        lines = f.readlines()
        f.close

    i = 0
    for linex in lines:
        if 'Affine Matrix ' in linex:
            m11 = float(lines[i + 2].split()[0])
            m12 = float(lines[i + 2].split()[1])
            m21 = float(lines[i + 3].split()[0])
            m22 = float(lines[i + 3].split()[1])
            t1  = float(lines[i + 7].split()[0])
            t2  = float(lines[i + 7].split()[1])
            break
        i += 1    

    formatx = inps.format.upper()
    if formatx in ['CPX']:
        methodx = 'Sinc'
    elif formatx in ['BYTE']:
        methodx = 'NN' #Nearest Neighbor
    elif formatx in ['RMG', 'REAL', 'DOUBLE']:
        methodx = 'Bilinear'
    else:
        raise Exception('File type not supported!')

    rectInput = '''Input Image File Name  (-) = {inFile}        ! dimension of file to be rectified
Output Image File Name (-) = {outFile}       ! dimension of output
Input Dimensions       (-) = {inWidth} {inLength} ! across, down
Output Dimensions      (-) = {outWidth} {outLength} ! across, down
Affine Matrix Row 1    (-) = {m11} {m12}      ! a b
Affine Matrix Row 2    (-) = {m21} {m22}      ! c d
Affine Offset Vector   (-) = {t1} {t2}        ! e f
Looks of the Offsets   (-) = {rlks} {alks}        ! lac, ldn
Looks of the Image File   (-) = {rlks0} {alks0}        ! lac0, ldn0
File Type              (-) = {format}        ! [RMG, COMPLEX]
Interpolation Method   (-) = {method}        ! [NN, Bilinear, Sinc]
'''.format(
    inFile = inps.inimage,
    outFile = inps.outimage,
    inWidth = inWidth, inLength = inLength,
    outWidth = outWidth, outLength = outLength,
    m11 = m11, m12 = m12,
    m21 = m21, m22 = m22,
    t1 = t1, t2 = t2,
    rlks = inps.rlks, alks = inps.alks,
    rlks0 = inps.rlks0, alks0 = inps.alks0,
    format = formatx,
    method = methodx
    )

    rectInputFile = 'rect.in'
    with open(rectInputFile, 'w') as ff:
        ff.write(rectInput)

    #run fitoff here
    cmd = '$INSAR_ZERODOP_BIN/rect_with_looks {}'.format(rectInputFile)
    #print("{}".format(cmd))
    runCmd(cmd)


    #get xml file
    shutil.copy(inps.inimage + '.xml', inps.outimage + '.xml')
    #fix file name
    img, dataName, metaName = IML.loadImage(inps.outimage + '.xml')
    img.filename = inps.outimage
    #img.setAccessMode('READ')
    img.renderHdr()


#./rect.py -i 141018-141123_rg.off -o rect_141018-141123_rg.off -c ampsim_16rlks_16alks.aff -f real -r 16 -a 16





    
















