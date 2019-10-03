#!/usr/bin/env python3

#Cunren Liang, 24-Apr-2017
#JPL/Caltech

########################################################################################
#after subswath mosaicking, the starting range of the mosaicked pickle is the same as
#starting swath, but the starting time may be different from the starting swath.
#the frame offset is estimated using starting swath images, so the azimuth frame offset
#need to be updated.
########################################################################################

import os
import re
import sys
import argparse
import datetime
import pickle
import numpy as np

import isce
import isceobj


def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="update frame offsets\
                    \n\
                    ")
    parser.add_argument('-offset', dest='offset', type=str, required=True,
            help = 'subswath offset file')
    parser.add_argument('-f', action='append', dest='frame', default=[],
            help = 'starting subswath frame pickle files of master IN ORDER')
    parser.add_argument('-f2', action='append', dest='frame2', default=[],
            help = 'mosaicked frame pickle files of master IN ORDER')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':

    inps = cmdLineParse()

    frames = inps.frame
    frames2 = inps.frame2

    if not (len(frames) == len(frames2)):
        raise Exception('the numbers of files are not equal!')
    if len(frames) < 2:
        raise Exception('at least 2 frames should be specified!')

    nframe = len(frames)

    ############################
    pcks = []
    pcks2= []

    sensingStart = []
    sensingStart2= []
    azimuthLineInterval = []
    ############################

    for frame in frames:
        with open(frame, 'rb') as f:
            pcks.append(pickle.load(f))
    for frame in frames2:
        with open(frame, 'rb') as f:
            pcks2.append(pickle.load(f))

    for pck in pcks:
        sensingStart.append(pck.getSensingStart())
        azimuthLineInterval.append(1.0/pck.PRF)

    for pck in pcks2:
        sensingStart2.append(pck.getSensingStart())


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


    if not (nframe == len(rangeOffset1)):
        raise Exception('number of frames is not equal to number of frame offsets')


    #update azimuth offset
    for i in range(1, nframe):
        azTimeDiff_0 = (sensingStart[i-1] - sensingStart2[i-1]).total_seconds()
        azTimeDiff_1 = (sensingStart[i] - sensingStart2[i]).total_seconds()
        if azTimeDiff_0 != 0.0:
            azimuthOffset1[i] -= azTimeDiff_0 / azimuthLineInterval[i-1]
            azimuthOffset2[i] -= azTimeDiff_0 / azimuthLineInterval[i-1]

        if azTimeDiff_1 != 0.0:
            azimuthOffset1[i] += azTimeDiff_1 / azimuthLineInterval[i]
            azimuthOffset2[i] += azTimeDiff_1 / azimuthLineInterval[i]

    #write new offsets
    os.remove(inps.offset)
    #######################################################################################
    offsetComp = "\n\ncomparision of offsets:\n\n"
    offsetComp += "offset type       i     geometrical           match      difference\n"
    offsetComp += "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
    for i, (offset1, offset2) in enumerate(zip(rangeOffset1, rangeOffset2)):
        offsetComp += "range offset     {:2d}   {:13.3f}   {:13.3f}   {:13.3f}\n".format(i, offset1, offset2, offset1 - offset2)
    for i, (offset1, offset2) in enumerate(zip(azimuthOffset1, azimuthOffset2)):
        offsetComp += "azimuth offset   {:2d}   {:13.3f}   {:13.3f}   {:13.3f}\n".format(i, offset1, offset2, offset1 - offset2)

    #write and report offsets
    with open(inps.offset, 'w') as f:
        f.write(offsetComp)
    print("{}".format(offsetComp))


