#!/usr/bin/env python3

#Cunren Liang, Aug. 15, 2016
#JPL/Caltech


import os
import sys
import pickle
import shutil
import argparse
import datetime
import numpy as np

import isce
import isceobj

from crlpac import runCmd


def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(description="calulate raw burst start times of SLC focused by full-aperture")
    parser.add_argument('-s', dest='slc', type=str, required=True,
            help='SLC focused by full-aperture algorithm')
    parser.add_argument('-f', dest='frame', type=str, required=True,
            help='frame object of the SLC saved in pickle file')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':
    '''
    Main driver.
    '''

    inps = cmdLineParse()

    #number of lines from start and end of file
    delta_line = 15000

    start_line1 = delta_line
    cmd = '$INSAR_ZERODOP_SCR/burst_time.py -s {} -f {} -l {} -c {}'.format(
        inps.slc, 
        inps.frame, 
        start_line1, 
        1000)
    runCmd(cmd)

    with open(inps.frame, 'rb') as fid:
        frame1 = pickle.load(fid)
    frame = frame1

    #number of burst cycles
    num_nc = np.around((frame1.getNumberOfLines() - delta_line*2) / frame1.ncraw)

    start_line2 = np.around(start_line1 + num_nc * frame1.ncraw)
    cmd = '$INSAR_ZERODOP_SCR/burst_time.py -s {} -f {} -l {:.0f} -c {}'.format(
        inps.slc, 
        inps.frame, 
        start_line2, 
        1000)
    runCmd(cmd)

    with open(inps.frame, 'rb') as fid:
        frame2 = pickle.load(fid)

    LineDiffIndex = 0
    LineDiffMin = np.fabs(frame1.burstStartLineEstimated + frame1.ncraw * LineDiffIndex - frame2.burstStartLineEstimated)
    for i in range(0, 100000):
        LineDiffMinx = np.fabs(frame1.burstStartLineEstimated + frame1.ncraw * i - frame2.burstStartLineEstimated)
        if LineDiffMinx <= LineDiffMin:
            LineDiffMin = LineDiffMinx
            LineDiffIndex = i

    #correct burst cycle value
    frame.ncraw -= (frame1.burstStartLineEstimated + frame1.ncraw * LineDiffIndex - frame2.burstStartLineEstimated) / LineDiffIndex
    os.remove(inps.frame)
    with open(inps.frame, 'wb') as f:
        pickle.dump(frame, f)

    #use correct frame1.ncraw to do final estimation
    start_line1 = np.around(frame1.getNumberOfLines()/2.0)
    cmd = '$INSAR_ZERODOP_SCR/burst_time.py -s {} -f {} -l {:.0f} -c {}'.format(
        inps.slc, 
        inps.frame, 
        start_line1, 
        1000)
    runCmd(cmd)


