#!/usr/bin/env python3

#program for cropping slc keeping only the overlap area
#works with different PRFs and range sampling rates
#Cunren Liang, JPL/Caltech

import os
import sys
import pickle
import datetime
import argparse
import numpy as np

import isce
import isceobj
from isceobj.Constants import SPEED_OF_LIGHT
from zerodop.baseline.Baseline import Baseline

from crlpac import create_xml


def cmdLineParse():
    '''
    Command line parser.
    '''
    parser = argparse.ArgumentParser( description='crop slc keeping only the overlap area.')
    parser.add_argument('-slc', dest='slc', type=str, required=True,
            help = 'slc file')
    parser.add_argument('-frame', dest='frame', type=str, required=True,
            help = 'frame pickle file')
    parser.add_argument('-frame2', dest='frame2', type=str, required=True,
            help = 'the other frame pickle file')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':

    inps = cmdLineParse()

    #frame
    with open(inps.frame, 'rb') as f:
        frame = pickle.load(f)
    #the other frame
    with open(inps.frame2, 'rb') as f:
        frame2 = pickle.load(f)


    #take frame2 as master
    bObj = Baseline()
    bObj.configure()
    bObj.baselineLocation = 'top'
    bObj.wireInputPort(name='masterFrame', object=frame2)
    bObj.wireInputPort(name='slaveFrame', object=frame)
    azoff0, rgoff0 = bObj.baseline()

    azoff = int(azoff0[0])
    rgoff = int(rgoff0[0])

    #from now all column and line numbers are based on the original 'frame'

    #here we use the convention that range/azimuth sample index starting from 0
    #first line of frame2 in frame
    fl2 = 0 + azoff
    #last line of frame2 in frame
    ll2 = round((frame2.getNumberOfLines()-1.0) / frame2.instrument.PRF * frame.instrument.PRF) + azoff
    #first column of frame2 in frame
    fc2 = 0 + rgoff
    #last column of frame2 in frame
    lc2 = round((frame2.getNumberOfSamples()-1.0) / frame2.instrument.rangeSamplingRate * frame.instrument.rangeSamplingRate) + rgoff

    #frame
    fl = 0
    ll = frame.getNumberOfLines()-1
    fc = 0
    lc = frame.getNumberOfSamples()-1

    #frame3, output
    fl3 = max(fl, fl2)
    ll3 = min(ll, ll2)
    fc3 = max(fc, fc2)
    lc3 = min(lc, lc2)

    #check if there is overlap
    if ll3 - fl3 < 1000 or lc3 - fc3 < 1000:
        raise Exception('there is not enough overlap area. nothing is done by this program.\n')

    #check if need to crop
    if abs(fl3-fl) < 100 and abs(ll3-ll) < 100 and abs(fc3-fc) < 100 and abs(lc3-lc) < 100:
        print('no need to crop {}. nothing is done by this program.'.format(inps.slc))
    else:
        #read and crop data
        with open(inps.slc, 'rb') as f:
            f.seek((fl3 - fl) * frame.getNumberOfSamples() * np.dtype(np.complex64).itemsize, 0)
            data = np.fromfile(f, dtype=np.complex64, count=(ll3 - fl3 + 1) * frame.getNumberOfSamples()).reshape((ll3 - fl3 + 1),frame.getNumberOfSamples())
            data2 = data[:, fc3:lc3+1]

        #remove original
        os.remove(inps.slc)
        #os.remove(inps.slc + '.par.xml')
        os.remove(inps.slc + '.xml')
        os.remove(inps.slc + '.vrt')
        os.remove(inps.frame)

        #write new
        with open(inps.slc, 'wb') as f:
            data2.astype(np.complex64).tofile(f)

        #update parameters
        outWidth = lc3 - fc3 + 1
        outLength = ll3 - fl3 + 1

        index2  = np.array(range(outWidth))
        index = np.array(range(frame.getNumberOfSamples()))
        dop = frame.dopCoeff[0] + frame.dopCoeff[1] * index + frame.dopCoeff[2] * index**2 + frame.dopCoeff[3] * index**3
        p = np.polyfit(index2, dop[fc3:lc3+1], 3) 
        dopCoeff=list(np.flipud(p))

        ka = frame.fmrateCoeff[0] + frame.fmrateCoeff[1] * index**1 + frame.fmrateCoeff[2] * index**2
        p = np.polyfit(index2, ka[fc3:lc3+1], 2)
        fmrateCoeff=list(np.flipud(p))

        #set output frame
        outframe = frame
        outframe.setNumberOfSamples(outWidth)
        outframe.setNumberOfLines(outLength)
        outframe.setStartingRange(frame.startingRange + (fc3-fc)*0.5 * SPEED_OF_LIGHT / frame.rangeSamplingRate)
        outframe.setFarRange(outframe.startingRange + (outWidth-1)*0.5 * SPEED_OF_LIGHT / frame.rangeSamplingRate)
        outframe.setSensingStart(frame.getSensingStart() + datetime.timedelta(seconds=(fl3 - fl) / frame.instrument.PRF))
        outframe.setSensingStop(outframe.getSensingStart() + datetime.timedelta(seconds=(outLength - 1.0) / frame.instrument.PRF))
        outframe.setSensingMid(outframe.getSensingStart() + datetime.timedelta(seconds=(outLength - 1.0) / frame.instrument.PRF/2.0))
        outframe.dopCoeff = dopCoeff
        outframe.fmrateCoeff = fmrateCoeff
        if hasattr(frame, 'burstStartLines'):
            d=np.asarray(frame.burstStartLines) - (fl3 - fl)
            outframe.burstStartLines = d.tolist()
        if hasattr(frame, 'burstStartLineEstimated'):
            outframe.burstStartLineEstimated = frame.burstStartLineEstimated - (fl3 - fl)

        create_xml(inps.slc, outWidth, outLength, 'slc')
        #sensing mid
        outframe.setSensingMid(outframe.getSensingStart() + datetime.timedelta(seconds=(outframe.getSensingStop()-outframe.getSensingStart()).total_seconds()/2.0))
        with open(inps.frame, 'wb') as f:
            pickle.dump(outframe, f)













