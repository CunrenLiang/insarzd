#!/usr/bin/env python3

import os
import sys
import ntpath
import pickle
import argparse
import datetime

import isce
import isceobj
from isceobj.Constants import SPEED_OF_LIGHT


def cmdLineParse():
    '''
    Command line parser.
    '''
    parser = argparse.ArgumentParser( description='update frame for mosaicked swath')
    parser.add_argument('-s', '--startframe', dest='startframe', type=str, required=True,
            help = 'starting subswath frame')
    parser.add_argument('-e', '--endframe', dest='endframe', type=str, required=True,
            help = 'ending subswath frame')
    parser.add_argument('-o', '--outputframe', dest='outputframe', type=str, required=True,
            help = 'output frame name in each subswath directory. A parameter file is also generated at the same time.')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':

    inps = cmdLineParse()

    #here 1 represents the starting subswath, may not actually be subswath 1
    #5 represents the ending subswath, may not actually be subswath 5
    pickleFile1 = inps.startframe
    pickleFile5 = inps.endframe
    
    with open(pickleFile5, 'rb') as fid:
        frame5 = pickle.load(fid)
        rangePixelSize5 = 0.5 * SPEED_OF_LIGHT / frame5.instrument.rangeSamplingRate
        startingRange5 = frame5.getStartingRange()
        numberOfSamples5 = frame5.getNumberOfSamples()
        sensingStart5 = frame5.getSensingStart()
        sensingStop5 = frame5.getSensingStop()

    with open(pickleFile1, 'rb') as fid:
        frame1 = pickle.load(fid)
        rangePixelSize1 = 0.5 * SPEED_OF_LIGHT / frame1.instrument.rangeSamplingRate
        startingRange1 = frame1.getStartingRange()
        sensingStart1 = frame1.getSensingStart()
        sensingStop1 = frame1.getSensingStop()

        numberOfLines1 = frame1.getNumberOfLines()
        azimuthLineInterval1 = 1.0/frame1.PRF



    #update file width here
    numberOfSamples = int(  (startingRange5 + numberOfSamples5 * rangePixelSize5 - startingRange1) / rangePixelSize1  )
    print('updating file width as: {}'.format(numberOfSamples))
    frame1.setNumberOfSamples(numberOfSamples)

    if sensingStart5 < sensingStart1:
        print('updating sensingStart: replace {} with {}'.format(sensingStart1, sensingStart5))
        frame1.setSensingStart(sensingStart5)
        numberOfLines = numberOfLines1 + int( (sensingStart1 - sensingStart5).total_seconds()/azimuthLineInterval1 ) + 1
        print('updating numberOfLines: replace {} with {}'.format(numberOfLines1, numberOfLines))
        frame1.setNumberOfLines(numberOfLines)

    if sensingStop5 > sensingStop1:
        print('updating sensingStop: replace {} with {}'.format(sensingStop1, sensingStop5))
        frame1.setSensingStop(sensingStop5)
        numberOfLines = numberOfLines1 + int( (sensingStart5 - sensingStart1).total_seconds()/azimuthLineInterval1 ) + 1
        print('updating numberOfLines: replace {} with {}'.format(numberOfLines1, numberOfLines))
        frame1.setNumberOfLines(numberOfLines)


    pickleFileNew = inps.outputframe
    #sensing mid
    frame1.setSensingMid(frame1.getSensingStart() + datetime.timedelta(seconds=(frame1.getSensingStop()-frame1.getSensingStart()).total_seconds()/2.0))
    with open(pickleFileNew, 'wb') as f:
        pickle.dump(frame1, f)












