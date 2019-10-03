#!/usr/bin/env python3

#Cunren Liang, Apr. 20, 2017
#JPL/Caltech


import os
import sys
import pickle
import argparse


import isce
import isceobj

from crlpac import renderParXml


def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(description="create parameter xml file")
    parser.add_argument('-frame', dest='frame', type=str, required=True,
            help='frame pickle file used to create parameter xml file')
    parser.add_argument('-par', dest='par', type=str, required=True,
            help='output parameter xml file')

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

    with open(inps.frame, 'rb') as f:
        frame = pickle.load(f)

    renderParXml(frame, os.path.splitext(inps.par)[0])

