#!/usr/bin/env python3

#Cunren Liang, JPL/Caltech


import os
import sys
import glob
import argparse

def cmdLineParse():
    '''
    Command line parser.
    '''
    parser = argparse.ArgumentParser( description='delete file')
    parser.add_argument('-f', dest='f', type=str, required=True,
            help = 'file to be deleted')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':

    inps = cmdLineParse()
   
    files = sorted(glob.glob(inps.f)) 
    #print(files)
    #print(len(files))

    if len(files) == 0:
        print('no file is to be deleted regarding to: {}'.format(inps.f))
    else:
        for file in files:
            os.remove(file)


