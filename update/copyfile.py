#!/usr/bin/env python3

#Cunren Liang, Aug. 15, 2016
#JPL/Caltech


import os
import sys
import glob
import shutil
import argparse


def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(description="copy all the files in this update directory to isce directory.\
             run this command under this directory")
    parser.add_argument('-idir', dest='idir', type=str, required=True,
            help='isce directory, under which there are folders applications, components, configuration, contrib...')

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

    i = 0
    for root, dirnames, filenames in os.walk('.'):
        for x in filenames:
            if x not in ['.DS_Store', 'copyfile.py', 'readme.txt']:
                src = os.path.join(root, x)
                dst = os.path.join(os.path.abspath(inps.idir), os.path.relpath(src, './'))

                dst_dir = os.path.split(dst)[0]
                if not os.path.exists(dst_dir):
                    print('{} does not exist, creat...'.format(dst_dir))
                    os.makedirs(dst_dir)

                print('copying from: {}'.format(src))
                print('          to: {}\n'.format(dst))

                shutil.copy(src, dst)
                i += 1
    print('total number of files copied: {}'.format(i))
        


#mkdir directory recursively
#mkdir -p 




