#!/usr/bin/env python3

import os
import sys
import shutil
import argparse

from crlpac import runCmd


def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(description='replace the magnitude of unwrapped interferogram with magnitude of original interferogram')
    parser.add_argument('-unw', type=str, dest='unw', required=True,
            help='unwrapped interferogram')
    parser.add_argument('-inf', type=str, dest='inf', required=True,
            help='wrapped interferogram')
    parser.add_argument('-msk', dest='msk', type=int, default=0,
            help = 'whether mask the area not unwrapped? 0: yes (default). 1: no')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':

    inps = cmdLineParse()

    tmpf = 'tmp.unw'
    os.rename(inps.unw, tmpf)
    os.rename(inps.unw + '.xml', tmpf + '.xml')
    os.rename(inps.unw + '.vrt', tmpf + '.vrt')

    if inps.msk == 0:
        cmd = "imageMath.py -e='(a_0>0)*abs(b);a_1' --a={} --b={} -t float -s BIL -o={}".format(
          tmpf,
          inps.inf,
          inps.unw
        )
    else:
        cmd = "imageMath.py -e='abs(b);a_1' --a={} --b={} -t float -s BIL -o={}".format(
          tmpf,
          inps.inf,
          inps.unw
        )
    runCmd(cmd)

    os.remove(tmpf)
    os.remove(tmpf+'.xml')
    os.remove(tmpf+'.vrt')

