#!/usr/bin/env python3

import os
import sys
import shutil
import argparse
import logging

import isce
import isceobj
from mroipac.filter.Filter import Filter
from mroipac.icu.Icu import Icu

from crlpac import getWidth
from crlpac import getLength
from crlpac import runCmd
from crlpac import create_xml

logger = logging.getLogger('isce.insar.runFilter')


def runFilter(inps):
    logger.info("Applying power-spectral filter")

    #get width from the header file of input interferogram
    width = getWidth(inps.intf + '.xml')
    length = getLength(inps.intf + '.xml')

    if shutil.which('psfilt1') != None:
        cmd = "psfilt1 {int} {filtint} {width} {filterstrength} 64 16".format(
               int = inps.intf,
               filtint = inps.fintf,
               width = width,
               filterstrength = inps.alpha
               )
        runCmd(cmd)

        #get xml file for interferogram
        create_xml(inps.fintf, width, length, 'int')
    else:
        #create flattened interferogram
        intImage = isceobj.createIntImage()
        intImage.setFilename(inps.intf)
        intImage.setWidth(width)
        intImage.setAccessMode('read')
        intImage.createImage()

        #create the filtered interferogram
        filtImage = isceobj.createIntImage()
        filtImage.setFilename(inps.fintf)
        filtImage.setWidth(width)
        filtImage.setAccessMode('write')
        filtImage.createImage()

        #create filter and run it
        objFilter = Filter()
        objFilter.wireInputPort(name='interferogram',object=intImage)
        objFilter.wireOutputPort(name='filtered interferogram',object=filtImage)
        objFilter.goldsteinWerner(alpha=inps.alpha)

        intImage.finalizeImage()
        filtImage.finalizeImage()
        del filtImage


    #recreate filt image to read
    filtImage = isceobj.createIntImage()
    filtImage.setFilename(inps.fintf)
    filtImage.setWidth(width)
    filtImage.setAccessMode('read')
    filtImage.createImage()

    #create amplitude image
    ampImage = isceobj.createAmpImage()
    ampImage.setFilename(inps.amp)
    ampImage.setWidth(width)
    ampImage.setAccessMode('read')
    ampImage.createImage()

    #create phase sigma correlation file here
    phsig_tmp = 'tmp.phsig'
    phsigImage = isceobj.createImage()
    phsigImage.setFilename(phsig_tmp)
    phsigImage.setWidth(width)
    phsigImage.dataType='FLOAT'
    phsigImage.bands = 1
    phsigImage.setImageType('cor')#the type in this case is not for mdx.py displaying but for geocoding method
    phsigImage.setAccessMode('write')
    phsigImage.createImage()

    #create icu and run it
    icuObj = Icu(name='insarapp_filter_icu')
    icuObj.configure()
    icuObj.unwrappingFlag = False
    icuObj.icu(intImage = filtImage, ampImage=ampImage, phsigImage=phsigImage)
    
    phsigImage.renderHdr()

    filtImage.finalizeImage()
    ampImage.finalizeImage()
    phsigImage.finalizeImage()

    # #add an amplitude channel to phsig file
    # cmd = "imageMath.py -e='sqrt(a_0*a_1)*(b!=0);b' --a={amp} --b={phsig_tmp} -o {phsig} -s BIL".format(
    #        amp = inps.amp,
    #        phsig_tmp = phsig_tmp,
    #        phsig = inps.phsig
    #        )

    #add an amplitude channel to phsig file
    cmd = "imageMath.py -e='sqrt(a_0*a_1)*(b!=0);b' --a={amp} --b={phsig_tmp} -o {phsig} -s BIL".format(
           amp = inps.amp,
           phsig_tmp = phsig_tmp,
           phsig = inps.phsig
           )


    runCmd(cmd)

    #remove the original phsig file
    os.remove(phsig_tmp)
    os.remove(phsig_tmp + '.xml')
    os.remove(phsig_tmp + '.vrt')

    
    #rename original filtered interferogram
    filt_tmp = 'filt_tmp.int'
    os.rename(inps.fintf, filt_tmp)
    os.rename(inps.fintf + '.xml', filt_tmp + '.xml')
    os.rename(inps.fintf + '.vrt', filt_tmp + '.vrt')
    
    #do the numpy calculations
    #replace the magnitude of the filtered interferogram with magnitude of original interferogram
    #mask output file using layover mask (values 2 and 3).
    # cmd = "imageMath.py -e='a/(abs(a)+(abs(a)==0))*abs(b)*(c<2)' --a={0} --b={1} --c={2} -t CFLOAT -o={3}".format(
    #       filt_tmp,
    #       inps.intf,
    #       inps.msk,
    #       inps.fintf
    #     )

    #replacing magnitude is not good for phase unwrapping using snaphu
    cmd = "imageMath.py -e='a*(b<2)' --a={0} --b={1} -t CFLOAT -o={2}".format(
          filt_tmp,
          inps.msk,
          inps.fintf
        )

    runCmd(cmd)
    
    #remove the original filtered interferogram
    os.remove(filt_tmp)
    os.remove(filt_tmp + '.xml')
    os.remove(filt_tmp + '.vrt')


def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(description="Filter interferogram using Richard M. Goldstein and Charles L. Werner's power spectrum filter, \
compute phase sigma for the filtered interferogram, and mask out the layerover, shadow, etc areas")

    parser.add_argument('-i','--intf', dest='intf', type=str,
            required=True, help='input interferogram to be filtered')
    parser.add_argument('-m','--amp', dest='amp', type=str,
            required=True, help='input amplitude image')
    parser.add_argument('-k','--msk', dest='msk', type=str,
            required=True, help='layerover, shadow, etc mask')
    parser.add_argument('-f','--fintf', dest='fintf', type=str,
            required=True, help='ouput filtered interferogram')
    parser.add_argument('-p', '--phsig', dest='phsig', type=str,
            required=True, help='output phase sigma correlation file')  
    parser.add_argument('-a', '--alpha', dest='alpha', type=float, default=0.75,
            required=False, help='filtering strength. (Default: 0.75)')

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
    runFilter(inps)


#./filter.py -i diff_20130927-20141211_16rlks_16alks.int -m 20130927-20141211_16rlks_16alks.amp -f filt_diff_20130927-20141211_16rlks_16alks.int -p filt_diff_20130927-20141211_16rlks_16alks.phsig -a 0.5

#new example after revision. 28-AUG-2015
#filter.py -i diff_150221-150502_3rlks_6alks.int -m 150221-150502_3rlks_6alks.amp -f filt_diff_150221-150502_3rlks_6alks.int -p 150221-150502_3rlks_6alks.phsig -a 0.4 -k msk_3_6.rdr -s 1.3









