#!/usr/bin/env python3

#Cunren Liang, 9-APR-2015
#JPL/Caltech
#updated, 04-APR-2017

import os
import sys
import pickle
import datetime
import argparse
import numpy as np

import isce
import isceobj
from contrib.Snaphu.Snaphu import Snaphu
from isceobj.Planet.Planet import Planet
from isceobj.Constants import SPEED_OF_LIGHT

from crlpac import getWidth
from crlpac import runCmd


def runUnwrap(inps, costMode = None,initMethod = None, defomax = None, initOnly = None):

    if costMode is None:
        costMode   = 'DEFO'
    
    if initMethod is None:
        initMethod = 'MST'
    
    if  defomax is None:
        defomax = 4.0
    
    if initOnly is None:
        initOnly = False
    
    #get earth radius and altitude using master's orbit
    with open(inps.mframe, 'rb') as f:
        mframe = pickle.load(f)

    azti = 1.0/mframe.PRF
    tmid = mframe.getSensingStart() + datetime.timedelta(seconds=azti*(mframe.getNumberOfLines()/2.0))
    orbit = mframe.getOrbit()
    peg = orbit.interpolateOrbit(tmid, method='hermite')

    refElp = Planet(pname='Earth').ellipsoid
    llh = refElp.xyz_to_llh(peg.getPosition())
    hdg = orbit.getHeading(tmid)
    #change to the following after updated to the newest version of isce
    #hdg = orbit.getENUHeading(tmid)
    refElp.setSCH(llh[0], llh[1], hdg)
    earthRadius = refElp.pegRadCur
    altitude = llh[2]

        
    wrapName = inps.inf
    unwrapName = inps.unw
    corrfile = inps.cor

    width = getWidth(wrapName + '.xml')
    wavelength = mframe.getInstrument().getRadarWavelength()
    rangeLooks = inps.rlks
    azimuthLooks = inps.alks

    #calculate azimuth resolution and pixel size
    azres = mframe.platform.antennaLength / 2.0
    vel = np.sqrt(peg.velocity[0]*peg.velocity[0]+peg.velocity[1]*peg.velocity[1]+peg.velocity[2]*peg.velocity[2])
    hgt = np.sqrt(peg.position[0]*peg.position[0]+peg.position[1]*peg.position[1]+peg.position[2]*peg.position[2])
    azimuthPixelSpacing = earthRadius / hgt * vel * azti
    azfact = azres / azimuthPixelSpacing

    #calculate range resolution and pixel size
    rBW = mframe.instrument.pulseLength * mframe.instrument.chirpSlope
    rgres = abs(SPEED_OF_LIGHT / (2.0 * rBW))
    slantRangePixelSpacing = 0.5 * SPEED_OF_LIGHT / mframe.rangeSamplingRate
    rngfact = rgres / slantRangePixelSpacing

    print('azfact: {}, rngfact: {}'.format(azfact, rngfact))
    corrLooks = azimuthLooks * rangeLooks / (azfact * rngfact)

    maxComponents = 20

    snp = Snaphu()
    snp.setInitOnly(initOnly) # follow
    snp.setInput(wrapName)
    snp.setOutput(unwrapName)
    snp.setWidth(width)
    snp.setCostMode(costMode) # follow
    snp.setEarthRadius(earthRadius)
    snp.setWavelength(wavelength)
    snp.setAltitude(altitude)
    snp.setCorrfile(corrfile)
    snp.setInitMethod(initMethod) # follow
    snp.setCorrLooks(corrLooks)
    snp.setMaxComponents(maxComponents)
    snp.setDefoMaxCycles(defomax) # follow
    snp.setRangeLooks(rangeLooks)
    snp.setAzimuthLooks(azimuthLooks)
    #snp.setCorFileFormat('FLOAT_DATA')
    snp.prepare()
    snp.unwrap()

    ######Render XML
    outImage = isceobj.Image.createUnwImage()
    outImage.setFilename(unwrapName)
    outImage.setWidth(width)
    outImage.setAccessMode('read')
    outImage.renderVRT()
    outImage.createImage()
    outImage.finalizeImage()
    outImage.renderHdr()

    #####Check if connected components was created
    if snp.dumpConnectedComponents:
        connImage = isceobj.Image.createImage()
        connImage.setFilename(unwrapName+'.conncomp')
        #At least one can query for the name used
        #self._insar.connectedComponentsFilename = unwrapName+'.conncomp'
        connImage.setWidth(width)
        connImage.setAccessMode('read')
        connImage.setDataType('BYTE')
        connImage.renderVRT()
        connImage.createImage()
        connImage.finalizeImage()
        connImage.renderHdr()

    return


def runUnwrapMcf(inps):
    runUnwrap(inps, costMode = 'SMOOTH',initMethod = 'MCF', defomax = 2, initOnly = True)
    return


def maskUnwrapInterf(inps):
    #mask unwrapped interferogram using connected components
    cmd = "imageMath.py -e='a_0*(b>0);a_1*(b>0)' --a={0} --b={1} -s BIL -o={2}".format(inps.unw, inps.unw+'.conncomp', inps.munw)
    runCmd(cmd)

    #also remove some wired things in original unwrapped interferograms
    os.rename(inps.unw, 'tmp.unw')
    os.rename(inps.unw+'.xml', 'tmp.unw.xml')
    os.rename(inps.unw+'.vrt', 'tmp.unw.vrt')
    cmd = "imageMath.py -e='a_0*(abs(b)!=0);a_1*(abs(b)!=0)' --a={0} --b={1} -s BIL -o={2}".format('tmp.unw', inps.inf, inps.unw)
    runCmd(cmd)
    os.remove('tmp.unw')
    os.remove('tmp.unw.xml')
    os.remove('tmp.unw.vrt')

    return


def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser( description='Unwrap Interferogram Using Snaphu (MCF)')
    parser.add_argument('-mframe', dest='mframe', type=str, required=True,
            help = 'master frame')
    parser.add_argument('-inf', dest='inf', type=str, required=True,
            help = 'wrapped Interferogram')
    parser.add_argument('-cor', dest='cor', type=str, required=True,
            help = 'coherence file')
    parser.add_argument('-unw', dest='unw', type=str, required=True,
            help = 'unwrapped Interferogram')
    parser.add_argument('-munw', dest='munw', type=str, required=True,
            help = 'unwrapped and masked Interferogram')
    parser.add_argument('-rlks', dest='rlks', type=int, default=1,
            help = 'Number of range looks')
    parser.add_argument('-alks', dest='alks', type=int, default=1,
            help = 'Number of azimuth looks')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':

    inps = cmdLineParse()

    runUnwrapMcf(inps)
    maskUnwrapInterf(inps)
