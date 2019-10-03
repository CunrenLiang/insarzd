#!/usr/bin/env python3

#This program computes azimuth FM rate using frame pickle file.
#The computed azimuth FM rate and range sample number (start with index 0)
#are fitted to a second polynomial.
#The coefficents of the polynomial is saved to the frame.
#The old frame file is then replaced with the new frame file.

#Cunren Liang, 21-OCT-2015


import os
import re
import sys
import pickle
import argparse
import datetime
import subprocess
import numpy as np

import isce
import isceobj
import stdproc
from isceobj.Location.Peg import Peg #not sure if we need this
from isceobj.Orbit.Orbit import Orbit
from isceobj.Constants import SPEED_OF_LIGHT
from isceobj.Planet.Planet import Planet
from stdproc.orbit.pegManipulator import averagePeg
from iscesys.StdOEL.StdOELPy import create_writer

from crlpac import runCmd


def runCmdCaptureOutput(cmd):
    try:
        output0 = subprocess.check_output(cmd)
    except subprocess.CalledProcessError as errorCode:                                                                                                   
        raise Exception('error when running:\n{}\nerror code: {}, {}\n'.format(cmd, errorCode.returncode, errorCode.output))

    output = output0.decode('utf-8')

    return output


#The following two functions are not used in this program
#But I just put it here, because they are useful.
def pulseTiming(frame, firstLine, lastLine):

    if lastLine < firstLine:
        raise Exception('lastLine < firstLine!')
    
    print("Do pulse timing for "+frame.getSensingStart().strftime('%Y-%m-%d'))
    
    
    startTime = frame.getSensingStart()
    orbit = frame.getOrbit()
    
    pulseOrbit = Orbit()
    for i in range(firstLine, lastLine + 1):
        dt = (i - 1) * pri
        time = startTime + datetime.timedelta(seconds=dt)
        sv = orbit.interpolateOrbit(time, method='hermite')
        pulseOrbit.addStateVector(sv)

    return pulseOrbit


def getPegpoint(planet, orbit):

    stdWriter = create_writer("log", "", True, filename="pegpoint.log")
    stdWriter.setFileTag("getpeg", "log")
    stdWriter.setFileTag("getpeg", "err")
    stdWriter.setFileTag("getpeg", "out")

    objGetpeg = stdproc.createGetpeg()
    objGetpeg.wireInputPort(name='planet', object=planet)
    objGetpeg.wireInputPort(name='Orbit', object=orbit)
    objGetpeg.setStdWriter(stdWriter)

    objGetpeg.estimatePeg()

    return objGetpeg



def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser( description='program for computing azimuth FM rate')
    parser.add_argument('-f', dest='frame', type=str, required=True,
            help = 'frame')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':

    inps = cmdLineParse()

    BIDIR=os.environ.get('INSAR_ZERODOP_BIN')

    with open(inps.frame, 'rb') as f:
        frame = pickle.load(f)


    #interpolate at center of the scene to get the peg point
    pri = 1.0 / frame.getInstrument().getPulseRepetitionFrequency()
    dt = frame.getNumberOfLines() / 2.0 * pri
    time = frame.getSensingStart() + datetime.timedelta(seconds=dt)
    sv = frame.getOrbit().interpolateOrbit(time, method='hermite')

    planet=Planet(pname='Earth')
    earthGM=planet.get_GM()
    earthSpin=planet.get_spin()
    #Scott Hensley's value from get_peg_info.f
    earthSpin=7.29211573052e-5


    pos = sv.getPosition()
    vel = sv.getVelocity()

    #cmd = "['{}/get_peg', '{}', '{}', '{}', '{}', '{}', '{}', '{}', '{}']".format(BIDIR, pos[0], pos[1], pos[2], vel[0], vel[1], vel[2], earthGM, earthSpin)
    cmd = "{}/get_peg {} {} {} {} {} {} {} {}".format(BIDIR, pos[0], pos[1], pos[2], vel[0], vel[1], vel[2], earthGM, earthSpin)
    cmd = re.split('\s+', cmd)

    output = runCmdCaptureOutput(cmd)
    for linex in output.split("\n"):
        if 'Platform SCH Velocity' in linex:
            velocity = linex.split(':')[1].strip(' ').strip('\n')
            velocity = re.split('\s+', velocity)
            velocity = [float(velocity[0]), float(velocity[1]), float(velocity[2])]
            #print(velocity)
        if 'Platform SCH Acceleration' in linex:
            acceleration = linex.split(':')[1].strip(' ').strip('\n')
            acceleration = re.split('\s+', acceleration)
            acceleration = [float(acceleration[0]), float(acceleration[1]), float(acceleration[2])]
            #print(acceleration)
        if 'Peg Lat/Lon , H' in linex:
            peg = linex.split('=')[1].strip(' ').strip('\n')
            peg = re.split('\s+', peg)
            peg = [float(peg[0]), float(peg[1]), float(peg[2])]
            #print(peg)
        if 'Radius Curvature' in linex:
            earthRadius = linex.split('=')[1].strip(' ').strip('\n')
            earthRadius = float(earthRadius)
            #print(earthRadius)

    centerVel = velocity
    centerAcc = acceleration
    avghgt = peg[2]
    radiusOfCurvature = earthRadius



    fmrate = []
    width = frame.getNumberOfSamples()
    lookSide = frame.getInstrument().getPlatform().pointingDirection # right: -1 left: 1
    centerVelNorm=np.sqrt(centerVel[0]**2 + centerVel[1]**2 + centerVel[2]**2)
    for i in range(width):
        rg = frame.startingRange + i * 0.5 * SPEED_OF_LIGHT / frame.rangeSamplingRate
        dop = frame.dopCoeff[0] + frame.dopCoeff[1] * i + frame.dopCoeff[2] * i * i + frame.dopCoeff[3] * i * i * i
        dop *= frame.PRF
        
        th = np.arccos(((avghgt+radiusOfCurvature)**2 + rg**2 - radiusOfCurvature**2) / (2.0 * (avghgt+radiusOfCurvature) * rg))
        thaz = np.arcsin(((frame.radarWavelegth*dop/(2.0*np.sin(th)))+(centerVel[2]/np.tan(th))) / np.sqrt(centerVel[0]**2+centerVel[1]**2))-lookSide*np.arctan(centerVel[1]/centerVel[0])

        lookVec = []
        lookVec.append(np.sin(th)*np.sin(thaz))
        lookVec.append(np.sin(th)*np.cos(thaz)*lookSide)
        lookVec.append(-np.cos(th))

        vdotl = lookVec[0] * centerVel[0] + lookVec[1] * centerVel[1] + lookVec[2] * centerVel[2]
        adotl = lookVec[0] * centerAcc[0] + lookVec[1] * centerAcc[1] + lookVec[2] * centerAcc[2]

        fmratex = 2.0*(adotl + (vdotl**2 - centerVelNorm**2)/rg)/(frame.radarWavelegth)
        fmrate.append(fmratex)

    #check computing results
    #for fmratex in fmrate:
    #    print("{} {}".format(i, fmratex))

###############################################################################
#test result for sensor: ALOS2, date: 150725, frame: f540, track: a157
#                            near                 far
#value from leader file:     -581.086931          -547.5056850000001
#value computed:             -581.070931381       -549.246329448

#test result for sensor: Sentinel-1A, date: 141123, over Napa
#                            near                 far
#value from xml file:        -2526.21927231       -2414.21524984
#value computed:             -2526.1041501223403  -2414.143480982222
###############################################################################

    print('fit the azimuth FM rate to a 2nd polynomial')

    #fit the azimuth FM rate to a 2nd polynomial
    rangeBins = list(range(width))
    p = np.polyfit(rangeBins, fmrate, 2)
    fmrateCoeff = list(p)
    fmrateCoeff.reverse()

    #check fitting results
    print('\n  rgbin    ori      fit        dif')
    print('**************************************')
    for i in range(width):
        testi = fmrateCoeff[2] * i**2 + fmrateCoeff[1] * i**1 + fmrateCoeff[0]
        if i % 1000 == 0:
            print("%7d %8.2f %8.2f %8.2f" % (i, fmrate[i], testi, fmrate[i] - testi))

###############################################################################

    frame.fmrateCoeff = fmrateCoeff
    #remove the old frame pickle file
    os.remove(inps.frame)
    #create the new frame pickle file with fmrateCoeff added
    with open(inps.frame, 'wb') as f:
        pickle.dump(frame, f)















