#!/usr/bin/env python3

#Cunren Liang, Apr. 20, 2015
#JPL/Caltech

import os
import sys
import pickle
import argparse
from xml.etree.ElementTree import ElementTree

import isce
import isceobj
from isceobj.Sensor import createSensor
from isceobj.Doppler import createDoppler

from crlpac import runCmd
from crlpac import getInfo


def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(description='extract data from standard products.')
    parser.add_argument('-s', '--sensor', dest='sensor', type=str, required=True,
            help = 'sensor name: ALOS, ALOS2, COSMO_SKYMED, COSMO_SKYMED_SLC, ENVISAT, ERS, KOMPSAT5, RADARSAT1, RADARSAT2, ROI_PAC, TERRASARX, RISAT1, UAVSAR_RPI, UAVSAR_STACK, SENTINEL1A')
    parser.add_argument('-i', '--input', dest='input', type=str, required=True,
            help = 'input xml configure file')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


def getS1aDoppler(xmlfile):
    from isceobj.Planet.AstronomicalHandbook import Const
    import math
    import numpy
    
    xmlfp = None
    try:
        xmlfp = open(xmlfile,'r')
        print('getting Doppler centroid from: {0}'.format(xmlfile))
        xmlx = ElementTree(file=xmlfp).getroot()

        # read some parameters
        azimuthTimeInterval = float(xmlx.find('imageAnnotation/imageInformation/azimuthTimeInterval').text)
        prf = 1.0 / azimuthTimeInterval
        prf2 = float(xmlx.find('generalAnnotation/downlinkInformationList/downlinkInformation/prf').text)
        #This is starting range
        slantRangeTime = float(xmlx.find('imageAnnotation/imageInformation/slantRangeTime').text)        
        rangePixelSize = float(xmlx.find('imageAnnotation/imageInformation/rangePixelSpacing').text)
        rangeSamplingRate = Const.c/(2.0*rangePixelSize)
        #number of samples
        samples = int(xmlx.find('imageAnnotation/imageInformation/numberOfSamples').text)

        #summarize parameters here
        #print("\n**********************************")
        print('')
        print('PRF: {}, {}'.format(prf, prf2))
        print('starting range time: {}'.format(slantRangeTime))
        print('range sampling rate: {}'.format(rangeSamplingRate))
        print('number of samples: {}\n'.format(samples))

       
        dcEstimateList = xmlx.find('dopplerCentroid/dcEstimateList')
        estimationCount = int(dcEstimateList.get('count'))
        print('number of DC estimation: {}'.format(estimationCount))
        estimationCount = len(dcEstimateList.findall('dcEstimate'))
        print('number of DC estimation: {}'.format(estimationCount))
        
        dc=[]
        print("list of DC polynomials:")
        for dcEstimate in dcEstimateList.findall('dcEstimate'):
            t0 = float(dcEstimate.find('t0').text)
            dataDcPolynomial = dcEstimate.find('dataDcPolynomial').text.split()
            #convert to float
            dataDcPolynomialCoeff = []
            for kk in dataDcPolynomial:
                dataDcPolynomialCoeff.append(float(kk))
            dataDcPolynomial = dataDcPolynomialCoeff
            print(t0, dataDcPolynomial)

            dci = []
            for ii in range(samples):
                dci.append(0.0)
                xx = (slantRangeTime + ii / rangeSamplingRate - t0)
                for jj in range(len(dataDcPolynomial)):
                    dci[ii] += dataDcPolynomial[jj] * math.pow(xx, jj)
                #print("{}".format(dci[ii]))
            dc.append(dci)

        #average the Doppler
        dcAverage = []
        for ii in range(samples):
            dcAverage.append(0.0)
            for jj in range(estimationCount):
                dcAverage[ii] += dc[jj][ii] / estimationCount
            #if(ii / 1000 == 0)
            #    print("{}  {}  {}  {}  {}".format(ii, dc[0][ii], dc[1][ii], dc[2][ii], dcAverage[ii]))
        
        #fitting averaged doppler
        rangeBins = list(range(samples)) #[0:samples-1]
        #fitting degree: 3
        p = numpy.polyfit(rangeBins, dcAverage, 3)
        dopCoeff = []
        dopCoeff.append(p[3] / prf)
        dopCoeff.append(p[2] / prf)
        dopCoeff.append(p[1] / prf)
        dopCoeff.append(p[0] / prf)

        #just print out for debug
        ######################################################################################################################
        DEBUG = None
        if DEBUG != None:
            test =[]
            print('\nrgbin    dcl[1]   dcl[2]   dcl[3]    dcave    dcfit')
            print('******************************************************')
            for ii in range(samples):
                testi = dopCoeff[3] * math.pow(ii, 3) + dopCoeff[2] * math.pow(ii, 2) + dopCoeff[1] * math.pow(ii, 1) + dopCoeff[0]
                testi *= prf
                if ii % 1000 == 0:
                    print("%7d %8.2f %8.2f %8.2f %8.2f %8.2f" % (ii, dc[0][ii], dc[1][ii], dc[2][ii], dcAverage[ii], testi))
        ######################################################################################################################
    except (IOError, OSError) as strerr:
        print("IOError: %s" % strerr)
        return None
    finally:
        if xmlfp is not None:
            xmlfp.close()
    return dopCoeff


def getAlos2Doppler(leader, frame):
    import math
    import numpy
    from isceobj.Constants import SPEED_OF_LIGHT

    with open(leader, 'rb') as f:
        #skip SAR Leader file descriptor records
        f.seek(720, 0)
        #skip front Data set summary records
        f.seek(1734, 1)
        a = float(f.read(16))
        b = float(f.read(16))
        print(a, b)
    
    width = frame.getNumberOfSamples()
    startingRange = frame.startingRange
    slantRangePixelSpacing = 0.5 * SPEED_OF_LIGHT / frame.rangeSamplingRate
    prf = frame.PRF

    dc = []
    for i in range(width):
        rangex = (startingRange + slantRangePixelSpacing * i) / 1000.0
        dc.append(a + b * rangex)

    rangeBins = list(range(width)) #[0:samples-1]
    p = numpy.polyfit(rangeBins, dc, 1)
    dopCoeff = []
    dopCoeff.append(p[1] / prf)
    dopCoeff.append(p[0] / prf)
    dopCoeff.append(0.)
    dopCoeff.append(0.)

    #just print out for debug
    ######################################################################################################################
    DEBUG = None
    if DEBUG != None:
        print('\n   rgbin     dc     dcfit')
        print('*****************************')
        for ii in range(width):
            testi = dopCoeff[3] * math.pow(ii, 3) + dopCoeff[2] * math.pow(ii, 2) + dopCoeff[1] * math.pow(ii, 1) + dopCoeff[0]
            testi *= prf
            if ii % 1000 == 0:
                print("%7d %8.2f %8.2f" % (ii, dc[ii], testi))
    ######################################################################################################################
    return dopCoeff






def make_raw(sensor, doppler):
    from make_raw import make_raw
    objMakeRaw = make_raw()
    objMakeRaw(sensor=sensor, doppler=doppler)
    return objMakeRaw

def setTSX(info):
    #set sensor
    objSensor = createSensor('TERRASARX')
    objSensor.xml = getInfo(info.input, 'xml')
    objSensor.output = getInfo(info.input, 'output')

    #set doppler
    objDop =  createDoppler('USEDOPTSX')
    #objDop.wireInputPort(name='sensor', object=objSensor)

    return (objSensor, objDop)


def setS1A(info):
    #set sensor
    objSensor = createSensor('SENTINEL1A')
    objSensor.xml = getInfo(info.input, 'xml')
    objSensor.tiff = getInfo(info.input, 'tiff')
    objSensor.orbitfile = getInfo(info.input, 'orbitfile')
    objSensor.output = getInfo(info.input, 'output')
    
    #set doppler
    objDop =  createDoppler('useDEFAULT')
    #objDop.wireInputPort(name='sensor', object=objSensor)

    return (objSensor, objDop)


def setCSKS(info):
    #set sensor
    objSensor = createSensor('COSMO_SKYMED_SLC')
    objSensor.hdf5 = getInfo(info.input, 'HDF5')
    objSensor.output = getInfo(info.input, 'output')

    #set doppler
    #objDop =  createDoppler('useDOPCSKSLC')
    objDop =  createDoppler('useDEFAULT')
    #objDop.wireInputPort(name='sensor', object=objSensor)

    return (objSensor, objDop)

def setALOS2(info):
    #set sensor
    objSensor = createSensor('ALOS2')
    objSensor._leaderFile = getInfo(info.input, 'LEADERFILE')
    objSensor._imageFile = getInfo(info.input, 'IMAGEFILE')
    objSensor.output = getInfo(info.input, 'output')

    #set doppler
    objDop =  createDoppler('useDEFAULT')
    #objDop.wireInputPort(name='sensor', object=objSensor)

    return (objSensor, objDop)

def setRS2(info):
    #set sensor
    objSensor = createSensor('RADARSAT2')
    objSensor.xml = getInfo(info.input, 'xml')
    objSensor.tiff = getInfo(info.input, 'tiff')
    objSensor.output = getInfo(info.input, 'output')
    
    #set doppler
    objDop =  createDoppler('useDEFAULT')
    #objDop.wireInputPort(name='sensor', object=objSensor)

    return (objSensor, objDop)


if __name__ == '__main__':
    '''
    Main driver.
    '''

    inps = cmdLineParse()
    
    if inps.sensor.upper() == 'TERRASARX':
        (objSensor, objDop) = setTSX(inps)
    elif inps.sensor.upper() == 'SENTINEL1A':
        (objSensor, objDop) = setS1A(inps)
    elif inps.sensor.upper() == 'COSMO_SKYMED_SLC':
        (objSensor, objDop) = setCSKS(inps)
    elif inps.sensor.upper() == 'ALOS2':
        (objSensor, objDop) = setALOS2(inps)
    elif inps.sensor.upper() == 'RADARSAT2':
        (objSensor, objDop) = setRS2(inps)
    else:
        raise Exception("\ninstrument not supported yet!\n")

    #extract image
    sarMeta = make_raw(objSensor, objDop)

    #get info
    frame = sarMeta.frame
    frame.dopCoeff = sarMeta.getDopplerValues().getDopplerCoefficients(inHz=False) #make doppler simple list

    #get doppler for s1a
    if inps.sensor.upper() == 'SENTINEL1A':
        #get doppler
        dopCoeff = getS1aDoppler(getInfo(inps.input, 'xml'))
        frame.dopCoeff = dopCoeff

    #update doppler for alos2
    if inps.sensor.upper() == 'ALOS2':
        dopCoeff = getAlos2Doppler(getInfo(inps.input, 'LEADERFILE'), frame)
        frame.dopCoeff = dopCoeff

    #change some information for RADARSAT2
    #if inps.sensor.upper() == 'RADARSAT2':
        #make PRF two times of its original value, still don't know why
        #frame.instrument.setPulseRepetitionFrequency(frame.PRF * 2.0)
        #change doppler frequency coefficients accordingly
        #frame.dopCoeff = [x/2.0 for x in frame.dopCoeff]


    #output information
    pickName = getInfo(inps.input, 'output') + '.pck'
    with open(pickName, 'wb') as f:
        pickle.dump(frame, f)
    frame.getImage().renderHdr()

    #calculate fmrate
    cmd = '$INSAR_ZERODOP_SCR/fmrate.py -f {}'.format(pickName)
    runCmd(cmd)

    #calculate burst start times
    if hasattr(frame, 'nbraw') and hasattr(frame, 'ncraw'):
        #cmd = '$INSAR_ZERODOP_SCR/burst_time.py -s {} -f {} -l {} -c {}'.format(getInfo(inps.input, 'output'), pickName, int(frame.getNumberOfLines()/2.0), 1000)
        cmd = '$INSAR_ZERODOP_SCR/burst_time2.py -s {} -f {}'.format(getInfo(inps.input, 'output'), pickName)
        runCmd(cmd)





#./makeData.py -s SENTINEL1A -i 20140807.xml
#./makeData.py -s TERRASARX -i 20130927_tsx.xml
#./makeData.py -s COSMO_SKYMED_SLC -i 20090412.xml














#sensor list
# SENSORS = {'ALOS' : createALOS,
#            'ALOS2' : createALOS2,
#            'COSMO_SKYMED' : createCOSMO_SkyMed,
#            'COSMO_SKYMED_SLC' : createCOSMO_SkyMed_SLC,
#            'ENVISAT' : createEnviSAT,
#            'ERS' : createERS,
#            'KOMPSAT5' : createKOMPSAT5,
#            'RADARSAT1' : createRadarsat1,
#            'RADARSAT2' : createRadarsat2,
#            'ROI_PAC' : createROI_PAC,
#            'TERRASARX' : createTerraSARX,
#            'RISAT1' : createRisat1,
#            'UAVSAR_RPI' : createUAVSAR_RPI,
#            'UAVSAR_STACK' : createUAVSAR_Stack,
#            'SENTINEL1A' : createSentinel1A}

