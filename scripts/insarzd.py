#!/usr/bin/env python3

#Cunren Liang, Apr. 21, 2015
#JPL/Caltech
#update, Apr. 2017

import os
import sys
import pickle
import argparse
import datetime
from xml.etree.ElementTree import ElementTree

import isce
import isceobj

from crlpac import Dummy
from crlpac import get_content
from crlpac import run_record_cmd


def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="InSAR processing in zero Doppler geometry\nPROCESSING STEPS: read_data, coreg_rdr2topo, coreg_topo2rdr, \
resamp, interf, coreg_ampsim, rect, diff, look, coherence, filter, unwrap, geocode\n\
                    ")
    parser.add_argument('-i', '--input', dest='input', type=str, required=True,
            help = 'input xml configure file')
    parser.add_argument('-s','--start', dest='start', type=str, default='read_data',
            help = 'starting step')
    parser.add_argument('-e','--end', dest='end', type=str, default='geocode',
            help = 'ending step')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


def read_insar_arg(xmlfile):
    
    #################################
    #define processing settings
    #################################
    insar_arg = Dummy()

    #mandatory
    insar_arg.sensor = None #str
    insar_arg.masterXml = None # str
    insar_arg.slaveXml = None # str
    insar_arg.dem = None #str DEM dir and its name, this DEM is used for geometrical calculations, except geocoding
    insar_arg.demGeo = None #str DEM dir and its name, this DEM is used for geocoding

    #optional
    insar_arg.ampcor = 'False' #str
    insar_arg.nrlks = 16  #int
    insar_arg.nalks = 16  #int
    insar_arg.nrlksMatch = insar_arg.nrlks  #int
    insar_arg.nalksMatch = insar_arg.nalks  #int
    insar_arg.filterStrength = 0.5 #float
    insar_arg.bbox = 'useMasterDefaultBbox' #str
    insar_arg.deleteFlag = 'False' #str


    #################################
    #read processing settings
    #################################
    #mandatory
    insar_arg.sensor = get_content(xmlfile, 'component.property.value', 'insar.sensor name.')
    if insar_arg.sensor == None:
        raise Exception('sensor name not set!')
    
    insar_arg.masterXml = get_content(xmlfile, 'component.component.catalog', 'insar.master.')
    if insar_arg.masterXml == None:
        raise Exception('master not set!')

    insar_arg.slaveXml = get_content(xmlfile, 'component.component.catalog', 'insar.slave.')
    if insar_arg.slaveXml == None:
        raise Exception('slave not set!')

    insar_arg.dem = get_content(xmlfile, 'component.component.catalog', 'insar.dem.')
    if insar_arg.dem == None:
        raise Exception('dem not set!')

    insar_arg.demGeo = get_content(xmlfile, 'component.component.catalog', 'insar.dem for geocoding.')
    if insar_arg.demGeo == None:
        raise Exception('dem for geocoding not set!')


    #optional
    optionalParameter = get_content(xmlfile, 'component.property.value', 'insar.whether do ampcor matching.')
    if optionalParameter != None:
        insar_arg.ampcor = optionalParameter

    optionalParameter = get_content(xmlfile, 'component.property.value', 'insar.number of range looks.')
    if optionalParameter != None:
        insar_arg.nrlks = int(optionalParameter)

    optionalParameter = get_content(xmlfile, 'component.property.value', 'insar.number of azimuth looks.')
    if optionalParameter != None:
        insar_arg.nalks = int(optionalParameter)

    optionalParameter = get_content(xmlfile, 'component.property.value', 'insar.number of range looks for matching between radar and simulated amplitudes.')
    if optionalParameter != None:
        insar_arg.nrlksMatch = int(optionalParameter)

    optionalParameter = get_content(xmlfile, 'component.property.value', 'insar.number of azimuth looks for matching between radar and simulated amplitudes.')
    if optionalParameter != None:
        insar_arg.nalksMatch = int(optionalParameter)
    
    optionalParameter = get_content(xmlfile, 'component.property.value', 'insar.filter strength.')
    if optionalParameter != None:
        insar_arg.filterStrength = float(optionalParameter)

    optionalParameter = get_content(xmlfile, 'component.property.value', 'insar.geocode bounding box.')
    if optionalParameter != None:
        insar_arg.bbox = optionalParameter.strip().strip('[').strip(']')
        #re-format bbox
        insar_arg.bbox = [float(val) for val in insar_arg.bbox.split(',')]
        if len(insar_arg.bbox) != 4:
            raise Exception('geocode bounding box should contain 4 floating point values')
        insar_arg.bbox = '{}/{}/{}/{}'.format(insar_arg.bbox[0], insar_arg.bbox[1], insar_arg.bbox[2], insar_arg.bbox[3])

    optionalParameter = get_content(xmlfile, 'component.property.value', 'insar.whether delete unnecessary large files.')
    if optionalParameter != None:
        insar_arg.deleteFlag = optionalParameter
    

    #################################
    #print processing settings
    #################################
    print('\ninput parameters')
    print('=====================================')
    for attr, value in vars(insar_arg).items():
        print("{}: {}".format(attr, value))
    print('=====================================\n')


    return insar_arg


def set_filename(insar_arg):

    fn = Dummy()

    #basic names
    fn.master = os.path.splitext((os.path.basename(get_content(insar_arg.masterXml, 'property.value', 'OUTPUT.'))))[0]
    if fn.master == None:
        raise Exception('master output not set!')
    fn.slave = os.path.splitext((os.path.basename(get_content(insar_arg.slaveXml, 'property.value', 'OUTPUT.'))))[0]
    if fn.slave == None:
        raise Exception('slave output not set!')
    ms = fn.master + '-' + fn.slave
    ml = '_{}rlks_{}alks'.format(insar_arg.nrlks, insar_arg.nalks) #do not use str(), which will change the variable type.

    #names from settings
    fn.dem = insar_arg.dem
    fn.demGeo = insar_arg.demGeo

    #output file names
    fn.masterSlc = fn.master + '.slc'
    fn.slaveSlc = fn.slave + '.slc'
    fn.resampleSlaveSlc = fn.slave + '_resamp.slc'
    fn.masterFrame = fn.master + '.slc.pck'
    fn.slaveFrame = fn.slave + '.slc.pck'
    fn.masterParameter = fn.master + '.slc.par.xml'
    fn.slaveParameter = fn.slave + '.slc.par.xml'

    fn.basline = ms + '_baseline' + '.xml'
    fn.interferogram = ms + '.int'
    fn.amplitude = ms + '.amp'
    fn.coherence = ms + '.cor' #actually not used, as we will calculate coherence after taking looks
    fn.differentialInterferogram = 'diff_' + ms + '.int'
    fn.multilookDifferentialInterferogram = 'diff_' + ms + ml + '.int'
    fn.multilookAmplitude = ms + ml + '.amp'
    fn.multilookCoherence = ms + ml + '.cor'
    fn.multilookPhsig = ms + ml + '.phsig'

    fn.filteredInterferogram = 'filt_' + 'diff_' + ms + ml + '.int'
    fn.unwrappedInterferogram =  'filt_' + 'diff_' + ms + ml + '.unw'
    fn.unwrappedMaskedInterferogram =  'filt_' + 'diff_' + ms + ml + '_msk.unw'

    #auxiliary
    fn.latitude = ms + '.lat'
    fn.longitude = ms + '.lon'
    fn.height = ms + '.hgt'
    fn.los = ms + '.los'
    fn.sim = ms + '.sim'
    fn.msk = ms + '.msk'
    fn.inc = ms + '.inc'
    fn.rangeOffset = ms + '_rg.off'
    fn.azimuthOffset = ms + '_az.off'
    fn.multilookLos = ms + ml + '.los'
    fn.multilookMsk = ms + ml + '.msk'
    fn.multilookLatitude = ms + ml + '.lat'
    fn.multilookLongitude = ms + ml + '.lon'
    fn.multilookHeight = ms + ml + '.hgt'

    #for matching
    fn.affineTransformation = 'ampsim.aff'
    fn.rectRangeOffset = 'rect_' + ms + '_rg.off'

    #geo
    fn.geoUnwrappedInterferogram = fn.unwrappedInterferogram + '.geo'
    fn.geoUnwrappedMaskedInterferogram = fn.unwrappedMaskedInterferogram + '.geo'
    fn.geoMultilookCoherence = fn.multilookCoherence + '.geo'
    fn.geoMultilookLos = fn.multilookLos + '.geo'


    #################################
    #print file names
    #################################
    print('\noutput file names')
    print('=====================================')
    for attr, value in vars(fn).items():
        print("{}: {}".format(attr, value))
    print('=====================================\n')


    return fn


def get_cmd(insar_arg, fn):

    #set your steps here
    cmd_all = [['read_data'], ['coreg_rdr2topo'], ['coreg_topo2rdr'], ['resamp'], ['interf'], ['coreg_ampsim'], ['rect'], ['diff'], ['look'], ['coherence'], ['filter'], ['unwrap'], ['geocode']]


#############################################################################
#   STEP. read_data
#############################################################################
    step_i = 0

    cmd = "$INSAR_ZERODOP_SCR/readData.py -s {} -i {}".format(
        insar_arg.sensor, 
        insar_arg.masterXml)
    cmd_all[step_i].append(cmd)

    cmd = "$INSAR_ZERODOP_SCR/readData.py -s {} -i {}".format(
        insar_arg.sensor, 
        insar_arg.slaveXml)
    cmd_all[step_i].append(cmd)

    #create parameter file
    cmd = "$INSAR_ZERODOP_SCR/create_parxml.py -frame {} -par {}".format(fn.masterFrame, fn.masterParameter)
    cmd_all[step_i].append(cmd)
    cmd = "$INSAR_ZERODOP_SCR/create_parxml.py -frame {} -par {}".format(fn.slaveFrame, fn.slaveParameter)
    cmd_all[step_i].append(cmd)

    cmd = "$INSAR_ZERODOP_SCR/calBaseline.py -m {} -s {} -o {}".format(
        fn.masterFrame, 
        fn.slaveFrame, 
        fn.basline)
    cmd_all[step_i].append(cmd)


#############################################################################
#   STEP. coreg_rdr2topo
#############################################################################
    step_i += 1

    cmd = "$INSAR_ZERODOP_SCR/topo.py -m {} -d {} -a {} -o {} -z {} -l {} -s {} -i {} -k {} -rlks 1 -alks 1".format(
        fn.masterFrame, 
        fn.dem, 
        fn.latitude, 
        fn.longitude, 
        fn.height, 
        fn.los, 
        fn.sim, 
        fn.inc, 
        fn.msk)
    cmd_all[step_i].append(cmd)


#############################################################################
#   STEP. coreg_topo2rdr
#############################################################################
    step_i += 1

    cmd = "$INSAR_ZERODOP_SCR/geo2rdr.py -s {} -a {} -o {} -z {} -r {} -i {} -rlks 1 -alks 1".format(
        fn.slaveFrame, 
        fn.latitude, 
        fn.longitude, 
        fn.height, 
        fn.rangeOffset, 
        fn.azimuthOffset)
    cmd_all[step_i].append(cmd)


#############################################################################
#   STEP. resamp
#############################################################################
    step_i += 1

    if insar_arg.ampcor.upper() == 'FALSE':
        cmd = "$INSAR_ZERODOP_SCR/resample.py -m {} -s {} -r {} -a {} -o {}".format(
            fn.masterSlc, 
            fn.slaveSlc, 
            fn.rangeOffset, 
            fn.azimuthOffset, 
            fn.resampleSlaveSlc)
        cmd_all[step_i].append(cmd)
    else:
        cmd = "$INSAR_ZERODOP_SCR/match_resamp_isce.py -m {} -s {} -int {} -amp {} -nr 30 -na 30 -rfw 64 -afw 64 -rsm 32 -asm 32 -rlks 1 -alks 1".format(
            fn.masterSlc, 
            fn.slaveSlc, 
            fn.interferogram, 
            fn.amplitude
            )
        cmd_all[step_i].append(cmd)


#############################################################################
#   STEP. interf
#############################################################################
    step_i += 1

    if insar_arg.ampcor.upper() == 'FALSE':
        cmd = "$INSAR_ZERODOP_SCR/interf.py -m {} -s {} -i {} -a {}".format(
            fn.masterSlc, 
            fn.resampleSlaveSlc, 
            fn.interferogram, 
            fn.amplitude)
        cmd_all[step_i].append(cmd)


#############################################################################
#   STEP. coreg_ampsim
#############################################################################
    step_i += 1
    cmd = "$INSAR_ZERODOP_SCR/match_ampsim.py -m {} -s {} -f {} -r {} -a {}".format(
        fn.amplitude, 
        fn.sim, 
        fn.affineTransformation, 
        insar_arg.nrlksMatch, 
        insar_arg.nalksMatch)
    cmd_all[step_i].append(cmd)


#############################################################################
#   STEP. rect
#############################################################################
    step_i += 1
    cmd = "$INSAR_ZERODOP_SCR/rect.py -i {} -o {} -t {} -f {} -r {} -a {} -r0 1 -a0 1".format(
        fn.rangeOffset, 
        fn.rectRangeOffset, 
        fn.affineTransformation, 
        'real', 
        insar_arg.nrlksMatch, 
        insar_arg.nalksMatch)
    cmd_all[step_i].append(cmd)


#############################################################################
#   STEP. diff
#############################################################################
    step_i += 1
    cmd = "$INSAR_ZERODOP_SCR/flat.py -m {} -i {} -r {} -d {} -rlks 1".format(
        fn.masterFrame, 
        fn.interferogram, 
        fn.rectRangeOffset, 
        fn.differentialInterferogram)
    cmd_all[step_i].append(cmd)


#############################################################################
#   STEP. look
#############################################################################
    step_i += 1
    if insar_arg.nrlks != 1 or insar_arg.nalks != 1:
        cmd = "$INSAR_ZERODOP_SCR/look.py -i {} -o {} -r {} -a {}".format(fn.differentialInterferogram, fn.multilookDifferentialInterferogram, insar_arg.nrlks, insar_arg.nalks) 
        cmd_all[step_i].append(cmd)

        cmd = "$INSAR_ZERODOP_SCR/look.py -i {} -o {} -r {} -a {}".format(fn.amplitude, fn.multilookAmplitude, insar_arg.nrlks, insar_arg.nalks)
        cmd_all[step_i].append(cmd)

        cmd = "$INSAR_ZERODOP_SCR/look.py -i {} -o {} -r {} -a {}".format(fn.msk, fn.multilookMsk, insar_arg.nrlks, insar_arg.nalks)
        cmd_all[step_i].append(cmd)

        cmd = "$INSAR_ZERODOP_SCR/look.py -i {} -o {} -r {} -a {}".format(fn.latitude, fn.multilookLatitude, insar_arg.nrlks, insar_arg.nalks)
        cmd_all[step_i].append(cmd)

        cmd = "$INSAR_ZERODOP_SCR/look.py -i {} -o {} -r {} -a {}".format(fn.longitude, fn.multilookLongitude, insar_arg.nrlks, insar_arg.nalks)
        cmd_all[step_i].append(cmd)

        cmd = "$INSAR_ZERODOP_SCR/look.py -i {} -o {} -r {} -a {}".format(fn.height, fn.multilookHeight, insar_arg.nrlks, insar_arg.nalks)
        cmd_all[step_i].append(cmd)

        #using isce's look program, tested correctness
        #tested file types: 1. BIP, 2 bands
        #                   2. BIL, 2 bands (the following case)
        cmd = "looks.py -i {} -o {} -r {} -a {}".format(fn.los, fn.multilookLos, insar_arg.nrlks, insar_arg.nalks)
        cmd_all[step_i].append(cmd)

    else:
        cmd = "cp {} {}".format(fn.differentialInterferogram, fn.multilookDifferentialInterferogram)
        cmd_all[step_i].append(cmd)
        cmd = "cp {} {}".format(fn.differentialInterferogram + '.xml', fn.multilookDifferentialInterferogram + '.xml')
        cmd_all[step_i].append(cmd) 
        cmd = "cp {} {}".format(fn.differentialInterferogram + '.vrt', fn.multilookDifferentialInterferogram + '.vrt')
        cmd_all[step_i].append(cmd) 

        cmd = "cp {} {}".format(fn.amplitude, fn.multilookAmplitude)
        cmd_all[step_i].append(cmd)
        cmd = "cp {} {}".format(fn.amplitude + '.xml', fn.multilookAmplitude + '.xml')
        cmd_all[step_i].append(cmd)
        cmd = "cp {} {}".format(fn.amplitude + '.vrt', fn.multilookAmplitude + '.vrt')
        cmd_all[step_i].append(cmd)
             
        cmd = "cp {} {}".format(fn.msk, fn.multilookMsk)
        cmd_all[step_i].append(cmd)
        cmd = "cp {} {}".format(fn.msk + '.xml', fn.multilookMsk + '.xml')
        cmd_all[step_i].append(cmd)
        cmd = "cp {} {}".format(fn.msk + '.vrt', fn.multilookMsk + '.vrt')
        cmd_all[step_i].append(cmd)

        cmd = "cp {} {}".format(fn.latitude, fn.multilookLatitude)
        cmd_all[step_i].append(cmd)
        cmd = "cp {} {}".format(fn.latitude + '.xml', fn.multilookLatitude + '.xml')
        cmd_all[step_i].append(cmd)
        cmd = "cp {} {}".format(fn.latitude + '.vrt', fn.multilookLatitude + '.vrt')
        cmd_all[step_i].append(cmd)

        cmd = "cp {} {}".format(fn.longitude, fn.multilookLongitude)
        cmd_all[step_i].append(cmd)
        cmd = "cp {} {}".format(fn.longitude + '.xml', fn.multilookLongitude + '.xml')
        cmd_all[step_i].append(cmd)
        cmd = "cp {} {}".format(fn.longitude + '.vrt', fn.multilookLongitude + '.vrt')
        cmd_all[step_i].append(cmd)

        cmd = "cp {} {}".format(fn.height, fn.multilookHeight)
        cmd_all[step_i].append(cmd)
        cmd = "cp {} {}".format(fn.height + '.xml', fn.multilookHeight + '.xml')
        cmd_all[step_i].append(cmd)
        cmd = "cp {} {}".format(fn.height + '.vrt', fn.multilookHeight + '.vrt')
        cmd_all[step_i].append(cmd)

        cmd = "cp {} {}".format(fn.los, fn.multilookLos)
        cmd_all[step_i].append(cmd)
        cmd = "cp {} {}".format(fn.los + '.xml', fn.multilookLos + '.xml')
        cmd_all[step_i].append(cmd)
        cmd = "cp {} {}".format(fn.los + '.vrt', fn.multilookLos + '.vrt')
        cmd_all[step_i].append(cmd)


#############################################################################
#   STEP. coherence
#############################################################################
    step_i += 1

    #method phase_gradient does not work, there must be bug in it.
    cmd = "$INSAR_ZERODOP_SCR/coherence.py -i {} -a {} -c {}".format(
        fn.multilookDifferentialInterferogram, 
        fn.multilookAmplitude, 
        fn.multilookCoherence)
    cmd_all[step_i].append(cmd)


#############################################################################
#   STEP. filter
#############################################################################
    step_i += 1

    cmd = "$INSAR_ZERODOP_SCR/filter.py -i {} -m {} -k {} -f {} -p {} -a {}".format(
        fn.multilookDifferentialInterferogram, 
        fn.multilookAmplitude, 
        fn.multilookMsk, 
        fn.filteredInterferogram, 
        fn.multilookPhsig, 
        insar_arg.filterStrength)
    cmd_all[step_i].append(cmd)


#############################################################################
#   STEP. unwrap
#############################################################################
    step_i += 1

    cmd = "$INSAR_ZERODOP_SCR/snaphuMCF.py -mframe {} -inf {} -cor {} -unw {} -munw {} -rlks {} -alks {}".format(
        fn.masterFrame, 
        fn.filteredInterferogram, 
        fn.multilookPhsig, 
        fn.unwrappedInterferogram, 
        fn.unwrappedMaskedInterferogram, 
        insar_arg.nrlks, 
        insar_arg.nalks)
    cmd_all[step_i].append(cmd)

    cmd = "$INSAR_ZERODOP_SCR/replace_mag.py -unw {} -inf {} -msk 0".format(
        fn.unwrappedMaskedInterferogram, 
        fn.multilookDifferentialInterferogram)
    cmd_all[step_i].append(cmd)

    cmd = "$INSAR_ZERODOP_SCR/replace_mag.py -unw {} -inf {} -msk 0".format(
        fn.unwrappedInterferogram, 
        fn.multilookDifferentialInterferogram)
    cmd_all[step_i].append(cmd)


#############################################################################
#   STEP. geocode
#############################################################################
    step_i += 1

    # geocode considering offsets between radar amplitude and simulated amplitude
    #cmd = "$INSAR_ZERODOP_SCR/geo_off.py -m {} -d {} -i {} -o {} -b {} -r {} -a {} -t {} -x {} -y {}".format(fn.masterFrame, fn.demGeo, fn.unwrappedInterferogram, fn.geoUnwrappedInterferogram, insar_arg.bbox, insar_arg.nrlks, insar_arg.nalks, fn.affineTransformation, insar_arg.nrlksMatch, insar_arg.nalksMatch)
    cmd = "$INSAR_ZERODOP_SCR/geo.py -m {} -d {} -i {} -o {} -b {} -r {} -a {}".format(fn.masterFrame, fn.demGeo, fn.unwrappedInterferogram, fn.geoUnwrappedInterferogram, insar_arg.bbox, insar_arg.nrlks, insar_arg.nalks)
    cmd_all[step_i].append(cmd)

    #cmd = "$INSAR_ZERODOP_SCR/geo_off.py -m {} -d {} -i {} -o {} -b {} -r {} -a {} -t {} -x {} -y {}".format(fn.masterFrame, fn.demGeo, fn.unwrappedMaskedInterferogram, fn.geoUnwrappedMaskedInterferogram, insar_arg.bbox, insar_arg.nrlks, insar_arg.nalks, fn.affineTransformation, insar_arg.nrlksMatch, insar_arg.nalksMatch)
    cmd = "$INSAR_ZERODOP_SCR/geo.py -m {} -d {} -i {} -o {} -b {} -r {} -a {}".format(fn.masterFrame, fn.demGeo, fn.unwrappedMaskedInterferogram, fn.geoUnwrappedMaskedInterferogram, insar_arg.bbox, insar_arg.nrlks, insar_arg.nalks)
    cmd_all[step_i].append(cmd)

    #cmd = "$INSAR_ZERODOP_SCR/geo_off.py -m {} -d {} -i {} -o {} -b {} -r {} -a {} -t {} -x {} -y {}".format(fn.masterFrame, fn.demGeo, fn.multilookCoherence, fn.geoMultilookCoherence, insar_arg.bbox, insar_arg.nrlks, insar_arg.nalks, fn.affineTransformation, insar_arg.nrlksMatch, insar_arg.nalksMatch)
    cmd = "$INSAR_ZERODOP_SCR/geo.py -m {} -d {} -i {} -o {} -b {} -r {} -a {}".format(fn.masterFrame, fn.demGeo, fn.multilookCoherence, fn.geoMultilookCoherence, insar_arg.bbox, insar_arg.nrlks, insar_arg.nalks)
    cmd_all[step_i].append(cmd)

    #cmd = "$INSAR_ZERODOP_SCR/geo_off.py -m {} -d {} -i {} -o {} -b {} -r {} -a {} -t {} -x {} -y {}".format(fn.masterFrame, fn.demGeo, fn.multilookLos, fn.geoMultilookLos, insar_arg.bbox, insar_arg.nrlks, insar_arg.nalks, fn.affineTransformation, insar_arg.nrlksMatch, insar_arg.nalksMatch)
    cmd = "$INSAR_ZERODOP_SCR/geo.py -m {} -d {} -i {} -o {} -b {} -r {} -a {}".format(fn.masterFrame, fn.demGeo, fn.multilookLos, fn.geoMultilookLos, insar_arg.bbox, insar_arg.nrlks, insar_arg.nalks)
    cmd_all[step_i].append(cmd)

    return cmd_all


def del_file(insar_arg, fn, end):
    
    #users can define their own files to be deleted.
    uselessFiles=(fn.amplitude,
                 fn.latitude,
                 fn.longitude,
                 fn.height,
                 fn.sim)

    if end == 'geocode' and insar_arg.deleteFlag.upper() == 'TRUE':
        print('deleting useless files')
        for filex in uselessFiles:
            os.remove(filex)


if __name__ == '__main__':
    '''
    Main driver.
    '''

    inps = cmdLineParse()

    insar_arg = read_insar_arg(inps.input)
    fn = set_filename(insar_arg)
    cmd_all = get_cmd(insar_arg, fn)

    run_record_cmd(cmd_all, inps.start, inps.end, os.path.splitext(os.path.basename(__file__))[0] + '.log')

    del_file(insar_arg, fn, inps.end)






