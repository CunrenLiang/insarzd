#!/usr/bin/env python3

#Cunren Liang, Apr. 21, 2015
#JPL/Caltech
#update, Apr. 2017


import os
import sys
import glob
import argparse

from crlpac import run_record_cmd
from crlpac import writeSARconfig_ALOS2
from crlpac import read_insar_arg
from crlpac import set_filename


def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="Interferometric Processing of ALOS-2 Data\nPROCESSING STEPS: read_data, prep_slc, form_interf, mosaic_subswath, mosaic_frame, \ncoreg_rdr2topo, coreg_topo2rdr, coreg_ampsim, rect, diff, look, coherence, filter, \nunwrap, geocode\n\
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


def prep_inputfile(insar_arg, fn):
#prepare input files

    for i, (frame_m, frame_s) in enumerate(zip(insar_arg.masterFrames, insar_arg.slaveFrames)):
        frame_dir = 'f{}'.format(frame_m)
        os.mkdir(frame_dir)
        os.chdir(frame_dir)

################################################################################################
# there may be same frames, we add sequence number to frame. 01-MAY-2017
# only three files are changed due to this: funcs.py, alos2app.py, alos2app_burst.py
# frames are used to do the following two things:
# 1. find original data files: if insar_arg.masterFrames and insar_arg.slaveFrames are changed
#                              the associated things need to be changed accordingly
# 2. handle frame directories: if insar_arg.masterFrames and insar_arg.slaveFrames are changed
#                              the associated things do not need to be changed, as frame
#                              directory can use any name.
        frame_m = frame_m.split('_')[1]
        frame_s = frame_s.split('_')[1]
################################################################################################

        if insar_arg.combination == 0:
            ldr_m = sorted(glob.glob(insar_arg.masterDir + "/" + "LED-ALOS2*{}-*-*".format(frame_m)))[0]
            img_m = sorted(glob.glob(insar_arg.masterDir + "/" + "IMG-{}-ALOS2*{}-*-*".format(insar_arg.masterPolarization, frame_m)))[0]
            writeSARconfig_ALOS2(fn.master + '.xml', ldr_m, img_m, fn.master + '.slc')

            ldr_s = sorted(glob.glob(insar_arg.slaveDir + "/" + "LED-ALOS2*{}-*-*".format(frame_s)))[0]
            img_s = sorted(glob.glob(insar_arg.slaveDir + "/" + "IMG-{}-ALOS2*{}-*-*".format(insar_arg.slavePolarization, frame_s)))[0]
            writeSARconfig_ALOS2(fn.slave + '.xml', ldr_s, img_s, fn.slave + '.slc')
        else:
            for j in range(insar_arg.startingSubswath, insar_arg.endingSubswath+1):
                subswath_dir = 's{}'.format(j)
                os.mkdir(subswath_dir)
                os.chdir(subswath_dir)

                ldr_m = sorted(glob.glob(insar_arg.masterDir + "/" + "LED-ALOS2*{}-*-*".format(frame_m)))[0]
                img_m = sorted(glob.glob(insar_arg.masterDir + "/" + "IMG-{}-ALOS2*{}-*-*-F{}".format(insar_arg.masterPolarization, frame_m, j)))[0]
                writeSARconfig_ALOS2(fn.master + '.xml', ldr_m, img_m, fn.master + '.slc')    


                ldr_s = sorted(glob.glob(insar_arg.slaveDir + "/" + "LED-ALOS2*{}-*-*".format(frame_s)))[0]
                if insar_arg.combination == 1:
                    img_s = sorted(glob.glob(insar_arg.slaveDir + "/" + "IMG-{}-ALOS2*{}-*-*".format(insar_arg.slavePolarization, frame_s)))[0]
                else:
                    img_s = sorted(glob.glob(insar_arg.slaveDir + "/" + "IMG-{}-ALOS2*{}-*-*-F{}".format(insar_arg.slavePolarization, frame_s, j)))[0]
                writeSARconfig_ALOS2(fn.slave + '.xml', ldr_s, img_s, fn.slave + '.slc')

                os.chdir('../')
        os.chdir('../')


def get_cmd(insar_arg, fn):

    #subswath mosaicking directory
    frame_mosaic_dir = 'mosaic'
    #final insar processing directory
    insar_dir = 'insar'

    #set your steps here
    cmd_all = [['read_data'], ['prep_slc'], ['form_interf'], ['mosaic_subswath'], ['mosaic_frame'], ['coreg_rdr2topo'], ['coreg_topo2rdr'], ['coreg_ampsim'], ['rect'], ['diff'], ['look'], ['coherence'], ['filter'], ['unwrap'], ['geocode']]

#############################################################################
#   STEP. read_data
#############################################################################
    step_i = 0

    #commands
    cmd1 = "$INSAR_ZERODOP_SCR/readData.py -s {} -i {}".format(
        insar_arg.sensor, 
        fn.master + '.xml')
    cmd2 = "$INSAR_ZERODOP_SCR/readData.py -s {} -i {}".format(
        insar_arg.sensor, 
        fn.slave + '.xml')
    #create parameter file
    cmd3 = "$INSAR_ZERODOP_SCR/create_parxml.py -frame {} -par {}".format(fn.masterFrame, fn.masterParameter)
    cmd4 = "$INSAR_ZERODOP_SCR/create_parxml.py -frame {} -par {}".format(fn.slaveFrame, fn.slaveParameter)

    #loop over frame and subswath directories
    for frame in insar_arg.masterFrames:
        cmd_all[step_i].append("cd f{}".format(frame))
        if insar_arg.combination == 0:
            cmd_all[step_i].extend([cmd1, cmd2, cmd3, cmd4])
        else:
            for subswath in range(insar_arg.startingSubswath, insar_arg.endingSubswath+1):
                cmd_all[step_i].append("cd s{}".format(subswath))
                cmd_all[step_i].extend([cmd1, cmd2, cmd3, cmd4])
                cmd_all[step_i].append("cd ../")
        cmd_all[step_i].append("cd ../")


#############################################################################
#   STEP. prep_slc
#############################################################################
    step_i += 1

    #crop the huge ScanSAR image for ScanSAR-stripmap interferometry
    cmd0 = "$INSAR_ZERODOP_SCR/crop_slc.py -slc {} -frame {} -frame2 {}".format(
        fn.masterSlc,
        fn.masterFrame,
        fn.slaveFrame)

    cmd1 = "$INSAR_ZERODOP_SCR/range_filter.py -mslc {} -mframe {} -sslc {} -sframe {}".format(
        fn.masterSlc,
        fn.masterFrame,
        fn.slaveSlc,
        fn.slaveFrame)

    cmd2 = "$INSAR_ZERODOP_SCR/equalize_sample_size.py -slc {} -frame {} -frame2 {}".format(
        fn.slaveSlc, 
        fn.slaveFrame, 
        fn.masterFrame)

    cmd3 = "$INSAR_ZERODOP_SCR/burst_sync.py -m {} -s {} -d {}".format(
        fn.masterFrame,
        fn.slaveFrame,
        fn.dem)
    if insar_arg.combination == 1:
        cmd3 += " -c 1 -i 1"

    s1_sycnfile = '../s{}/burst_synchronization.txt'.format(insar_arg.startingSubswath)
    cmd4 = "$INSAR_ZERODOP_SCR/unsync_removal.py -mslc {} -mframe {} -sslc {} -sframe {} -syncf {} -syncf2 {}".format(
        fn.masterSlc,
        fn.masterFrame,
        fn.slaveSlc,
        fn.slaveFrame,
        'burst_synchronization.txt',
        s1_sycnfile)
    if insar_arg.combination == 1:
        #                                             trick for always doing mbf filtering
        #                                             for ScanSAR-stripmap interferometry
        cmd4 += " -synthr {} -mfilt 0 -sfilt 1".format(101.0)
    elif insar_arg.combination == 2:
        cmd4 += " -synthr {} -mfilt 1 -sfilt 1".format(insar_arg.burstOverlapThreshhold)

    #loop over frame and subswath directories
    for frame in insar_arg.masterFrames:
        cmd_all[step_i].append("cd f{}".format(frame))
        if insar_arg.combination == 0:
            cmd_all[step_i].extend([cmd1, cmd2])
        else:
            for subswath in range(insar_arg.startingSubswath, insar_arg.endingSubswath+1):
                cmd_all[step_i].append("cd s{}".format(subswath))
                #crop the huge ScanSAR image for ScanSAR-stripmap interferometry
                if insar_arg.combination == 1:
                    cmd_all[step_i].append(cmd0)
                cmd_all[step_i].extend([cmd1, cmd2, cmd3, cmd4])
                cmd_all[step_i].append("cd ../")
        cmd_all[step_i].append("cd ../")


#############################################################################
#   STEP. form_interf
#############################################################################
    step_i += 1

    if insar_arg.combination == 0:
        nrgoff = 30
        nazoff = 30
        rgfftw = 64
        azfftw = 64
        rgsmax = 32
        azsmax = 32
    else:
        nrgoff = 40
        nazoff = 100
        rgfftw = 64
        azfftw = 512
        rgsmax = 32
        azsmax = 32

    cmd = "$INSAR_ZERODOP_SCR/match_resamp_isce.py -m {} -s {} -int {} -amp {} -nr {} -na {} -rfw {} -afw {} -rsm {} -asm {} -rlks {} -alks {}".format(
        fn.masterSlc, 
        fn.slaveSlc, 
        fn.interferogram, 
        fn.amplitude, 
        nrgoff,
        nazoff,
        rgfftw,
        azfftw,
        rgsmax,
        azsmax,
        insar_arg.nrlks0,
        insar_arg.nalks0)

    #loop over frame and subswath directories
    for frame in insar_arg.masterFrames:
        cmd_all[step_i].append("cd f{}".format(frame))
        if insar_arg.combination == 0:
            cmd_all[step_i].append(cmd)
        else:
            for subswath in range(insar_arg.startingSubswath, insar_arg.endingSubswath+1):
                cmd_all[step_i].append("cd s{}".format(subswath))
                cmd_all[step_i].append(cmd)
                cmd_all[step_i].append("cd ../")
        cmd_all[step_i].append("cd ../")


#############################################################################
#   STEP. mosaic_subswath
#############################################################################
    step_i += 1

    #offset
    cmd1 = "$INSAR_ZERODOP_SCR/subswath_offset.py"
    for i in range(insar_arg.startingSubswath, insar_arg.endingSubswath+1):
        cmd1 += " -s ../s{}/{}".format(i, fn.masterSlc)
    for i in range(insar_arg.startingSubswath, insar_arg.endingSubswath+1):
        cmd1 += " -f ../s{}/{}".format(i, fn.masterFrame)
    cmd1 += " -o subswath_offset.txt -rlks 1 -alks 10"

    #mosaic
    cmd2 = "$INSAR_ZERODOP_SCR/subswath_mosaic.py"
    for i in range(insar_arg.startingSubswath, insar_arg.endingSubswath+1):
        cmd2 += " -f ../s{}/{}".format(i, fn.masterFrame)
    for i in range(insar_arg.startingSubswath, insar_arg.endingSubswath+1):
        cmd2 += " -i ../s{}/{}".format(i, fn.interferogram)
    for i in range(insar_arg.startingSubswath, insar_arg.endingSubswath+1):
        cmd2 += " -a ../s{}/{}".format(i, fn.amplitude)
    cmd2 += " -offset subswath_offset.txt -io {} -ao {} -fo {}".format(fn.interferogram, fn.amplitude, fn.masterFrame)
    cmd2 += " -imethod 0 -amethod 0 -rlks {} -alks {} -phc 0 -drlks 1 -dalks 4".format(insar_arg.nrlks0, insar_arg.nalks0)

    #slave frame pickle
    cmd3 = "$INSAR_ZERODOP_SCR/subswath_frame.py -s ../s{}/{} -e ../s{}/{} -o {}".format(insar_arg.startingSubswath, fn.slaveFrame, insar_arg.endingSubswath, fn.slaveFrame, fn.slaveFrame)

    #loop overlp frames
    ns = insar_arg.endingSubswath - insar_arg.startingSubswath + 1
    if (insar_arg.combination == 1 or insar_arg.combination == 2) and (ns > 1):
        for frame in insar_arg.masterFrames:
            cmd_all[step_i].append("cd f{}".format(frame))
            cmd_all[step_i].append("mkdir {}".format(frame_mosaic_dir))
            cmd_all[step_i].append("cd {}".format(frame_mosaic_dir))
            cmd_all[step_i].extend([cmd1, cmd2, cmd3])
            cmd_all[step_i].append("cd ../")
            cmd_all[step_i].append("cd ../")


#############################################################################
#   STEP. mosaic_frame
#############################################################################
    step_i += 1

    cmd_all[step_i].append("mkdir {}".format(insar_dir))
    cmd_all[step_i].append("cd {}".format(insar_dir))

    #processing according to different cases
    #####################################################################
    ns = insar_arg.endingSubswath - insar_arg.startingSubswath + 1
    if len(insar_arg.masterFrames) == 1:
        if insar_arg.combination == 0:
           dirsm = ''
        else:
            if ns == 1:
                dirsm = 's{}/'.format(insar_arg.startingSubswath)
            else:
                dirsm = '{}/'.format(frame_mosaic_dir)
        #copy small files and move big files
        cmd_all[step_i].append("cp ../f{}/{}{} ./".format(insar_arg.masterFrames[0], dirsm, fn.masterFrame))
        cmd_all[step_i].append("cp ../f{}/{}{} ./".format(insar_arg.masterFrames[0], dirsm, fn.slaveFrame))
        cmd_all[step_i].append("mv ../f{}/{}{} ./".format(insar_arg.masterFrames[0], dirsm, fn.interferogram))
        cmd_all[step_i].append("mv ../f{}/{}{} ./".format(insar_arg.masterFrames[0], dirsm, fn.amplitude))
        cmd_all[step_i].append("cp ../f{}/{}{} ./".format(insar_arg.masterFrames[0], dirsm, fn.interferogram+'.xml'))
        cmd_all[step_i].append("cp ../f{}/{}{} ./".format(insar_arg.masterFrames[0], dirsm, fn.amplitude+'.xml'))
        cmd_all[step_i].append("cp ../f{}/{}{} ./".format(insar_arg.masterFrames[0], dirsm, fn.interferogram+'.vrt'))
        cmd_all[step_i].append("cp ../f{}/{}{} ./".format(insar_arg.masterFrames[0], dirsm, fn.amplitude+'.vrt'))
    else:
        if insar_arg.combination == 0:
           dirsm = ''
           file_for_match = fn.masterSlc
           frame_for_match = fn.masterFrame
           cmd1_1 = "echo no need to update frame_offset.txt"
        else:
            if ns == 1:
                dirsm = 's{}/'.format(insar_arg.startingSubswath)
                file_for_match = fn.masterSlc
                frame_for_match = fn.masterFrame
                cmd1_1 = "echo no need to update frame_offset.txt"
            else:
                dirsm = '{}/'.format(frame_mosaic_dir)
                file_for_match = '../s{}/{}'.format(insar_arg.startingSubswath, fn.masterSlc)
                frame_for_match = '../s{}/{}'.format(insar_arg.startingSubswath, fn.masterFrame)

                #we use starting subswath slc for frame matching, but the azimuth starting time may be
                #changed in previous subswath mosaicking. we need to update frame offset here
                cmd1_1 = "$INSAR_ZERODOP_SCR/frame_offset_update.py -offset frame_offset.txt"
                for frame in insar_arg.masterFrames:
                    cmd1_1 += " -f ../f{}/s{}/{}".format(frame, insar_arg.startingSubswath, fn.masterFrame)
                for frame in insar_arg.masterFrames:
                    cmd1_1 += " -f2 ../f{}/{}/{}".format(frame, frame_mosaic_dir, fn.masterFrame)

        #offset
        cmd1 = "$INSAR_ZERODOP_SCR/frame_offset.py"
        for frame in insar_arg.masterFrames:
            cmd1 += " -s ../f{}/{}{}".format(frame, dirsm, file_for_match)
        for frame in insar_arg.masterFrames:
            cmd1 += " -f ../f{}/{}{}".format(frame, dirsm, frame_for_match)
        cmd1 += " -o frame_offset.txt -c 0"

        #mosaic amplitude
        cmd2 = "$INSAR_ZERODOP_SCR/frame_mosaic.py"
        for frame in insar_arg.masterFrames:
            cmd2 += " -f ../f{}/{}{}".format(frame, dirsm, fn.masterFrame)
        for frame in insar_arg.masterFrames:
            cmd2 += " -i ../f{}/{}{}".format(frame, dirsm, fn.amplitude)
        cmd2 += " -offset frame_offset.txt -io {} -fo {} -rlks {} -alks {} -phc 0".format(fn.amplitude, fn.masterFrame, insar_arg.nrlks0, insar_arg.nalks0)

        #mosaic interferogram
        cmd3 = "$INSAR_ZERODOP_SCR/frame_mosaic.py"
        for frame in insar_arg.masterFrames:
            cmd3 += " -f ../f{}/{}{}".format(frame, dirsm, fn.masterFrame)
        for frame in insar_arg.masterFrames:
            cmd3 += " -i ../f{}/{}{}".format(frame, dirsm, fn.interferogram)
        cmd3 += " -offset frame_offset.txt -io {} -fo {} -rlks {} -alks {} -phc 1".format(fn.interferogram, fn.masterFrame, insar_arg.nrlks0, insar_arg.nalks0)

        #slave frame pickle
        cmd4 = "$INSAR_ZERODOP_SCR/frame_frame.py"
        for frame in insar_arg.masterFrames:
            cmd4 += " -f ../f{}/{}{}".format(frame, dirsm, fn.slaveFrame)
        cmd4 += " -o {}".format(fn.slaveFrame)

        cmd_all[step_i].extend([cmd1, cmd1_1, cmd2, cmd3, cmd4])
    #####################################################################

    #create parameter file
    cmd_all[step_i].append("$INSAR_ZERODOP_SCR/create_parxml.py -frame {} -par {}".format(fn.masterFrame, fn.masterParameter))
    cmd_all[step_i].append("$INSAR_ZERODOP_SCR/create_parxml.py -frame {} -par {}".format(fn.slaveFrame, fn.slaveParameter))

    #calculate baseline
    cmd3 = "$INSAR_ZERODOP_SCR/calBaseline.py -m {} -s {} -o {}".format(
        fn.masterFrame, 
        fn.slaveFrame, 
        fn.basline)
    cmd_all[step_i].append(cmd3)

    #go back to original directory
    cmd_all[step_i].append("cd ../")

    #delete frame/subswath interferograms and amplitudes at this point.
    if (insar_arg.deleteFlag).lower() == 'true':
        cmd1 = "$INSAR_ZERODOP_SCR/delete_file.py -f f*/*.int"
        cmd2 = "$INSAR_ZERODOP_SCR/delete_file.py -f f*/*.amp"
        cmd3 = "$INSAR_ZERODOP_SCR/delete_file.py -f f*/*/*.int"
        cmd4 = "$INSAR_ZERODOP_SCR/delete_file.py -f f*/*/*.amp"
        cmd_all[step_i].extend([cmd1, cmd2, cmd3, cmd4])


#############################################################################
#   STEP. coreg_rdr2topo
#############################################################################
    step_i += 1

    cmd = "$INSAR_ZERODOP_SCR/topo.py -m {} -d {} -a {} -o {} -z {} -l {} -s {} -i {} -k {} -rlks {} -alks {}".format(
        fn.masterFrame, 
        fn.dem, 
        fn.latitude, 
        fn.longitude, 
        fn.height, 
        fn.los, 
        fn.sim, 
        fn.inc, 
        fn.msk,
        insar_arg.nrlks0,
        insar_arg.nalks0)

    cmd_all[step_i].extend(["cd {}".format(insar_dir), cmd, "cd ../"])


#############################################################################
#   STEP. coreg_topo2rdr
#############################################################################
    step_i += 1

    cmd = "$INSAR_ZERODOP_SCR/geo2rdr.py -s {} -a {} -o {} -z {} -r {} -i {} -rlks {} -alks {}".format(
        fn.slaveFrame, 
        fn.latitude, 
        fn.longitude, 
        fn.height, 
        fn.rangeOffset, 
        fn.azimuthOffset,
        insar_arg.nrlks0,
        insar_arg.nalks0)

    cmd_all[step_i].extend(["cd {}".format(insar_dir), cmd, "cd ../"])


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

    cmd_all[step_i].extend(["cd {}".format(insar_dir), cmd, "cd ../"])


#############################################################################
#   STEP. rect
#############################################################################
    step_i += 1

    cmd = "$INSAR_ZERODOP_SCR/rect.py -i {} -o {} -t {} -f {} -r {} -a {} -r0 {} -a0 {}".format(
        fn.rangeOffset, 
        fn.rectRangeOffset, 
        fn.affineTransformation, 
        'real', 
        insar_arg.nrlks0 * insar_arg.nrlksMatch,
        insar_arg.nalks0 * insar_arg.nalksMatch,
        insar_arg.nrlks0,
        insar_arg.nalks0)

    cmd_all[step_i].extend(["cd {}".format(insar_dir), cmd, "cd ../"])


#############################################################################
#   STEP. diff
#############################################################################
    step_i += 1

    cmd = "$INSAR_ZERODOP_SCR/flat.py -m {} -i {} -r {} -d {} -rlks {}".format(
        fn.masterFrame, 
        fn.interferogram, 
        fn.rectRangeOffset, 
        fn.differentialInterferogram,
        insar_arg.nrlks0)
    cmd_all[step_i].extend(["cd {}".format(insar_dir), cmd, "cd ../"])


#############################################################################
#   STEP. look
#############################################################################
    step_i += 1

    if insar_arg.nrlks != 1 or insar_arg.nalks != 1:
        cmd1 = "$INSAR_ZERODOP_SCR/look.py -i {} -o {} -r {} -a {}".format(fn.differentialInterferogram, fn.multilookDifferentialInterferogram, insar_arg.nrlks, insar_arg.nalks) 
        cmd2 = "$INSAR_ZERODOP_SCR/look.py -i {} -o {} -r {} -a {}".format(fn.amplitude, fn.multilookAmplitude, insar_arg.nrlks, insar_arg.nalks)
        cmd3 = "$INSAR_ZERODOP_SCR/look.py -i {} -o {} -r {} -a {}".format(fn.msk, fn.multilookMsk, insar_arg.nrlks, insar_arg.nalks)
        cmd4 = "$INSAR_ZERODOP_SCR/look.py -i {} -o {} -r {} -a {}".format(fn.latitude, fn.multilookLatitude, insar_arg.nrlks, insar_arg.nalks)
        cmd5 = "$INSAR_ZERODOP_SCR/look.py -i {} -o {} -r {} -a {}".format(fn.longitude, fn.multilookLongitude, insar_arg.nrlks, insar_arg.nalks)
        cmd6 = "$INSAR_ZERODOP_SCR/look.py -i {} -o {} -r {} -a {}".format(fn.height, fn.multilookHeight, insar_arg.nrlks, insar_arg.nalks)
        #using isce's look program, tested correctness
        #tested file types: 1. BIP, 2 bands
        #                   2. BIL, 2 bands (the following case)
        cmd7 = "looks.py -i {} -o {} -r {} -a {}".format(fn.los, fn.multilookLos, insar_arg.nrlks, insar_arg.nalks)
        cmd_all[step_i].extend(["cd {}".format(insar_dir), cmd1, cmd2, cmd3, cmd4, cmd5, cmd6, cmd7, "cd ../"])
    else:
        cmd1 = "cp {} {}".format(fn.differentialInterferogram, fn.multilookDifferentialInterferogram)
        cmd2 = "cp {} {}".format(fn.differentialInterferogram + '.xml', fn.multilookDifferentialInterferogram + '.xml')
        cmd3 = "cp {} {}".format(fn.differentialInterferogram + '.vrt', fn.multilookDifferentialInterferogram + '.vrt')

        cmd4 = "cp {} {}".format(fn.amplitude, fn.multilookAmplitude)
        cmd5 = "cp {} {}".format(fn.amplitude + '.xml', fn.multilookAmplitude + '.xml')
        cmd6 = "cp {} {}".format(fn.amplitude + '.vrt', fn.multilookAmplitude + '.vrt')

        cmd7 = "cp {} {}".format(fn.msk, fn.multilookMsk)
        cmd8 = "cp {} {}".format(fn.msk + '.xml', fn.multilookMsk + '.xml')
        cmd9 = "cp {} {}".format(fn.msk + '.vrt', fn.multilookMsk + '.vrt')

        cmd10 = "cp {} {}".format(fn.latitude, fn.multilookLatitude)
        cmd11 = "cp {} {}".format(fn.latitude + '.xml', fn.multilookLatitude + '.xml')
        cmd12 = "cp {} {}".format(fn.latitude + '.vrt', fn.multilookLatitude + '.vrt')

        cmd13 = "cp {} {}".format(fn.longitude, fn.multilookLongitude)
        cmd14 = "cp {} {}".format(fn.longitude + '.xml', fn.multilookLongitude + '.xml')
        cmd15 = "cp {} {}".format(fn.longitude + '.vrt', fn.multilookLongitude + '.vrt')

        cmd16 = "cp {} {}".format(fn.height, fn.multilookHeight)
        cmd17 = "cp {} {}".format(fn.height + '.xml', fn.multilookHeight + '.xml')
        cmd18 = "cp {} {}".format(fn.height + '.vrt', fn.multilookHeight + '.vrt')

        cmd19 = "cp {} {}".format(fn.los, fn.multilookLos)
        cmd20 = "cp {} {}".format(fn.los + '.xml', fn.multilookLos + '.xml')
        cmd21 = "cp {} {}".format(fn.los + '.vrt', fn.multilookLos + '.vrt')

        cmd_all[step_i].extend(["cd {}".format(insar_dir), cmd1, cmd2, cmd3, cmd4, cmd5, cmd6, cmd7, cmd8, cmd9, cmd10, cmd11, cmd12, cmd13, cmd14, cmd15, cmd16, cmd17, cmd18, cmd19, cmd20, cmd21, "cd ../"])


#############################################################################
#   STEP. coherence
#############################################################################
    step_i += 1

    #for interferometry with different modes, we need to equalize master and slave amplitudes
    #for stripmap mode, the amplitude and interferogram amplitude are too big at this stage, so scale them
    #This program makes the slave magnitude equal to master magnitude, and scale the magnitude and interferogram using user-specified value
    if insar_arg.combination == 0:
        ratio = 100000.0
    else:
        ratio = 1.0
    cmd1 = "$INSAR_ZERODOP_SCR/scale_int_amp.py -inf {} -amp {} -ratio {}".format(
        fn.multilookDifferentialInterferogram, 
        fn.multilookAmplitude,
        ratio)

    #method phase_gradient does not work, there must be bug in it.
    cmd2 = "$INSAR_ZERODOP_SCR/coherence.py -i {} -a {} -c {}".format(
        fn.multilookDifferentialInterferogram, 
        fn.multilookAmplitude, 
        fn.multilookCoherence)
    cmd_all[step_i].extend(["cd {}".format(insar_dir), cmd1, cmd2, "cd ../"])


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
    cmd_all[step_i].extend(["cd {}".format(insar_dir), cmd, "cd ../"])


#############################################################################
#   STEP. unwrap
#############################################################################
    step_i += 1

    cmd1 = "$INSAR_ZERODOP_SCR/snaphuMCF.py -mframe {} -inf {} -cor {} -unw {} -munw {} -rlks {} -alks {}".format(
        fn.masterFrame, 
        fn.filteredInterferogram, 
        fn.multilookPhsig, 
        fn.unwrappedInterferogram, 
        fn.unwrappedMaskedInterferogram, 
        insar_arg.nrlks0 * insar_arg.nrlks, 
        insar_arg.nalks0 * insar_arg.nalks)
    cmd2 = "$INSAR_ZERODOP_SCR/replace_mag.py -unw {} -inf {} -msk 0".format(
        fn.unwrappedMaskedInterferogram, 
        fn.multilookDifferentialInterferogram)
    cmd3 = "$INSAR_ZERODOP_SCR/replace_mag.py -unw {} -inf {} -msk 0".format(
        fn.unwrappedInterferogram, 
        fn.multilookDifferentialInterferogram)

    cmd_all[step_i].extend(["cd {}".format(insar_dir), cmd1, cmd2, cmd3, "cd ../"])


#############################################################################
#   STEP. geocode
#############################################################################
    step_i += 1

    # geocode considering offsets between radar amplitude and simulated amplitude
    #cmd1 = "$INSAR_ZERODOP_SCR/geo_off.py -m {} -d {} -i {} -o {} -b {} -r {} -a {} -t {} -x {} -y {}".format(fn.masterFrame, fn.demGeo, fn.unwrappedInterferogram, fn.geoUnwrappedInterferogram, insar_arg.bbox, insar_arg.nrlks0*insar_arg.nrlks, insar_arg.nalks0*insar_arg.nalks, fn.affineTransformation, insar_arg.nrlks0*insar_arg.nrlksMatch, insar_arg.nalks0*insar_arg.nalksMatch)
    #cmd2 = "$INSAR_ZERODOP_SCR/geo_off.py -m {} -d {} -i {} -o {} -b {} -r {} -a {} -t {} -x {} -y {}".format(fn.masterFrame, fn.demGeo, fn.unwrappedMaskedInterferogram, fn.geoUnwrappedMaskedInterferogram, insar_arg.bbox, insar_arg.nrlks0*insar_arg.nrlks, insar_arg.nalks0*insar_arg.nalks, fn.affineTransformation, insar_arg.nrlks0*insar_arg.nrlksMatch, insar_arg.nalks0*insar_arg.nalksMatch)
    #cmd3 = "$INSAR_ZERODOP_SCR/geo_off.py -m {} -d {} -i {} -o {} -b {} -r {} -a {} -t {} -x {} -y {}".format(fn.masterFrame, fn.demGeo, fn.multilookCoherence, fn.geoMultilookCoherence, insar_arg.bbox, insar_arg.nrlks0*insar_arg.nrlks, insar_arg.nalks0*insar_arg.nalks, fn.affineTransformation, insar_arg.nrlks0*insar_arg.nrlksMatch, insar_arg.nalks0*insar_arg.nalksMatch)
    #cmd4 = "$INSAR_ZERODOP_SCR/geo_off.py -m {} -d {} -i {} -o {} -b {} -r {} -a {} -t {} -x {} -y {}".format(fn.masterFrame, fn.demGeo, fn.multilookLos, fn.geoMultilookLos, insar_arg.bbox, insar_arg.nrlks0*insar_arg.nrlks, insar_arg.nalks0*insar_arg.nalks, fn.affineTransformation, insar_arg.nrlks0*insar_arg.nrlksMatch, insar_arg.nalks0*insar_arg.nalksMatch)

    cmd1 = "$INSAR_ZERODOP_SCR/geo.py -m {} -d {} -i {} -o {} -b={} -r {} -a {}".format(fn.masterFrame, fn.demGeo, fn.unwrappedInterferogram, fn.geoUnwrappedInterferogram, insar_arg.bbox, insar_arg.nrlks0*insar_arg.nrlks, insar_arg.nalks0*insar_arg.nalks)
    cmd2 = "$INSAR_ZERODOP_SCR/geo.py -m {} -d {} -i {} -o {} -b={} -r {} -a {}".format(fn.masterFrame, fn.demGeo, fn.unwrappedMaskedInterferogram, fn.geoUnwrappedMaskedInterferogram, insar_arg.bbox, insar_arg.nrlks0*insar_arg.nrlks, insar_arg.nalks0*insar_arg.nalks)
    cmd3 = "$INSAR_ZERODOP_SCR/geo.py -m {} -d {} -i {} -o {} -b={} -r {} -a {}".format(fn.masterFrame, fn.demGeo, fn.multilookCoherence, fn.geoMultilookCoherence, insar_arg.bbox, insar_arg.nrlks0*insar_arg.nrlks, insar_arg.nalks0*insar_arg.nalks)
    cmd4 = "$INSAR_ZERODOP_SCR/geo.py -m {} -d {} -i {} -o {} -b={} -r {} -a {}".format(fn.masterFrame, fn.demGeo, fn.multilookLos, fn.geoMultilookLos, insar_arg.bbox, insar_arg.nrlks0*insar_arg.nrlks, insar_arg.nalks0*insar_arg.nalks)

    cmd_all[step_i].extend(["cd {}".format(insar_dir), cmd1, cmd2, cmd3, cmd4, "cd ../"])


    return cmd_all


if __name__ == '__main__':
    '''
    Main driver.
    '''

    inps = cmdLineParse()

    insar_arg = read_insar_arg(inps.input)
    fn = set_filename(insar_arg)
    cmd_all = get_cmd(insar_arg, fn)


    if inps.start == 'read_data':
        prep_inputfile(insar_arg, fn)

    run_record_cmd(cmd_all, inps.start, inps.end, os.path.splitext(os.path.basename(__file__))[0] + '.log')








