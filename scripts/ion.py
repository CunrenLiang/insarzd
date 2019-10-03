#!/usr/bin/env python3

#Cunren Liang
#JPL/Caltech

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
            epilog="estimate ionospheric phase\nPROCESSING STEPS: sub_band, mosaic_subswath, mosaic_frame, diff, look, coherence, unwrap, cal_ion\n\
                    ")
    parser.add_argument('-i', '--input', dest='input', type=str, required=True,
            help = 'input xml configure file')
    parser.add_argument('-s','--start', dest='start', type=str, default='sub_band',
            help = 'starting step')
    parser.add_argument('-e','--end', dest='end', type=str, default='cal_ion',
            help = 'ending step')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


def prep_dir(insar_arg):
#prepare directory

    ion_dir = 'ion'
    sub_band_dir = ['lower', 'upper']

    os.mkdir(ion_dir)
    os.chdir(ion_dir)
    for sub_band_dir_i in sub_band_dir:
        os.mkdir(sub_band_dir_i)
        os.chdir(sub_band_dir_i)
        for i, (frame_m, frame_s) in enumerate(zip(insar_arg.masterFrames, insar_arg.slaveFrames)):
            frame_dir = 'f{}'.format(frame_m)
            os.mkdir(frame_dir)
            os.chdir(frame_dir)
            if insar_arg.combination == 1 or insar_arg.combination == 2:
                for j in range(insar_arg.startingSubswath, insar_arg.endingSubswath+1):
                    subswath_dir = 's{}'.format(j)
                    os.mkdir(subswath_dir)
            os.chdir('../')
        os.chdir('../')
    os.chdir('../')


def get_cmd(insar_arg, fn):

    #subswath mosaicking directory
    frame_mosaic_dir = 'mosaic'
    #final insar processing directory
    insar_dir = 'insar'
    #ionospheric correction directory
    ion_dir = 'ion'
    #subband directory
    sub_band_dir = ['lower', 'upper']
    #final ionospheric phase calculation directory
    ion_cal_dir = 'ion_cal'


    #set your steps here
    cmd_all = [['sub_band'], ['mosaic_subswath'], ['mosaic_frame'], ['diff'], ['look'], ['coherence'], ['unwrap'], ['cal_ion']]


#############################################################################
#   STEP. sub_band
#############################################################################
    step_i = 0


    cmd_all[step_i].append("cd {}".format(ion_dir))
    for sub_band_dir_i in sub_band_dir:
        cmd_all[step_i].append("cd {}".format(sub_band_dir_i))

        #get subswath offset for both master and slave
        if sub_band_dir_i == 'lower':
            if (insar_arg.combination == 1 or insar_arg.combination == 2) and (insar_arg.endingSubswath - insar_arg.startingSubswath + 1 >= 2):
                #loop over frames
                for frame in insar_arg.masterFrames:
                    cmd_all[step_i].append("cd f{}".format(frame))
                    #1. offset of master has already been generated in InSAR processing
                    smatch_dir = 'smatch_{}'.format(fn.master)
                    cmd_all[step_i].append("mkdir {}".format(smatch_dir))
                    cmd_all[step_i].append("cp ../../../f{}/{}/subswath_offset.txt ./{}".format(frame, frame_mosaic_dir, smatch_dir))
                    #2. create a directory to do match for slave
                    #   no need to do this for ScanSAR-stripmap interferometry
                    if insar_arg.combination == 2:
                        smatch_dir = 'smatch_{}'.format(fn.slave)
                        cmd_all[step_i].append("mkdir {}".format(smatch_dir))
                        cmd_all[step_i].append("cd {}".format(smatch_dir))
                        #offset of slave
                        cmd = "$INSAR_ZERODOP_SCR/subswath_offset.py"
                        for i in range(insar_arg.startingSubswath, insar_arg.endingSubswath+1):
                            cmd += " -s ../../../../f{}/s{}/{}".format(frame, i, fn.slaveSlc)
                        for i in range(insar_arg.startingSubswath, insar_arg.endingSubswath+1):
                            cmd += " -f ../../../../f{}/s{}/{}".format(frame, i, fn.slaveFrame)
                        cmd += " -o subswath_offset.txt -rlks 1 -alks 10"
                        cmd_all[step_i].append(cmd)
                        cmd_all[step_i].append("cd ../")
                    cmd_all[step_i].append("cd ../")

        #get commands for filtering slc and forming interferogram
        #subband indicator
        if sub_band_dir_i == 'lower':
            lu = 0
        else:
            lu = 1
        #master
        cmd1_0 = "$INSAR_ZERODOP_SCR/ion/sub_band.py -oslc {oslc} -oframe {oframe} -lu {lu}".format(
            oslc = fn.masterSlc,
            oframe = fn.masterFrame,
            lu = lu)
        #slave
        cmd2_0 = "$INSAR_ZERODOP_SCR/ion/sub_band.py -oslc {oslc} -oframe {oframe} -lu {lu}".format(
            oslc = fn.slaveSlc,
            oframe = fn.slaveFrame,
            lu = lu)
        #match
        cmd3_0 = "$INSAR_ZERODOP_SCR/match_resamp_isce.py -m {master} -s {slave} -sframe {slaveframe} -int {interferogram} -amp {amplitude} -rlks {nrlks0} -alks {nalks0}".format(
            master = fn.masterSlc,
            slave = fn.slaveSlc,
            slaveframe = fn.slaveFrame, #it should be good to use the subband frame pickle file
            interferogram = fn.interferogram,
            amplitude = fn.amplitude,
            nrlks0 = insar_arg.nrlks0,
            nalks0 = insar_arg.nalks0)
        #no need to keep these slcs
        cmd4 = 'rm *.slc'

        cmd1 = cmd2 = cmd3 = ''
        #loop over frame and subswath directories
        for frame in insar_arg.masterFrames:
            cmd_all[step_i].append("cd f{}".format(frame))
            if insar_arg.combination == 0:
                ori_dir = '../../../f{}'.format(frame)
                cmd1 = cmd1_0 + ' -islc {}/{} -iframe {}/{}'.format(ori_dir, fn.masterSlc, ori_dir, fn.masterFrame)
                cmd2 = cmd2_0 + ' -islc {}/{} -iframe {}/{}'.format(ori_dir, fn.slaveSlc, ori_dir, fn.slaveFrame)
                cmd3 = cmd3_0 + ' -offset {}/cull.off'.format(ori_dir)
                cmd_all[step_i].extend([cmd1, cmd2, cmd3, cmd4])
            else:
                for subswath in range(insar_arg.startingSubswath, insar_arg.endingSubswath+1):
                    cmd_all[step_i].append("cd s{}".format(subswath))
                    ori_dir = '../../../../f{}/s{}'.format(frame, subswath)
                    cmd1 = cmd1_0 + ' -islc {}/{} -iframe {}/{}'.format(ori_dir, fn.masterSlc, ori_dir, fn.masterFrame)
                    cmd2 = cmd2_0 + ' -islc {}/{} -iframe {}/{}'.format(ori_dir, fn.slaveSlc, ori_dir, fn.slaveFrame)
                    cmd3 = cmd3_0 + ' -offset {}/cull.off'.format(ori_dir)

                    if sub_band_dir_i == 'lower':
                        offset_dir0 = ".."
                    else:
                        offset_dir0 = "../../../{}/f{}".format(sub_band_dir[0], frame)

                    #master
                    if insar_arg.endingSubswath - insar_arg.startingSubswath + 1 >= 2:
                        cmd1 += " -offset {}/smatch_{}/subswath_offset.txt -sn {}".format(offset_dir0, fn.master, subswath-insar_arg.startingSubswath)
                    #slave: no need to do this for stripmap (slave) in ScanSAR-stripmap interferometry
                    if insar_arg.combination == 2 and (insar_arg.endingSubswath - insar_arg.startingSubswath + 1 >= 2):
                        cmd2 += " -offset {}/smatch_{}/subswath_offset.txt -sn {}".format(offset_dir0, fn.slave, subswath-insar_arg.startingSubswath)
                    cmd_all[step_i].extend([cmd1, cmd2, cmd3, cmd4])
                    cmd_all[step_i].append("cd ../")
            cmd_all[step_i].append("cd ../")

        cmd_all[step_i].append("cd ../")
    cmd_all[step_i].append("cd ../")


#############################################################################
#   STEP. mosaic_subswath
#############################################################################
    step_i += 1

    cmd_all[step_i].append("cd {}".format(ion_dir))
    for sub_band_dir_i in sub_band_dir:
        cmd_all[step_i].append("cd {}".format(sub_band_dir_i))
        ############################################################
        #loop overlp frames
        ns = insar_arg.endingSubswath - insar_arg.startingSubswath + 1
        if (insar_arg.combination == 1 or insar_arg.combination == 2) and (ns > 1):
            for frame in insar_arg.masterFrames:
                cmd_all[step_i].append("cd f{}".format(frame))
                cmd_all[step_i].append("mkdir {}".format(frame_mosaic_dir))
                cmd_all[step_i].append("cd {}".format(frame_mosaic_dir))
                ori_dir = "../../../../f{}".format(frame)
                #mosaic
                cmd2 = "$INSAR_ZERODOP_SCR/subswath_mosaic.py"
                for i in range(insar_arg.startingSubswath, insar_arg.endingSubswath+1):
                    cmd2 += " -f {}/s{}/{}".format(ori_dir, i, fn.masterFrame)
                for i in range(insar_arg.startingSubswath, insar_arg.endingSubswath+1):
                    cmd2 += " -i ../s{}/{}".format(i, fn.interferogram)
                for i in range(insar_arg.startingSubswath, insar_arg.endingSubswath+1):
                    cmd2 += " -a ../s{}/{}".format(i, fn.amplitude)
                cmd2 += " -offset {}/{}/subswath_offset.txt -io {} -ao {} -fo {}".format(ori_dir, frame_mosaic_dir, fn.interferogram, fn.amplitude, fn.masterFrame)
                cmd2 += " -imethod 0 -amethod 0 -rlks {} -alks {} -phc 1 -drlks 1 -dalks 4".format(insar_arg.nrlks0, insar_arg.nalks0)
                cmd_all[step_i].append(cmd2)
                cmd_all[step_i].append("cd ../")
                cmd_all[step_i].append("cd ../")
        ############################################################
        cmd_all[step_i].append("cd ../")
    cmd_all[step_i].append("cd ../")


#############################################################################
#   STEP. mosaic_frame
#############################################################################
    step_i += 1

    cmd_all[step_i].append("cd {}".format(ion_dir))
    for sub_band_dir_i in sub_band_dir:
        cmd_all[step_i].append("cd {}".format(sub_band_dir_i))
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
            cmd_all[step_i].append("mv ../f{}/{}{} ./".format(insar_arg.masterFrames[0], dirsm, fn.interferogram))
            cmd_all[step_i].append("mv ../f{}/{}{} ./".format(insar_arg.masterFrames[0], dirsm, fn.amplitude))
            cmd_all[step_i].append("cp ../f{}/{}{} ./".format(insar_arg.masterFrames[0], dirsm, fn.interferogram+'.xml'))
            cmd_all[step_i].append("cp ../f{}/{}{} ./".format(insar_arg.masterFrames[0], dirsm, fn.amplitude+'.xml'))
            cmd_all[step_i].append("cp ../f{}/{}{} ./".format(insar_arg.masterFrames[0], dirsm, fn.interferogram+'.vrt'))
            cmd_all[step_i].append("cp ../f{}/{}{} ./".format(insar_arg.masterFrames[0], dirsm, fn.amplitude+'.vrt'))
        else:
            if insar_arg.combination == 0:
               dirsm = ''
            else:
                if ns == 1:
                    dirsm = 's{}/'.format(insar_arg.startingSubswath)
                else:
                    dirsm = '{}/'.format(frame_mosaic_dir)
            #mosaic amplitude
            cmd2 = "$INSAR_ZERODOP_SCR/frame_mosaic.py"
            for frame in insar_arg.masterFrames:
                cmd2 += " -f ../../../f{}/{}{}".format(frame, dirsm, fn.masterFrame)
            for frame in insar_arg.masterFrames:
                cmd2 += " -i ../f{}/{}{}".format(frame, dirsm, fn.amplitude)
            cmd2 += " -offset ../../../{}/frame_offset.txt -io {} -fo {} -rlks {} -alks {} -phc 0".format(insar_dir, fn.amplitude, fn.masterFrame, insar_arg.nrlks0, insar_arg.nalks0)
            #mosaic interferogram
            cmd3 = "$INSAR_ZERODOP_SCR/frame_mosaic.py"
            for frame in insar_arg.masterFrames:
                cmd3 += " -f ../../../f{}/{}{}".format(frame, dirsm, fn.masterFrame)
            for frame in insar_arg.masterFrames:
                cmd3 += " -i ../f{}/{}{}".format(frame, dirsm, fn.interferogram)
            cmd3 += " -offset ../../../{}/frame_offset.txt -io {} -fo {} -rlks {} -alks {} -phc 1".format(insar_dir, fn.interferogram, fn.masterFrame, insar_arg.nrlks0, insar_arg.nalks0)
            cmd_all[step_i].extend([cmd2, cmd3])
        #####################################################################
        cmd_all[step_i].append("cd ../")
        #delete frame/subswath interferograms and amplitudes at this point.
        if (insar_arg.deleteFlag).lower() == 'true':
            cmd1 = "$INSAR_ZERODOP_SCR/delete_file.py -f f*/*.int"
            cmd2 = "$INSAR_ZERODOP_SCR/delete_file.py -f f*/*.amp"
            cmd3 = "$INSAR_ZERODOP_SCR/delete_file.py -f f*/*/*.int"
            cmd4 = "$INSAR_ZERODOP_SCR/delete_file.py -f f*/*/*.amp"
            cmd_all[step_i].extend([cmd1, cmd2, cmd3, cmd4])
        cmd_all[step_i].append("cd ../")
    cmd_all[step_i].append("cd ../")


#############################################################################
#   STEP. diff
#############################################################################
    step_i += 1

    #define some variables here first
    lower_dir = "../{}/insar".format(sub_band_dir[0])
    upper_dir = "../{}/insar".format(sub_band_dir[1])
    rgoff = "../../insar/{}".format(fn.rectRangeOffset)
    ml = '{}rlks_{}alks'.format(insar_arg.nrlksIon, insar_arg.nalksIon)
    ml_ori = '{}rlks_{}alks'.format(insar_arg.nrlks, insar_arg.nalks)


    #processing
    cmd_all[step_i].append("mkdir {}/{}".format(ion_dir, ion_cal_dir))
    cmd_all[step_i].append("cd {}/{}".format(ion_dir, ion_cal_dir))
    cmd1 = '$INSAR_ZERODOP_SCR/ion/change_frequency.py -iframe ../../{}/{} -oframe master_lower.slc.pck -lu 0'.format(
        insar_dir,
        fn.masterFrame)
    cmd2 = '$INSAR_ZERODOP_SCR/ion/change_frequency.py -iframe ../../{}/{} -oframe master_upper.slc.pck -lu 1'.format(
        insar_dir,
        fn.masterFrame)
    cmd3 = '$INSAR_ZERODOP_SCR/flat.py -m master_lower.slc.pck -i {}/{} -r {} -d lower.int -rlks {}'.format(
        lower_dir,
        fn.interferogram,
        rgoff,
        insar_arg.nrlks0)
    cmd4 = '$INSAR_ZERODOP_SCR/flat.py -m master_upper.slc.pck -i {}/{} -r {} -d upper.int -rlks {}'.format(
        upper_dir,
        fn.interferogram,
        rgoff,
        insar_arg.nrlks0)
    cmd_all[step_i].extend([cmd1, cmd2, cmd3, cmd4])
    cmd_all[step_i].append("cd ../../")


#############################################################################
#   STEP. look
#############################################################################
    step_i += 1

    cmd_all[step_i].append("cd {}/{}".format(ion_dir, ion_cal_dir))
    if insar_arg.nrlks == 1 and insar_arg.nalks == 1:
        raise Exception('number of range looks: {} and number of azimuth looks: {} are too small for ionospheric correction'.format(insar_arg.nrlks, insar_arg.nalks))
    else:
        cmd1 = '$INSAR_ZERODOP_SCR/look.py -i lower.int -o lower_{}.int -r {} -a {}'.format(
            ml,
            insar_arg.nrlksIon,
            insar_arg.nalksIon)
        cmd2 = '$INSAR_ZERODOP_SCR/look.py -i upper.int -o upper_{}.int -r {} -a {}'.format(
            ml,
            insar_arg.nrlksIon,
            insar_arg.nalksIon)
        cmd3 = '$INSAR_ZERODOP_SCR/look.py -i {}/{} -o lower_{}.amp -r {} -a {}'.format(
            lower_dir,
            fn.amplitude,
            ml,
            insar_arg.nrlksIon,
            insar_arg.nalksIon)
        cmd4 = '$INSAR_ZERODOP_SCR/look.py -i {}/{} -o upper_{}.amp -r {} -a {}'.format(
            upper_dir,
            fn.amplitude,
            ml,
            insar_arg.nrlksIon,
            insar_arg.nalksIon)
        cmd_all[step_i].extend([cmd1, cmd2, cmd3, cmd4])
    cmd_all[step_i].append("cd ../../")


#############################################################################
#   STEP. coherence
#############################################################################
    step_i += 1

    cmd_all[step_i].append("cd {}/{}".format(ion_dir, ion_cal_dir))
    cmd1 = '$INSAR_ZERODOP_SCR/coherence.py -i lower_{ml}.int -a lower_{ml}.amp -c lower_{ml}.cor'.format(
        ml = ml)
    cmd2 = '$INSAR_ZERODOP_SCR/coherence.py -i upper_{ml}.int -a upper_{ml}.amp -c upper_{ml}.cor'.format(
        ml = ml)
    cmd_all[step_i].extend([cmd1, cmd2])
    cmd_all[step_i].append("cd ../../")


#############################################################################
#   STEP. unwrap
#############################################################################
    step_i += 1

    cmd_all[step_i].append("cd {}/{}".format(ion_dir, ion_cal_dir))

    cmd1 = "$INSAR_ZERODOP_SCR/snaphuMCF.py -mframe master_lower.slc.pck -inf lower_{ml}.int -cor lower_{ml}.cor -unw lower_{ml}.unw -munw lower_{ml}_msk.unw -rlks {rlks} -alks {alks}".format(
        ml = ml, 
        rlks = insar_arg.nrlks0 * insar_arg.nrlksIon, 
        alks = insar_arg.nalks0 * insar_arg.nalksIon)
    cmd2 = "$INSAR_ZERODOP_SCR/snaphuMCF.py -mframe master_upper.slc.pck -inf upper_{ml}.int -cor upper_{ml}.cor -unw upper_{ml}.unw -munw upper_{ml}_msk.unw -rlks {rlks} -alks {alks}".format(
        ml = ml, 
        rlks = insar_arg.nrlks0 * insar_arg.nrlksIon, 
        alks = insar_arg.nalks0 * insar_arg.nalksIon)
    cmd3 = '$INSAR_ZERODOP_SCR/ion/fix_unw_error.py -unwl lower_{ml}.unw -ccl lower_{ml}.unw.conncomp -unwu upper_{ml}.unw -ccu upper_{ml}.unw.conncomp -cor lower_{ml}.cor -unwc upper_{ml}_adj.unw -cot 0.1'.format(
        ml = ml)

    cmd_all[step_i].extend([cmd1, cmd2, cmd3])
    cmd_all[step_i].append("cd ../../")


#############################################################################
#   STEP. cal_ion
#############################################################################
    step_i += 1

    cmd_all[step_i].append("cd {}/{}".format(ion_dir, ion_cal_dir))

    #backup original first
    ori_int = os.path.splitext(fn.multilookDifferentialInterferogram)[0] + '_ori.int'
    cmd1 = "mv ../../insar/{} ../../insar/{}".format(fn.multilookDifferentialInterferogram, ori_int)
    cmd2 = "mv ../../insar/{} ../../insar/{}".format(fn.multilookDifferentialInterferogram+'.xml', ori_int+'.xml')
    cmd3 = "mv ../../insar/{} ../../insar/{}".format(fn.multilookDifferentialInterferogram+'.vrt', ori_int+'.vrt')

    #ionosphere
    cmd4 = '$INSAR_ZERODOP_SCR/ion/ion_filt.py -lframe master_lower.slc.pck -uframe master_upper.slc.pck -lower lower_{ml}.unw -upper upper_{ml}_adj.unw -cor lower_{ml}.cor -ion filt_{ml_ori}.ion -nrli {rlks} -nali {alks} -nrlo {rlks_ori} -nalo {alks_ori} -ion_fit 1 -order 4 -gsize_max {gsize_max} -gsize_min {gsize_min} -gsize_smt {gsize_smt}'.format(
        ml = ml,
        ml_ori = ml_ori,
        rlks = insar_arg.nrlks0 * insar_arg.nrlksIon,
        alks = insar_arg.nalks0 * insar_arg.nalksIon,        
        rlks_ori = insar_arg.nrlks0 * insar_arg.nrlks,
        alks_ori = insar_arg.nalks0 * insar_arg.nalks,
        gsize_max = 151,
        gsize_min = 31,
        gsize_smt = 31)

    #correct
    cmd5 = "imageMath.py -e='a*exp(-1.0*J*b)' --a={ori_int} --b=filt_{ml_ori}.ion -s BIP -t cfloat -o {new_int}".format(
        ori_int = '../../insar/'+ori_int,
        ml_ori = ml_ori,
        new_int = '../../insar/'+fn.multilookDifferentialInterferogram)

    cmd_all[step_i].extend([cmd1, cmd2, cmd3, cmd4, cmd5])
    cmd_all[step_i].append("cd ../../")


    return cmd_all


if __name__ == '__main__':
    '''
    Main driver.
    '''

    inps = cmdLineParse()

    insar_arg = read_insar_arg(inps.input)
    fn = set_filename(insar_arg)
    cmd_all = get_cmd(insar_arg, fn)


    if inps.start == 'sub_band':
        prep_dir(insar_arg)

    run_record_cmd(cmd_all, inps.start, inps.end, os.path.splitext(os.path.basename(__file__))[0] + '.log')








