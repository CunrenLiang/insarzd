#!/usr/bin/env python3

#Cunren Liang, Nov. 17, 2015


import os
import sys
import pickle
import shutil
import argparse
import datetime
import numpy as np

import isce
import isceobj


def create_lfm(ns, it, offset, k):
    # create linear FM signal

    # ns: number of samples
    # it: time interval of the samples
    # offset: offset
    # k: linear FM rate
    
    #offset-centered, this applies to both odd and even cases
    ht = (ns - 1) / 2.0
    t = np.arange(-ht, ht+1.0, 1)
    t = (t + offset) * it
    cj = np.complex64(1j)
    lfm = np.exp(cj * np.pi * k * t**2)

    return lfm


def next_pow2(a):
    x=2
    while x < a:
        x *= 2
    return x


def is_power2(num):

    'states if a number is a power of two'

    return num != 0 and ((num & (num - 1)) == 0)



def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(description="calulate raw burst start times of SLC focused by full-aperture")
    parser.add_argument('-s', dest='slc', type=str, required=True,
            help='SLC focused by full-aperture algorithm')
    parser.add_argument('-f', dest='frame', type=str, required=True,
            help='frame object of the SLC saved in pickle file')
    parser.add_argument('-l', dest='sl', type=int, default=500,
            help = 'start line (included, index start with 0. default: 500) used to estimate burst times')
    parser.add_argument('-c', dest='sc', type=int, default=500,
            help = 'start column (included, index start with 0. default: 500) used to estimate burst times')
    parser.add_argument('-p', dest='pow2', type=int, default=1,
            help = 'must be 1(default) or power of 2. azimuth fft length = THIS ARGUMENT * next of power of 2 of full-aperture length.')

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

    with open(inps.frame, 'rb') as fid:
        frame = pickle.load(fid)


#######################################################
#set these parameters
    width = frame.getNumberOfSamples()
    length = frame.getNumberOfLines()
    prf = frame.PRF
    nb = frame.nbraw
    nc = frame.ncraw
    fmrateCoeff = frame.fmrateCoeff
    sensing_start = frame.getSensingStart()
    p2 = inps.pow2
#######################################################    

#######################################################
#check parameters
    if not (p2 == 1 or is_power2(p2)):
        raise Exception('option -p must be 1 or power of 2\n')
#######################################################

    #calculate fmrate
    ka = np.zeros(width, dtype=np.float32)
    for i in range(width):
        ka[i] = fmrateCoeff[2] * i**2 + fmrateCoeff[1] * i**1 + fmrateCoeff[0]     
    #use the convention that ka > 0
    ka = -ka

    #area to be used for estimation
    saz = inps.sl #start line to be used (included)
    naz = int(np.round(nc)) #number of lines to be used.
    eaz = saz+naz-1 #ending line to be used (included)
    caz = int(np.round((saz+eaz)/2.0)) #central line of the lines used.
    caz_deramp = (saz+eaz)/2.0 #center of deramp signal (may be fractional line number)

    srg = inps.sc #start column to be used (included)
    nrg = 400 #number of columns to be used
    erg = srg+nrg-1 #ending column to be used (included)
    crg = int(np.round((srg+erg)/2.0)) #central column of the columns used.


#######################################################
#check parameters
    if not (saz >=0 and saz <= length - 1):
        raise Exception('wrong starting line\n')
    if not (eaz >=0 and eaz <= length - 1):
        raise Exception('wrong ending line\n')
    if not (srg >=0 and srg <= width - 1):
        raise Exception('wrong starting column\n')
    if not (erg >=0 and erg <= width - 1):
        raise Exception('wrong ending column\n')
#######################################################


    #number of lines of a full-aperture
    nfa = int(np.round(prf / ka[crg] / (1.0 / prf)))
    #use nfa to determine fft size. fft size can be larger than this
    nazfft = next_pow2(nfa) * p2

    #deramp signal  
    deramp = np.zeros((naz, nrg), dtype=np.complex64)
    for i in range(nrg):
        deramp[:, i] = create_lfm(naz, 1.0/prf, 0, -ka[i+srg])

    #read data
    with open(inps.slc, 'rb') as f:
        f.seek(saz * width * np.dtype(np.complex64).itemsize, 0)
        data = np.fromfile(f, dtype=np.complex64, count=naz*width).reshape(naz,width)
        data = data[:, srg:erg+1]

    #deramp data
    datadr = deramp * data

    #spectrum
    spec = np.fft.fft(datadr, n=nazfft, axis=0)

    #shift zero-frequency component to center of spectrum
    spec = np.fft.fftshift(spec, axes=0)

    specm=np.mean(np.absolute(spec), axis=1)

    #number of samples of the burst in frequncy domain
    nbs  = int(np.round(nb*(1.0/prf)*ka[crg]/prf*nazfft));
    #number of samples of the burst cycle in frequncy domain
    ncs  = int(np.round(nc*(1.0/prf)*ka[crg]/prf*nazfft));
    rect = np.ones(nbs, dtype=np.float32)

    #make sure orders of specm and rect are correct, so that peaks
    #happen in the same order as their corresponding bursts
    corr=np.correlate(specm, rect,'same')

    #find burst spectrum center
    ncs_rh = int(np.round((nazfft - ncs) / 2.0))
    #corr_bc = corr[ncs_rh: ncs_rh+ncs]
    #offset between spectrum center and center
    offset_spec = np.argmax(corr[ncs_rh: ncs_rh+ncs])+ncs_rh - (nazfft - 1.0) / 2.0
    #offset in number of azimuth lines
    offset_naz  = offset_spec / nazfft * prf / ka[crg] / (1.0/prf)
    

    #start line of burst (fractional line number)
    saz_burst = -offset_naz + caz_deramp - (nb - 1.0) / 2.0

    #find out the start line of all bursts (fractional line number,
    #line index start with 0, line 0 is the first SLC line)
    burstStartLines=[]
    burstStartTimes=[]
    for i in range(-100000, 100000):
        saz_burstx = saz_burst + nc * i
        st_burstx = sensing_start + datetime.timedelta(seconds=saz_burstx * (1.0/prf))
        if saz_burstx >= 0.0 and saz_burstx <= length:
            burstStartLines.append(saz_burstx)
            burstStartTimes.append(st_burstx)


    frame.burstStartLines = burstStartLines
    frame.burstStartTimes = burstStartTimes
    frame.burstStartLineEstimated = saz_burst

    #add the burst starting times to frame pickle
    os.remove(inps.frame)
    with open(inps.frame, 'wb') as f:
        pickle.dump(frame, f)

    #check result
    print('burst no     start line                 start time')
    print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    for i in range(len(burstStartLines)):
        print('{:4d} {:18.2f}          {}'.format(i, burstStartLines[i], burstStartTimes[i]))


    #output spectrum and correlation
    specm_corr = ''
    for i in range(nazfft):
        specm_corr += '{:6d}     {:f}     {:6d}     {:f}\n'.format(i, specm[i], i, corr[i])

    year = str(frame.getSensingStart().year)[2:]
    month = '%02d' % frame.getSensingStart().month
    day = '%02d' % frame.getSensingStart().day
    specm_corr_name = 'spec_corr_' + year + month + day + '.txt'
    with open(specm_corr_name, 'w') as f:
        f.write(specm_corr)


########################################################################################################
#example
#./burst_time.py -s 140809.slc -f 140809_withnbnc.slc.pck -l 100 -c 100



