#!/usr/bin/env python3

#program for filtering ionosphere phase
#Cunren Liang, APR. 2016
#Jet Propulsion Laboratory, California Institute of Technology

import os
import sys
import glob
import shutil
import pickle
import datetime
import argparse
import numpy as np
import numpy.matlib
import scipy.signal as ss

import isce
import isceobj
from imageMath import IML
from isceobj.Constants import SPEED_OF_LIGHT

from crlpac import runCmd
from crlpac import create_xml
from crlpac import gaussian
from crlpac import create_multi_index2
from crlpac import fit_surface
from crlpac import cal_surface

def weight_fitting(ionos, cor, width, length, nrli, nali, nrlo, nalo, order, coth):
    '''
    ionos:  input ionospheric phase
    cor:    coherence of the interferogram
    width:  file width
    length: file length
    nrli:   number of range looks of the input interferograms
    nali:   number of azimuth looks of the input interferograms
    nrlo:   number of range looks of the output ionosphere phase
    nalo:   number of azimuth looks of the ioutput ionosphere phase
    order:  the order of the polynomial for fitting ionosphere phase estimates
    coth:   coherence threshhold for ionosphere phase estimation
    '''

    lengthi = int(length/nali)
    widthi = int(width/nrli)
    lengtho = int(length/nalo)
    widtho = int(width/nrlo)

    #calculate output index
    rgindex = create_multi_index2(widtho, nrli, nrlo)
    azindex = create_multi_index2(lengtho, nali, nalo)

    #convert coherence to weight
    cor = cor**2/(1.009-cor**2)

    #look for data to use
    flag = (cor>coth)*(ionos!=0)
    point_index = np.nonzero(flag)
    m = point_index[0].shape[0]

    #calculate input index matrix
    x0=np.matlib.repmat(np.arange(widthi), lengthi, 1)
    y0=np.matlib.repmat(np.arange(lengthi).reshape(lengthi, 1), 1, widthi)

    x = x0[point_index].reshape(m, 1)
    y = y0[point_index].reshape(m, 1)
    z = ionos[point_index].reshape(m, 1)
    w = cor[point_index].reshape(m, 1)

    #convert to higher precision type before use
    x=np.asfarray(x,np.float64)
    y=np.asfarray(y,np.float64)
    z=np.asfarray(z,np.float64)
    w=np.asfarray(w,np.float64)
    coeff = fit_surface(x, y, z, w, order)

    #convert to higher precision type before use
    rgindex=np.asfarray(rgindex,np.float64)
    azindex=np.asfarray(azindex,np.float64)
    phase_fit = cal_surface(rgindex, azindex.reshape(lengtho, 1), coeff, order)

    #format: widtho, lengtho, single band float32
    return phase_fit


def adaptive_filtering_winsize(ionos, wgt, cor, size_max, size_min, sigma):
######################################
    #ionos: ionosphere
    #wgt: weight
    #cor: coherence
    #size_max: maximum window size
    #size_min: minimum window size
    #sigma: sigma
######################################

    length = (ionos.shape)[0]
    width = (ionos.shape)[1]

    #we fill the zero areas in cor and wgt, using only the zero areas of ionos as indication of zero or not.
    #the filtered ionosphere is not zero where ionos is not zero, and
    #the filtered ionosphere may be not zero where ionos is zero.    
    #fill the zero areas in cor and wgt
    cor_win_size = np.int32(np.around(length/2.0));
    scale = ss.fftconvolve((cor!=0), np.ones((cor_win_size, cor_win_size)), mode='same')
    fcor = ss.fftconvolve(cor, np.ones((cor_win_size, cor_win_size)), mode='same') / (scale + (scale==0))
    cor[np.nonzero(cor==0)] = fcor[np.nonzero(cor==0)]

    scale = ss.fftconvolve((wgt!=0), np.ones((cor_win_size, cor_win_size)), mode='same')
    fwgt = ss.fftconvolve(wgt, np.ones((cor_win_size, cor_win_size)), mode='same') / (scale + (scale==0))
    wgt[np.nonzero(wgt==0)] = fwgt[np.nonzero(wgt==0)]

    wgt *= (ionos!=0)
    cor *= (ionos!=0)

    cor_max = np.amax(cor[np.nonzero(cor!=0)])
    cor_min = np.amin(cor[np.nonzero(cor!=0)])
    print('coherence max: %6.3f, min: %6.3f'%(cor_max, cor_min))

    size_num = 100
    step = (size_max - size_min) / (size_num)

    ionos_all=np.zeros((length, width, size_num*2))
    size_num = 0
    for size in np.arange(size_min, size_max+step, step, dtype=np.float32):
        
        size = np.around(size)
        if size % 2 == 0:
            size += 1
        if size_num % 10 == 0:
            print('filtering with window size: %5d, sigma: %8.2f'%(size, sigma*size/size_max))

        g2d = gaussian(int(size), sigma*size/size_max, scale=1.0)
        scale = ss.fftconvolve(wgt, g2d, mode='same')
        ionos_all[:, :, size_num] = ss.fftconvolve(ionos*wgt, g2d, mode='same') / (scale + (scale==0))
        
        size_num += 1

    size_x = size_max - (cor-cor_min)/(cor_max - cor_min)*(size_max-size_min)
    #size_x = size_x * flag_cor + size_med * (flag_cor==0)
    size_index = np.int32(np.around((size_x - size_min)/step))
    size_index[np.nonzero(size_index>size_num-1)] = size_num-1
    size_index[np.nonzero(size_index<0)] = 0

    #get filtering result
    index = np.nonzero(np.ones((length, width))) + (size_index.reshape(1, length*width), )
    ionos_filt = ionos_all[index]

    return ionos_filt.reshape(length, width)


def cmdLineParse():
    '''
    Command line parser.
    '''
    parser = argparse.ArgumentParser( description='estimate ionosphere phase')
    parser.add_argument('-lframe', dest='lframe', type=str, required=True,
            help = 'master frame pickle file of lower band')
    parser.add_argument('-uframe', dest='uframe', type=str, required=True,
            help = 'master frame pickle file of upper band')
    parser.add_argument('-lower', dest='lower', type=str, required=True,
            help = 'unwrapped interferogram of lower band')
    parser.add_argument('-upper', dest='upper', type=str, required=True,
            help = 'unwrapped interferogram of upper band')
    parser.add_argument('-cor', dest='cor', type=str, required=True,
            help = 'coherence of the interferogram')
    parser.add_argument('-ion', dest='ion', type=str, required=True,
            help = 'ionosphere phase')
    parser.add_argument('-nrli', dest='nrli', type=int, default=5,
            help = 'number of range looks of the input interferograms')
    parser.add_argument('-nali', dest='nali', type=int, default=28,
            help = 'number of azimuth looks of the input interferograms')
    parser.add_argument('-nrlo', dest='nrlo', type=int, default=5,
            help = 'number of range looks of the output ionosphere phase')
    parser.add_argument('-nalo', dest='nalo', type=int, default=28,
            help = 'number of azimuth looks of the ioutput ionosphere phase')

    parser.add_argument('-ion_fit', dest='ion_fit', type=int, default=1,
            help = 'whether do fitting. 0: no. 1 yes (default)')
    parser.add_argument('-order', dest='order', type=int, default=4,
            help = 'order of the polynomial for fitting')
    parser.add_argument('-gsize_max', dest='gsize_max', type=int, default=151,
            help = 'maximum Gaussian filter window size, must be odd. Default: 151.')
    parser.add_argument('-gsize_min', dest='gsize_min', type=int, default=35,
            help = 'minimum Gaussian filter window size, must be odd. Default: 35.')
    parser.add_argument('-gsize_smt', dest='gsize_smt', type=int, default=-1,
            help = 'Gaussian filter window size for smoothing out artifacts due to coherence changes, must be odd. Default: -1 (no smoothing).')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':

    inps = cmdLineParse()

    size_max = inps.gsize_max
    size_min = inps.gsize_min
    sigma = inps.gsize_max / 2.0

    with open(inps.lframe, 'rb') as f:
        lframe = pickle.load(f)
    with open(inps.uframe, 'rb') as f:
        uframe = pickle.load(f)

    length = lframe.getNumberOfLines()
    width  = lframe.getNumberOfSamples()

    lengthi = int(length/inps.nali)
    widthi = int(width/inps.nrli)
    lengtho = int(length/inps.nalo)
    widtho = int(width/inps.nrlo)

    fl = SPEED_OF_LIGHT / lframe.radarWavelegth
    fu = SPEED_OF_LIGHT / uframe.radarWavelegth
    f0 = (fl + fu) / 2.0

    #get files
    lower = (np.fromfile(inps.lower, dtype=np.float32).reshape(lengthi*2, widthi))[1:lengthi*2:2, :]
    upper = (np.fromfile(inps.upper, dtype=np.float32).reshape(lengthi*2, widthi))[1:lengthi*2:2, :]
    cor = (np.fromfile(inps.cor, dtype=np.float32).reshape(lengthi*2, widthi))[1:lengthi*2:2, :]

    flag = (lower!=0)*(upper!=0)

    #ionosphere
    ionos = fl * fu * (lower * fu - upper * fl) / f0 / (fu**2 - fl**2)
    ionos *= flag

    #output original ionosphere phase of input size
    input_size_ionf_ori = 'ori.ion'
    ionos.astype(np.float32).tofile(input_size_ionf_ori)
    create_xml(input_size_ionf_ori, widthi, lengthi, 'float')

    #fit
    if inps.ion_fit == 1:
        ion_fit = weight_fitting(ionos, cor, width, length, inps.nrli, inps.nali, inps.nrli, inps.nali, inps.order, 0.05)
    else:
        ion_fit = np.ones((lengthi, widthi)) * np.mean(ionos[np.nonzero(flag)], dtype=np.float64)
    ion_fit *= flag

    #remove
    ionos -= ion_fit

    #filter
    #weight
    wgt = cor**2
    ionos_filt_0 = adaptive_filtering_winsize(ionos, wgt, cor, size_max, size_min, sigma)
    #smooth out the artifacts caused by coherence change
    if inps.gsize_smt != -1:
        g2d = gaussian(inps.gsize_smt, inps.gsize_smt/2.0, scale=1.0)
        scale = ss.fftconvolve((ionos_filt_0!=0), g2d, mode='same')
        ionos_filt_0 = ss.fftconvolve(ionos_filt_0, g2d, mode='same') / (scale + (scale==0))
    ionos_filt_0 *= flag

    #add back
    ionos_filt_0 += ion_fit

    #output ionosphere phase of input size
    input_size_ionf = 'input_size.ion'
    ionos_filt_0.astype(np.float32).tofile(input_size_ionf)
    create_xml(input_size_ionf, widthi, lengthi, 'float')

    #resample to required output size
    rec_in = '''
Input Image File Name  (-) = {inputfile}        ! dimension of file to be rectified
Output Image File Name (-) = {outputfile}       ! dimension of output
Input Dimensions       (-) = {width1} {length1} ! across, down
Output Dimensions      (-) = {width2} {length2} ! across, down
Affine Matrix Row 1    (-) = {m11} {m12}      ! a b
Affine Matrix Row 2    (-) = {m21} {m22}      ! c d
Affine Offset Vector   (-) = {t1} {t2}        ! e f
File Type              (-) = REAL        ! [RMG, COMPLEX]
Interpolation Method   (-) = Bilinear        ! [NN, Bilinear, Sinc]

'''.format(
    inputfile = input_size_ionf,
    outputfile = inps.ion,
    width1 = widthi, length1 = lengthi,
    width2 = widtho, length2 = lengtho,
    #m11 = inps.nrli/inps.nrlo, m12 = 0.0,
    #m21 = 0.0,                 m22 = inps.nali/inps.nalo,
    #t1 = (inps.nrli-inps.nrlo)/(2.0*inps.nrlo), t2 = (inps.nali-inps.nalo)/(2.0*inps.nalo)
    m11 = inps.nrlo/inps.nrli, m12 = 0.0,
    m21 = 0.0,                 m22 = inps.nalo/inps.nali,
    t1 = (inps.nrlo-inps.nrli)/(2.0*inps.nrli), t2 = (inps.nalo-inps.nali)/(2.0*inps.nali)
    )

    rectInputFile = 'rect.in'
    rectOutputFile = 'rect.out'
    with open(rectInputFile, 'w') as f:
        f.write(rec_in)
    cmd = '$INSAR_ZERODOP_BIN/rect {} >{}'.format(rectInputFile, rectOutputFile)
    runCmd(cmd)
    create_xml(inps.ion, widtho, lengtho, 'float')
    os.remove(rectInputFile)
    os.remove(rectOutputFile)
    os.remove(input_size_ionf)
    os.remove(input_size_ionf+'.vrt')
    os.remove(input_size_ionf+'.xml')