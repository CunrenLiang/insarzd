InSAR processor for ALOS-2 multi-mode SAR data and ionospheric correction
Cunren LIANG
Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA

##########################################
# What's special about this processor
##########################################
This package implements a number of algorithms discribed in the 'references'
section of this file.

ONE APP TO PROCESS THEM ALL
* Support all ALOS-2 acquisition modes (except spotlight, to be added shortly)
* Combinations include ScanSAR-ScanSAR, ScanSAR-stripmap and stripmap-stripmap
* Support making interferograms across all the three wavelengths of ALOS-2
* Automatic mosaicking of multiple subswaths and frames

KEY FEATURES
* Both full-aperture and burst-by-burst InSAR processing workflows implemented
* One app to process them all
* ScanSAR system parameters estimated from azimuth spectrum
* Automatic estimation of the start times of raw burst from azimuth spectrum
* High precision burst synchronization calculation
* MBF filter for removing non-overlap spectra caused by burst misalignment and Doppler centroid frequency difference
* High precision estimation of offsets between subswaths and frames
* High precision and automatic mosaic of subswaths and frames
* Automatic ionospheric correction

* Burst extraction from full-aperture images
* Burst-by-burst interferometric processing
* Azimuth offset from burst spectral diversity
* Azimuth FM rate error estimation 
* Ionospheric correction for burst-by-burst InSAR processing

The burst-by-burst InSAR processing will be added later.


##########################################
# installation
##########################################
1. ISCE
ISCE can be downloaded at: https://github.com/isce-framework/isce2.

In the update folder of insarzd.zip, you will find a number of updated files. They are located in the same
directories in ISCE. Replace the files in ISCE with the corresponding ones in 'update' in insarzd.zip. You
can do this by running copyfile.py under 'update' directory. Follow instructions of this command.

To install isce, follow the isce installing instructions.


2. insarzd
After unpacking insarzd.zip,

cd insarzd/src
./install.sh

set environment variables:

export PYTHONPATH=YOURPATH/insarzd/scripts/pac:$PYTHONPATH
export INSAR_ZERODOP_SCR=YOURPATH/insarzd/scripts
export INSAR_ZERODOP_BIN=YOURPATH/insarzd/bin
export PATH=$INSAR_ZERODOP_SCR:$INSAR_ZERODOP_BIN:$PATH

replace YOURPATH with your own path.


##########################################
#processing
##########################################
1. DEM
This processor does not download DEM. You need to download DEM by yourself. Here is how you
can download a DEM.

#only need to change this before running
SNWE="-4 4 -81 -73"

#1 arcsec for geometric coregistration
dem.py -a stitch -b ${SNWE} -k -s 1 -c -f
rm *.zip *.log *.dem *.dem.vrt *.dem.xml

#3 arcsec for geocoding
mkdir 3
cd 3
dem.py -a stitch -b ${SNWE} -k -s 3 -c -f
rm *.zip *.log *.dem *.dem.vrt *.dem.xml
cd ../

In *.dem.wgs84.xml, there are following lines

    <property name="file_name">
        <value>demLat_S01_N01_Lon_W081_W073.dem.wgs84</value>
        <doc>Name of the image file.</doc>
    </property>

replace DEM file name with full path + DEM file name; otherwise ISCE will not recognize the file.


2. ALOS-2 data and parameter file
For each date, if you have multiple frames, unpack them in the same directory. Example parameter files
are provided in examples.zip. Follow an example parameter file in examples/alos2 to setup your own
paramter file.

alos2app_scansar-stripmap.xml  ScanSAR-stripmap interferometry
alos2app_scansar.xml           ScanSAR-ScanSAR interferometry
alos2app_stripmap.xml          stripmap-stripmap interferometry

Choose the one you need to set the parameter file. For example, if you want to do ScanSAR-ScanSAR
interferometry, choose alos2app_scansar.xml.

Normally you only have to specify "master", "slave" and "dem". For ScanSAR-ScanSAR interferometry you 
may also want to specify "starting subswath" and "ending subswath", if you don't want to process all 
the subswaths.

When setting the parameter file, you may need the frame number. Here is an example of how you can find
it. Below is a JAXA SLC product,

0000168233_001001_ALOS2183010690-171012.zip

After you unpack the JAXA SLC product, you will find an image file like:

IMG-HH-ALOS2183010685-171012-FBDR1.1__A
                 ----
The number 0685 (with underline) is the frame number. DON'T use the frame number in the zip file name,
as it may be incorrect (like the above example). 

ATTENTION FOR SCANSAR DATA ONLY: When you order the data, you should order full-aperture product.
Currently this is the product that is usable for interferometry, and is supported by the processor.


3. processing
If you don't do ionospheric correction, run the following command to process ALOS-2 data:

alos2app.py -i alos2app_scansar.xml

This will do the whole insar processing. If you do ionospheric correction, run the following commands:

alos2app.py -i alos2app_scansar.xml -e coherence
ion.py -i alos2app_scansar.xml
alos2app.py -i alos2app_scansar.xml -s filter

In both cases change the parameter filename accordingly before running, if your parameter filename
is different. In the second case, the first command does the insar processing until the coherence
calculation step. The second command calculates ionospheric phase and apply ionospheric correction. 
The third command does the remaining insar processing steps.


4. results
After the processing you may find files and directories like these:

alos2app.log              commands called by alos2app.py
alos2app_scansar.xml      your paramter file
f*_*                      frame processing directory
insar                     insar processing directory after frame processing
ion                       ionosphere calculation directory
ion.log                   commands called by ion.py

The processor first processes each frame until the generation of the interferogram in the directory
f*_*. The remaining insar processing is done in the directory insar. The final processing result is
insar/filt_diff_*-*_5rlks_2alks_msk.unw.geo

For ionospheric correction, check the following interferograms
insar/diff_*-*_5rlks_2alks_ori.int    original interferogram
insar/diff_*-*_5rlks_2alks.int        ionosphere-corrected interferogram

Check the following files for SAR parameters:
insar/*.slc.par.xml

Check the following file for baseline:
insar/*_baseline.xml

For ScanSAR-ScanSAR interferometry, check the following file for burst synchronization:
f*/s*/burst_synchronization.txt


5. sign conventions of final products
Note here master is the earlier acquistion. If your master is the later acquistion, the sign will be
opposite.

For regular InSAR product, moving down is +. This is theoretically and experimentally verified.
(e.g. hawaii/alos2/a089/180508-180522/filt_diff_180508-180522_8rlks_16alks_msk.unw.geo),

For ScanSAR MAI, moving south is +. This is experimentally verified.
(e.g. 1. hawaii/alos2/d185/sd/azd_180120-180512_28rlks_8alks_msk.unw.geo, 2. iran_2017/d71/171004-171115_burst)


##########################################
#references
##########################################
This package implements a number of algorithms, which are described in the following papers.

1. InSAR processing:
C. Liang and E. J. Fielding, "Interferometry with ALOS-2 full-aperture ScanSAR data,"
IEEE Transactions on Geoscience and Remote Sensing, vol. 55, no. 5, pp. 2739-2750,
May 2017.

2. Ionospheric correction:
C. Liang and E. J. Fielding, "Measuring azimuth deformation with L-band ALOS-2
ScanSAR interferometry," IEEE Transactions on Geoscience and Remote Sensing, vol.
55, no. 5, pp. 2725-2738, May 2017.

C. Liang, Z. Liu, E. J. Fielding, and R. Bürgmann, “InSAR time series analysis 
of L-band wide-swath SAR data acquired by ALOS-2,” IEEE Transactions on Geoscience 
and Remote Sensing, vol. 56, no. 8, pp. 4492-4506, Aug. 2018.



##########################################
#Frequently asked questions
##########################################
1. What data does the processor support?
Basically, ALOS-2 can acquire data in spotlight, stripmap and ScanSAR modes. In each mode,
the data may have different sample size, coverage, look side (right/left), range band, 
azimuth band, wavelength, and so on. In addition, the two images used for interferometry
can be acquired in different acquistion modes such as ScanSAR-stripmap interferometry. 
As long as the master and slave images meet the following requirements:

acquired on the same track
have enough spatial coverage
have enough range overlap band
have enough azimuth overlap band

an interferogram can be created by the processor. For the processing of the various data types
, the only difference on the user side is to choose the appropariate parameter file. For example,
if you want to do ScanSAR-ScanSAR interferometry, choose alos2app_scansar.xml.

Currently, the processor does not support the processing of spotlight data. Programs have been
developed for processing COSMO-SkyMed spotlight data, but are not yet tested for ALOS-2 spotlight 
data. It may be supported in the next version. For interferometry with data acquired at different 
wavelengths, the resulting interferogram might have a residual range ramp. This is probably caused 
by the relative wavelength errors of the two wavelengths.


2. Issues with ionospheric correction?
According to our experience, ionospheric correction works for ~70% of the interferograms. Because
it relies on coherence and phase unwrapping, it does not work in some cases. These include:

data have low coherence
the majority of the imaged area is low coherence area like lake, ocean...
the imaged area is completely divided into several areas by low coherence area

In addition to the above issues, there are also data-mode-related issues.
(1) ScanSAR-ScanSAR interferometry. While you can process one single subswath, it's better to process
more than one subswath. This is good for ionospheric correction.

(2) ScanSAR-stripmap interferometry and interferometry with data of different range bands. Because
of the low effective number of looks and the possible small overlap of the two range bands, ionospheric
correciton is likely not working well.

(3) Range distortions in JAXA product. This mostly happens in stripmap-stripmap interferometry using data
not covering Japan. If you see very dense fringes in the corrected inteferogram, probably it is caused by
this problem. This has been reported to JAXA and JAXA is working on debugging the focusing program.

UPDATE: On November 20, 2018 (JST), JAXA updated the software for PALSAR-2 standard products. Therefore, if
your product is ordered after this time, you don't have this problem.


3. How do I improve ionospheric correction?
Sometimes you may find that the ionospheric phase automatically calculated using default parameters
are not good enough. In this case, you may want to adjust the parameters by yourself. To do this, open
ion.log, at the end of the file you will find some commands like this:

$INSAR_ZERODOP_SCR/ion/ion_filt.py -lframe master_lower.slc.pck -uframe master_upper.slc.pck -lower lower_40rlks_16alks.unw -upper upper_40rlks_16alks_adj.unw -cor lower_40rlks_16alks.cor -ion filt_5rlks_2alks.ion -nrli 40 -nali 224 -nrlo 5 -nalo 28 -ion_fit 1 -order 4 -gsize_max 151 -gsize_min 31 -gsize_smt 31
imageMath.py -e='a*exp(-1.0*J*b)' --a=../../insar/diff_150219-150402_5rlks_2alks_ori.int --b=filt_5rlks_2alks.ion -s BIP -t cfloat -o ../../insar/diff_150219-150402_5rlks_2alks.int

The first command calculates ionsopheric phase, and the second one do the correction. In the first 
command, you can first try changing the argument of -order with values such as 1~5. Then you can also 
try changing the argument of -gsize_max.


4. ISCE forum
Participate in discussion with the user/developer community at the following website:
http://earthdef.caltech.edu/projects/isce_forum/boards


















