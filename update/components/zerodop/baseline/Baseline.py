#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# copyright: 2010 to the present, california institute of technology.
# all rights reserved. united states government sponsorship acknowledged.
# any commercial use must be negotiated with the office of technology transfer
# at the california institute of technology.
# 
# this software may be subject to u.s. export control laws. by accepting this
# software, the user agrees to comply with all applicable u.s. export laws and
# regulations. user has the responsibility to obtain export licenses,  or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.
# 
# installation and use of this software is restricted by a license agreement
# between the licensee and the california institute of technology. it is the
# user's responsibility to abide by the terms of the license agreement.
#
# Author: Giangi Sacco
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



import math
import datetime
import logging
from iscesys.Component.Component import Component, Port
from isceobj.Orbit.Orbit import StateVector
import numpy as np 

BASELINE_LOCATION = Component.Parameter('baselineLocation',
        public_name = 'BASELINE_LOCATION',
        default = 'all',
        type=str,
        mandatory=False,
        doc = 'Location at which to compute baselines - "all" implies top, middle, bottom of master image, "top" implies near start of master image, "bottom" implies at bottom of master image, "middle" implies near middle of master image. To be used in case there is a large shift between images.')



class Baseline(Component):

    family = 'baseline'
    logging_name = 'isce.zerodop.baseline'

    parameter_list = (BASELINE_LOCATION,)

    # Calculate the baseline components between two frames
    def baseline(self):

        from isceobj.Util.geo.ellipsoid import Ellipsoid
        from isceobj.Planet.Planet import Planet
        for port in self.inputPorts:
            port()
            
        planet = Planet(pname='Earth')
        refElp = Ellipsoid(a=planet.ellipsoid.a, e2=planet.ellipsoid.e2, model='WGS84')


        if self.baselineLocation.lower() == 'all':
            print('Using entire span of image for estimating baselines')
            masterTime = [self.masterFrame.getSensingStart(),self.masterFrame.getSensingMid(),self.masterFrame.getSensingStop()]
        elif self.baselineLocation.lower() == 'middle':
            print('Estimating baselines around center of master image')
            masterTime = [self.masterFrame.getSensingMid() - datetime.timedelta(seconds=1.0), self.masterFrame.getSensingMid(), self.masterFrame.getSensingMid() + datetime.timedelta(seconds=1.0)]

        elif self.baselineLocation.lower() == 'top':
            print('Estimating baselines at top of master image')
            masterTime =  [self.masterFrame.getSensingStart(), self.masterFrame.getSensingStart() + datetime.timedelta(seconds=1.0), self.masterFrame.getSensingStart() + datetime.timedelta(seconds=2.0)]
        elif self.baselineLocation.lower() == 'bottom':
            print('Estimating baselines at bottom of master image')
            masterTime =  [self.masterFrame.getSensingStop() - datetime.timedelta(seconds=2.0), self.masterFrame.getSensingStop() - datetime.timedelta(seconds=1.0), self.masterFrame.getSensingStop()]
        else:
            raise Exception('Unknown baseline location: {0}'.format(self.baselineLocation))


        s = [0., 0., 0.]
        bpar = []
        bperp = []
        azoff = []
        rgoff = []

        for i in range(3):
            # Calculate the Baseline at the start of the scene, mid-scene, and the end of the scene
            # First, get the position and velocity at the start of the scene
            # Calculate the distance moved since the last baseline point
            s[i] = (masterTime[i] - masterTime[0]).total_seconds()
            

            masterSV = self.masterOrbit.interpolateOrbit(masterTime[i], method='hermite')
            rng = self.startingRange1
            target = self.masterOrbit.pointOnGround(masterTime[i], rng, side=self.masterFrame.getInstrument().getPlatform().pointingDirection)

            slaveTime, slvrng = self.slaveOrbit.geo2rdr(target) 
            slaveSV = self.slaveOrbit.interpolateOrbit(slaveTime, method='hermite')

            targxyz = np.array(refElp.LLH(target[0], target[1], target[2]).ecef().tolist())
            mxyz = np.array(masterSV.getPosition())
            mvel = np.array(masterSV.getVelocity())
            sxyz = np.array(slaveSV.getPosition())

            aa = np.linalg.norm(sxyz-mxyz)

            costheta = (rng*rng + aa*aa - slvrng*slvrng)/(2.*rng*aa)

#            print(aa, costheta)
            bpar.append(aa*costheta)

            perp = aa * np.sqrt(1 - costheta*costheta)
            direction = np.sign(np.dot( np.cross(targxyz-mxyz, sxyz-mxyz), mvel))
            bperp.append(direction*perp)

            ####Azimuth offset
            slvaz = (slaveTime - self.slaveFrame.sensingStart).total_seconds() * self.prf2
            masaz = s[i] * self.prf1
            azoff.append(slvaz - masaz)

            ####Range offset
            slvrg =  (slvrng - self.startingRange2)/self.rangePixelSize2
            masrg = (rng - self.startingRange1) / self.rangePixelSize1
            rgoff.append(slvrg - masrg)

       
#        print(bpar)
#        print(bperp)

        #Calculating baseline
        parBaselinePolynomialCoefficients = np.polyfit(s,bpar,2)
        perpBaselinePolynomialCoefficients = np.polyfit(s,bperp,2)
        
        # Populate class attributes 
        self.BparMean = parBaselinePolynomialCoefficients[-1]
        self.BparRate = parBaselinePolynomialCoefficients[1]
        self.BparAcc = parBaselinePolynomialCoefficients[0]
        self.BperpMean = perpBaselinePolynomialCoefficients[-1]
        self.BperpRate = perpBaselinePolynomialCoefficients[1]
        self.BperpAcc = perpBaselinePolynomialCoefficients[0]

        delta = (self.masterFrame.getSensingStart() - masterTime[0]).total_seconds()
        self.BparTop = np.polyval(parBaselinePolynomialCoefficients, delta)
        self.BperpTop = np.polyval(perpBaselinePolynomialCoefficients, delta)

        delta = (self.masterFrame.getSensingStop() - masterTime[0]).total_seconds()
        self.BparBottom = np.polyval(parBaselinePolynomialCoefficients, delta)
        self.BperpBottom = np.polyval(perpBaselinePolynomialCoefficients, delta)
       
        return azoff, rgoff
            
    def setMasterRangePixelSize(self,pixelSize):
        self.rangePixelSize1 = pixelSize
        return

    def setSlaveRangePixelSize(self,pixelSize):
        self.rangePixelSize2 = pixelSize
        return

    def setMasterStartingRange(self,range):
        self.startingRange1 = range
        return

    def setSlaveStartingRange(self,range):
        self.startingRange2 = range
        return

    def setMasterPRF(self,prf):
        self.prf1 = prf
        return

    def setSlavePRF(self,prf):
        self.prf2 = prf
        return
    
    def getHBaselineTop(self):
        return self.hBaselineTop

    def getHBaselineRate(self):
        return self.hBaselineRate

    def getHBaselineAcc(self):
        return self.hBaselineAcc

    def getVBaselineTop(self):
        return self.vBaselineTop

    def getVBaselineRate(self):
        return self.vBaselineRate

    def getVBaselineAcc(self):
        return self.vBaselineAcc

    def getPBaselineTop(self):
        return self.pBaselineTop

    def getPBaselineBottom(self):
        return self.pBaselineBottom



        
    def addMasterFrame(self):
        frame = self._inputPorts.getPort(name='masterFrame').getObject()
        self.startingRange1 = frame.getStartingRange()
        self.prf1 = frame.getInstrument().getPulseRepetitionFrequency()
        self.rangePixelSize1 = frame.getInstrument().getRangePixelSize()
        self.masterOrbit = frame.getOrbit()
        self.masterFrame = frame

    def addSlaveFrame(self):
        frame = self._inputPorts.getPort(name='slaveFrame').getObject()
        self.startingRange2 = frame.getStartingRange()
        self.slaveOrbit = frame.getOrbit()
        self.prf2 = frame.getInstrument().getPulseRepetitionFrequency()
        self.rangePixelSize2 = frame.getInstrument().getRangePixelSize()
        self.slaveFrame = frame
        
    def __init__(self, name=''):
        super(Baseline, self).__init__(family=self.__class__.family, name=name)
        self.masterOrbit = None
        self.slaveOrbit = None
        self.masterFrame = None
        self.slaveFrame = None
        self.rangePixelSize1 = None
        self.rangePixelSize2 = None
        self.startingRange1 = None
        self.startingRange2 = None
        self.prf1 = None
        self.prf2 = None
        self.lookSide = None
        self.BparMean = None
        self.BparRate = None
        self.BparAcc = None
        self.BperpMean = None
        self.BperpRate = None
        self.BperpAcc = None
        self.BperpTop = None
        self.BperpBottom = None
        self.BparTop = None
        self.BperpBottom = None
        self.logger = logging.getLogger('isce.zerodop.baseline')
        self.createPorts()
        
        # Satisfy the old Component
        self.dictionaryOfOutputVariables = {}        
        self.dictionaryOfVariables = {}        
        self.descriptionOfVariables = {}
        self.mandatoryVariables = []
        self.optionalVariables = []
        return None

    def createPorts(self):
        
        # Set input ports
        # It looks like we really need two orbits, a time, range and azimuth pixel sizes
        # the two starting ranges, a planet, and the two prfs
        # These provide the orbits
        # These provide the range and azimuth pixel sizes, starting ranges, 
        # satellite heights and times for the first lines
        masterFramePort = Port(name='masterFrame',method=self.addMasterFrame)  
        slaveFramePort = Port(name='slaveFrame',method=self.addSlaveFrame)       
        self._inputPorts.add(masterFramePort)
        self._inputPorts.add(slaveFramePort)
        return None

        
    def __str__(self):
        retstr = "Initial Baseline estimates \n"
        retlst = ()
        retstr += "Parallel Baseline Top: %s\n"
        retlst += (self.BparTop,)
        retstr += "Perpendicular Baseline Top: %s\n"
        retlst += (self.BperpTop,)
        retstr += "Parallel Baseline Bottom: %s\n"
        retlst += (self.BparBottom,)
        retstr += "Perpendicular Baseline Bottom: %s \n"
        retlst += (self.BperpBottom,)
        return retstr % retlst      
