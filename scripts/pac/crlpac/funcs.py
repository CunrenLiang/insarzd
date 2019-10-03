#!/usr/bin/env python3

class Dummy(object):
    pass


def getWidth(xmlfile):
    from xml.etree.ElementTree import ElementTree
    xmlfp = None
    try:
        xmlfp = open(xmlfile,'r')
        print('reading file width from: {0}'.format(xmlfile))
        xmlx = ElementTree(file=xmlfp).getroot()
        #width = int(xmlx.find("component[@name='coordinate1']/property[@name='size']/value").text)
        tmp = xmlx.find("component[@name='coordinate1']/property[@name='size']/value")
        if tmp == None:
            tmp = xmlx.find("component[@name='Coordinate1']/property[@name='size']/value")
        width = int(tmp.text)
        print("file width: {0}".format(width))
    except (IOError, OSError) as strerr:
        print("IOError: %s" % strerr)
        return []
    finally:
        if xmlfp is not None:
            xmlfp.close()
    return width


def getLength(xmlfile):
    from xml.etree.ElementTree import ElementTree
    xmlfp = None
    try:
        xmlfp = open(xmlfile,'r')
        print('reading file length from: {0}'.format(xmlfile))
        xmlx = ElementTree(file=xmlfp).getroot()
        #length = int(xmlx.find("component[@name='coordinate2']/property[@name='size']/value").text)
        tmp = xmlx.find("component[@name='coordinate2']/property[@name='size']/value")
        if tmp == None:
            tmp = xmlx.find("component[@name='Coordinate2']/property[@name='size']/value")
        length = int(tmp.text)
        print("file length: {0}".format(length))
    except (IOError, OSError) as strerr:
        print("IOError: %s" % strerr)
        return []
    finally:
        if xmlfp is not None:
            xmlfp.close()
    return length


def getInfo(xmlfile, propertyName):
    from xml.etree.ElementTree import ElementTree
    xmlfp = None
    content = None
    try:
        xmlfp = open(xmlfile,'r')
        xmlx = ElementTree(file=xmlfp).getroot()
        #search each possible propertyName
        propertyNameList = [propertyName, propertyName.lower(), propertyName.upper(), propertyName.capitalize()]
        for propertyNameX in propertyNameList:
            route = "property[@name='{}']/value".format(propertyNameX)
            contentAll = xmlx.find(route)
            if contentAll != None:
                content = contentAll.text
                content = content.strip("'") #remove the leading and trailing quote
                content = content.strip('"') #remove the leading and trailing quote
                #print("{}: {}".format(propertyName, content))
                break
    except (IOError, OSError) as strerr:
        print("IOError: %s" % strerr)
        return None
    finally:
        if xmlfp is not None:
            xmlfp.close()
    return content


def get_content(xmlfile, properties, names):
    from xml.etree.ElementTree import ElementTree

    xmlfp = None
    content = None

    propertylist = properties.split('.')
    namelist = names.split('.')

    #get route
    if len(propertylist) != len(namelist):
        raise Exception('property list length not equal to name list length!')
    else:
        route = ''
        for i in range(len(propertylist)):
            route += propertylist[i]
            if namelist[i] != '':
                route += "[@name='{}']".format(namelist[i])
            if i != len(propertylist) - 1:
                route += "/"

    #find content
    try:
        xmlfp = open(xmlfile,'r')
        root = ElementTree(file=xmlfp).getroot()
        content0 = root.find(route)
        if content0 != None:
            content = content0.text
            content = content.strip("'") #remove the possible leading and trailing quote
            content = content.strip('"') #remove the possible leading and trailing quote

    except (IOError, OSError) as strerr:
        print("IOError: %s" % strerr)
        return None
    finally:
        if xmlfp is not None:
            xmlfp.close()
    return content


def changeXmlName(xmlFile, xmlFileNew):
    import os
    import isce
    import isceobj
    from imageMath import IML    

    #currently, only sure about supporting SLC xml file
    xmlFileNew=xmlFileNew.rstrip('.xml')
    img, dataName, metaName = IML.loadImage(xmlFile)
    img.filename = xmlFileNew
    #img.setAccessMode('READ')
    os.remove(xmlFile)
    img.renderHdr()


def writeSARconfig_ALOS2(filename, leader, slc, output):
    a = '''<component name="sar">
    <property name="OUTPUT">
        <value>{output}</value>
    </property>
    <property name="LEADERFILE">
      <value>'{leader}'</value>
    </property>
    <property name="IMAGEFILE">
      <value>'{slc}'</value>
    </property>
</component>
    '''.format(
        output = output,
        leader = leader,
        slc = slc
        )

    with open(filename, 'w') as f:
        f.write(a)


def create_xml(fileName, width, length, fileType):
    import isce
    import isceobj

    if fileType == 'slc':
        image = isceobj.createSlcImage()
    elif fileType == 'int':
        image = isceobj.createIntImage()
    elif fileType == 'amp':
        image = isceobj.createAmpImage()
    elif fileType == 'rmg' or fileType == 'unw':
        image = isceobj.Image.createUnwImage()
    elif fileType == 'float':
        image = isceobj.createImage()
        image.setDataType('FLOAT')
    else:
        raise Exception('format not supported yet!\n')

    image.setFilename(fileName)
    image.setWidth(width)
    image.setLength(length)
        
    #image.setAccessMode('read')
    #image.createImage()
    image.renderHdr()
    #image.finalizeImage()


def renderParXml(frame, name):
    import math
    import numpy as np

    import isce
    import isceobj
    from isceobj.Constants import SPEED_OF_LIGHT

    catalog = isceobj.Catalog.createCatalog(name)
    par = 'frame'

    catalog.addItem('facility', frame.getProcessingFacility(), par + '.slcProcessingSoftware')
    catalog.addItem('system', frame.getProcessingSystem(), par + '.slcProcessingSoftware')
    catalog.addItem('version', frame.getProcessingSoftwareVersion(), par + '.slcProcessingSoftware')

    catalog.addItem('mission', frame.instrument.platform.getMission(), par)
    catalog.addItem('passDirection', frame.getPassDirection(), par)
    catalog.addItem('antennaLength', frame.instrument.platform.getAntennaLength(), par)
    if frame.getInstrument().getPlatform().pointingDirection == -1:
        pointingDirection = 'right'
    else:
        pointingDirection = 'left'
    catalog.addItem('antennaPointingDirection', pointingDirection, par)
    catalog.addItem('radarWavelegth', frame.instrument.radarWavelength, par)
    catalog.addItem('polarization', frame.getPolarization(), par)
    catalog.addItem('sensingStart', frame.getSensingStart(), par)    
    catalog.addItem('sensingStop', frame.getSensingStop(), par)
    catalog.addItem('startingRange', frame.getStartingRange(), par)
    catalog.addItem('numberOfLines', frame.getNumberOfLines(), par)
    catalog.addItem('numberOfSamples', frame.getNumberOfSamples(), par)
    catalog.addItem('rangeSamplingRate', frame.instrument.rangeSamplingRate, par)
    catalog.addItem('rangeBandWidth', math.fabs(frame.instrument.pulseLength * frame.instrument.chirpSlope), par)
    catalog.addItem('PRF', frame.instrument.PRF, par)

    width=frame.getNumberOfSamples()
    if hasattr(frame, 'dopCoeff'):
        catalog.addItem('DopplerCoefficients', frame.dopCoeff, par)
        doppler = []
        for i in [0, (width-1.0)/2.0, width-1.0]:
            doppler += [frame.PRF*(frame.dopCoeff[3] * math.pow(i, 3) + frame.dopCoeff[2] * math.pow(i, 2) + frame.dopCoeff[1] * math.pow(i, 1) + frame.dopCoeff[0])]
        catalog.addItem('DopplerNearMidFar', doppler, par)

    if hasattr(frame, 'fmrateCoeff'):
        catalog.addItem('azimuthFmRateCoefficients', frame.fmrateCoeff, par)
        fmrate = []
        for i in [0, (width-1.0)/2.0, width-1.0]:
            fmrate += [frame.fmrateCoeff[2] * i**2 + frame.fmrateCoeff[1] * i**1 + frame.fmrateCoeff[0]]
        catalog.addItem('azimuthFmRateNearMidFar', fmrate, par)

######################################################################################################
    #calculate pixel size
    orbit = frame.getOrbit()
    numberOfOrbit = len(orbit)
    orbitMid = orbit[int(numberOfOrbit/2)]
    vn = np.sqrt(orbitMid.velocity[0]**2 + orbitMid.velocity[1]**2 + orbitMid.velocity[2]**2)
    h = np.sqrt(orbitMid.position[0]**2 + orbitMid.position[1]**2 + orbitMid.position[2]**2)
    #earth radius in meters
    r = 6371 * 1000.0;
    #slant range pixel size
    slantRangePixelSize = 0.5 * SPEED_OF_LIGHT / frame.rangeSamplingRate
    #azimuth pixel size
    azi = 1.0 / frame.PRF
    #azimuth pixel size on track
    azimuthPixelSizeOnTrack = vn * azi
    #azimuth pixel size on ground
    azimuthPixelSizeOnGround = vn * azi * r / h

    catalog.addItem('slantRangePixelSize', slantRangePixelSize, par)
    catalog.addItem('azimuthPixelSizeOnTrack', azimuthPixelSizeOnTrack, par)
    catalog.addItem('azimuthPixelSizeOnGround', azimuthPixelSizeOnGround, par)
######################################################################################################
    orbitElementsDesignator = {'0':'preliminary',
                               '1':'decision',
                               '2':'high precision'}
    catalog.addItem('orbitQuality', orbit.orbitQuality, par)
    for i in range(numberOfOrbit):
        catalog.addItem('time', orbit[i].getTime(), par+'.orbit_{}'.format(i+1))
        catalog.addItem('position', orbit[i].getPosition(), par+'.orbit_{}'.format(i+1))
        catalog.addItem('velocity', orbit[i].getVelocity(), par+'.orbit_{}'.format(i+1))

    catalog.renderXml()


def runCmd(cmd, silent=0):
    import os

    if silent == 0:
        print("{}".format(cmd))
    status = os.system(cmd)
    if status != 0:
        raise Exception('error when running:\n{}\n'.format(cmd))


def run_record_cmd(cmd_all, start_step, end_step, cmdfile):
    import os
    import datetime

    #find starting and ending index    
    step = start_step
    start_step_index = [i for i, step_cmd in enumerate(cmd_all) if step in step_cmd]
    step = end_step
    end_step_index = [i for i, step_cmd in enumerate(cmd_all) if step in step_cmd]

    #check index
    if len(start_step_index) != 0:
        start_step_index = start_step_index[0]
    else:
        raise Exception('wrong start step')
    if len(end_step_index) != 0:
        end_step_index = end_step_index[0]
    else:
        raise Exception('wrong end step')
    if start_step_index > end_step_index:
        raise Exception('start step > end step')

    #record header
    with open(cmdfile, 'a') as f:
        header = '###################################################################\n'
        header = header + '#  processing commands started at: {}\n'.format(datetime.datetime.now())
        header = header + '###################################################################\n'
        f.write(header)

    #get absolute directory in order to always write to this file when change directory
    cmdfile = os.path.abspath(cmdfile)

    #run and record commands
    for i in range(start_step_index, end_step_index+1):
        with open(cmdfile, 'a') as f:
            f.write('\n##STEP: {}\n'.format(cmd_all[i][0]))
        for j in range(1, len(cmd_all[i])):
            with open(cmdfile, 'a') as f:
                f.write(cmd_all[i][j] + '\n')
            #if change directory, it only changes bash process's directory, it wont change python process's directory. so use python's os.chdir instead. NOT WORK WITH DIRECTORY NAME WITH SPACE YET!
            if cmd_all[i][j].split()[0]=='cd': 
                os.chdir(cmd_all[i][j].split()[1])
            else:
                runCmd(cmd_all[i][j])

    #add some blank lines
    with open(cmdfile, 'a') as f:
        f.write('\n\n\n')

#os.system cannot capture status of the primary commands. so it still has problems.
def run_record_cmd2(cmd_all, start_step, end_step, cmdfile, outputfile):
    import os
    import datetime

    #find starting and ending index    
    step = start_step
    start_step_index = [i for i, step_cmd in enumerate(cmd_all) if step in step_cmd]
    step = end_step
    end_step_index = [i for i, step_cmd in enumerate(cmd_all) if step in step_cmd]

    #check index
    if len(start_step_index) != 0:
        start_step_index = start_step_index[0]
    else:
        raise Exception('wrong start step')
    if len(end_step_index) != 0:
        end_step_index = end_step_index[0]
    else:
        raise Exception('wrong end step')
    if start_step_index > end_step_index:
        raise Exception('start step > end step')

    time = datetime.datetime.now()
    header  = '###################################################################\n'
    header += '#  processing commands started at: {}\n'.format(time)
    header += '###################################################################\n'

    #record header for command and output file
    with open(cmdfile, 'a') as f:
        f.write(header)
    with open(outputfile, 'a') as f:
        f.write(header)
    #get absolute directory in order to always write to this file when change directory
    cmdfile = os.path.abspath(cmdfile)
    outputfile = os.path.abspath(outputfile)

    #run and record commands
    for i in range(start_step_index, end_step_index+1):
        #record step
        header_step = '\n##STEP: {}\n'.format(cmd_all[i][0])
        with open(cmdfile, 'a') as f:
            f.write(header_step)
        with open(outputfile, 'a') as f:
            f.write(header_step)
        #run commands
        for j in range(1, len(cmd_all[i])):
            #record commands
            with open(cmdfile, 'a') as f:
                f.write(cmd_all[i][j] + '\n')
            print("{}".format(cmd_all[i][j]))

            #run commands and record output
            #if change directory, it only changes bash process's directory, it wont change python process's directory. so use python's os.chdir instead. NOT WORK WITH DIRECTORY NAME WITH SPACE YET!
            if cmd_all[i][j].split()[0]=='cd':
                os.chdir(cmd_all[i][j].split()[1])
            else:
                #os.system cannot capture status of the primary commands. so it still has problems.
                status = os.system('{} 2>&1 | tee -a {}'.format(cmd_all[i][j], outputfile))
                #status = os.system('{} >> {} 2>&1'.format(cmd_all[i][j], outputfile))
                if status != 0:
                    raise Exception('error when running:\n{}\n'.format(cmd_all[i][j]))

    #add some blank lines
    tail = '\n\n\n'
    with open(cmdfile, 'a') as f:
        f.write(tail)
    with open(outputfile, 'a') as f:
        f.write(tail)


def writeOffset(offset, fileName):

    offsetsPlain = ''
    for offsetx in offset:
        offsetsPlainx = "{}".format(offsetx)
        offsetsPlainx = offsetsPlainx.split()
        offsetsPlain = offsetsPlain + "{:8d} {:10.3f} {:8d} {:12.3f} {:11.5f} {:11.6f} {:11.6f} {:11.6f}\n".format(
            int(offsetsPlainx[0]),
            float(offsetsPlainx[1]),
            int(offsetsPlainx[2]),
            float(offsetsPlainx[3]),
            float(offsetsPlainx[4]),
            float(offsetsPlainx[5]),
            float(offsetsPlainx[6]),
            float(offsetsPlainx[7])
            )

    offsetFile = fileName
    with open(offsetFile, 'w') as f:
        f.write(offsetsPlain)


def meanOffset(offsets):

    rangeOffset = 0.0
    azimuthOffset = 0.0
    i = 0
    for offsetx in offsets:
        i += 1
        rangeOffset += offsetx.dx
        azimuthOffset += offsetx.dy

    rangeOffset /= i
    azimuthOffset /= i


    return [rangeOffset, azimuthOffset]


def cullOffsetbyMean(offsets, azimuthOffsetMean, threshold):
    import isce
    import isceobj
    from isceobj.Location.Offset import OffsetField,Offset

    culledOffsetField = OffsetField()
    i = 0
    for offset in offsets:
        if abs(offset.dy - azimuthOffsetMean) > threshold:
            i += 1
        else:
            culledOffsetField.addOffset(offset)

    print("{} offsets culled, with azimuth mean offset: {} and threshold: {}".format(i, azimuthOffsetMean, threshold))
    return culledOffsetField


def cullOffset(offsets, distances, numCullOffsetsLimits):
    import os
    import isce
    import isceobj
    from iscesys.StdOEL.StdOELPy import create_writer
    #offsets: offsets from ampcor
    #distances: tuple
    #numCullOffsetsLimits: tuple
    
    refinedOffsets = offsets
    for i, (distance, numCullOffsetsLimit) in enumerate(zip(distances, numCullOffsetsLimits)):

        cullOff = isceobj.createOffoutliers()
        cullOff.wireInputPort(name='offsets', object=refinedOffsets)
        cullOff.setSNRThreshold(2.0)
        cullOff.setDistance(distance)
        
        #set the tag used in the outfile. each message is precided by this tag
        #is the writer is not of "file" type the call has no effect
        logfile = "offoutliers.log"
        stdWriter = create_writer("log", "", True, filename=logfile)
        stdWriter.setFileTag("offoutliers", "log")
        stdWriter.setFileTag("offoutliers", "err")
        stdWriter.setFileTag("offoutliers", "out")
        cullOff.setStdWriter(stdWriter)

        try:
            cullOff.offoutliers()
            refinedOffsets = cullOff.getRefinedOffsetField()
            numLeft = len(refinedOffsets._offsets)
            print('Number of offsets left after %2dth culling: %5d'%(i, numLeft))
            if numLeft < numCullOffsetsLimit:
                print('*******************************************************')
                print('WARNING: Too few points left after culling: {} left'.format(numLeft))
                print('*******************************************************')
                return None
        except:
            print('*******************************************************')
            print('WARNING: unsuccessful offset culling')
            print('*******************************************************')
            return None

        os.remove(logfile)

    return refinedOffsets


def getOffset(offsets, offsetFile, cullOffsetFile, dumpFile):

    offsetsPlain = ''
    for offsetx in offsets:
        offsetsPlainx = "{}".format(offsetx)
        offsetsPlainx = offsetsPlainx.split()
        offsetsPlain = offsetsPlain + "{:8d} {:10.3f} {:8d} {:12.3f} {:11.5f} {:11.6f} {:11.6f} {:11.6f}\n".format(
            int(offsetsPlainx[0]),
            float(offsetsPlainx[1]),
            int(offsetsPlainx[2]),
            float(offsetsPlainx[3]),
            float(offsetsPlainx[4]),
            float(offsetsPlainx[5]),
            float(offsetsPlainx[6]),
            float(offsetsPlainx[7])
            )

    #offsetFile = 'offset_{}_{}.off'.format(i-1, i)
    with open(offsetFile, 'w') as f:
        f.write(offsetsPlain)

    breakFlag = 0
    for maxrms in [0.08,  0.16,  0.24]:
        #maxrms = maxrms * 0.01
        for nsig in [1.5,  1.4,  1.3,  1.2,  1.1,  1.0,  0.9]:
            #nsig = nsig * 0.1

            #cullOffsetFile = 'cull_{}_{}.off'.format(i-1, i)
            #dumpFile = 'fitoff_{}_{}.out'.format(i-1, i)
            #run fitoff here
            cmd = '$INSAR_ZERODOP_BIN/fitoff {} {} {} {} 50 > {}'.format(offsetFile, cullOffsetFile, nsig, maxrms, dumpFile)
            runCmd(cmd)

            #check number of matching points left
            with open(cullOffsetFile, 'r') as ff:
                numCullOffsets = sum(1 for linex in ff)
            if numCullOffsets < 50:
                print('offsets culling with nsig {} maxrms {}: {} left after culling, two few points'.format(nsig, maxrms, numCullOffsets))
            else:
                print('offsets culling with nsig {} maxrms {}: {} left after culling, success'.format(nsig, maxrms, numCullOffsets))
                breakFlag = 1
                break
        
        if breakFlag == 1:
            break

    if numCullOffsets < 50:
            print('*******************************************************')
            print('WARNING: Too few points left after culling: {} left'.format(numCullOffsets))
            print('*******************************************************')
            return None


    with open(dumpFile) as f:
        lines = f.readlines()

    i = 0
    for linex in lines:
        if 'Affine Matrix ' in linex:
            m11 = float(lines[i + 2].split()[0])
            m12 = float(lines[i + 2].split()[1])
            m21 = float(lines[i + 3].split()[0])
            m22 = float(lines[i + 3].split()[1])
            t1  = float(lines[i + 7].split()[0])
            t2  = float(lines[i + 7].split()[1])
            break
        i += 1    

    return [m11, m12, m21, m22, t1, t2]


# def cal_coherence(inf, win=5):
#     '''
#     Compute coherence using scipy convolve 2D.
#     '''

#     filt = np.ones((win,win))/ (1.0*win*win)

#     cJ = np.complex64(1.0j)
#     angle = np.exp(cJ * np.angle(inf))

#     cor = ss.convolve2d(angle, filt, mode='same')
#     cor[0:win-1,:] = 0.0
#     cor[-win+1:,:] = 0.0
#     cor[:,0:win-1] = 0.0
#     cor[:,-win+1:] = 0.0

#     cor = np.absolute(cor)
#     print(np.max(cor), np.min(cor))
#     #cor.astype(np.float32).tofile(f)

#     return cor

#better way to trim edges
def cal_coherence(inf, win=5):
    '''
    Compute coherence using scipy convolve 2D.
    '''
    import numpy as np
    import scipy.signal as ss


    filt = np.ones((win,win))/ (1.0*win*win)

    #calculate flag
    flag = ss.convolve2d((np.absolute(inf)!=0), filt, mode='same')

    cJ = np.complex64(1.0j)
    angle = np.exp(cJ * np.angle(inf))

    cor = ss.convolve2d(angle, filt, mode='same')
    #cor[0:win-1,:] = 0.0
    #cor[-win+1:,:] = 0.0
    #cor[:,0:win-1] = 0.0
    #cor[:,-win+1:] = 0.0

    cor = np.absolute(cor)
    cor[np.nonzero(flag < 0.999)] = 0.0;
    print(np.max(cor), np.min(cor))
    #cor.astype(np.float32).tofile(f)

    return cor


def overlapFrequency(centerfreq1, bandwidth1, centerfreq2, bandwidth2):
    startfreq1 = centerfreq1 - bandwidth1 / 2.0
    endingfreq1 = centerfreq1 + bandwidth1 / 2.0

    startfreq2 = centerfreq2 - bandwidth2 / 2.0
    endingfreq2 = centerfreq2 + bandwidth2 / 2.0

    overlapfreq = []
    if startfreq2 <= startfreq1 <= endingfreq2:
        overlapfreq.append(startfreq1)
    if startfreq2 <= endingfreq1 <= endingfreq2:
        overlapfreq.append(endingfreq1)
    
    if startfreq1 < startfreq2 < endingfreq1:
        overlapfreq.append(startfreq2)
    if startfreq1 < endingfreq2 < endingfreq1:
        overlapfreq.append(endingfreq2)

    if len(overlapfreq) != 2:
        #no overlap bandwidth
        return None
    else:
        startfreq = min(overlapfreq)
        endingfreq = max(overlapfreq)
        return [startfreq, endingfreq] 


def gaussian(size, sigma, scale = 1.0):
    import numpy as np
    import numpy.matlib

    if size % 2 != 1:
        raise Exception('size must be odd')
    hsize = (size - 1) / 2
    x = np.arange(-hsize, hsize + 1) * scale
    f = np.exp(-x**2/(2.0*sigma**2)) / (sigma * np.sqrt(2.0*np.pi))
    f2d=np.matlib.repmat(f, size, 1) * np.matlib.repmat(f.reshape(size, 1), 1, size)

    return f2d/np.sum(f2d)


def create_multi_index(width, rgl):
    import numpy as np
    #create index after multilooking
    #assuming original index start with 0
    #applies to both range and azimuth direction

    widthm = int(width/rgl)

    #create range index: This applies to both odd and even cases, "rgl = 1" case, and "rgl = 2" case
    start_rgindex = (rgl - 1.0) / 2.0
    rgindex0 = start_rgindex + np.arange(widthm) * rgl

    return rgindex0


def create_multi_index2(width2, l1, l2):
    import numpy as np
    #for number of looks of l1 and l2
    #calculate the correponding index number of l2 in the l1 array
    #applies to both range and azimuth direction

    return ((l2 - l1) / 2.0  + np.arange(width2) * l2) / l1


def fit_surface(x, y, z, wgt, order):
    import numpy as np
    import numpy.matlib
    # x: x coordinate, a column vector
    # y: y coordinate, a column vector
    # z: z coordinate, a column vector
    # wgt: weight of the data points, a column vector


    #number of data points
    m = x.shape[0]
    l = np.ones((m,1), dtype=np.float64)

#    #create polynomial
#    if order == 1:
#        #order of estimated coefficents: 1, x, y
#        a1 = np.concatenate((l, x, y), axis=1)
#    elif order == 2:
#        #order of estimated coefficents: 1, x, y, x*y, x**2, y**2
#        a1 = np.concatenate((l, x, y, x*y, x**2, y**2), axis=1)
#    elif order == 3:
#        #order of estimated coefficents: 1, x, y, x*y, x**2, y**2, x**2*y, y**2*x, x**3, y**3
#        a1 = np.concatenate((l, x, y, x*y, x**2, y**2, x**2*y, y**2*x, x**3, y**3), axis=1)
#    else:
#        raise Exception('order not supported yet\n')

    if order < 1:
        raise Exception('order must be larger than 1.\n')

    #create polynomial
    a1 = l;
    for i in range(1, order+1):
        for j in range(i+1):
            a1 = np.concatenate((a1, x**(i-j)*y**(j)), axis=1)

    #number of variable to be estimated
    n = a1.shape[1]

    #do the least squares
    a = a1 * np.matlib.repmat(np.sqrt(wgt), 1, n)
    b = z * np.sqrt(wgt)
    c = np.linalg.lstsq(a, b)[0]
    
    #type: <class 'numpy.ndarray'>
    return c


def cal_surface(x, y, c, order):
    import numpy as np
    import numpy.matlib
    #x: x coordinate, a row vector
    #y: y coordinate, a column vector
    #c: coefficients of polynomial from fit_surface
    #order: order of polynomial

    if order < 1:
        raise Exception('order must be larger than 1.\n')

    #number of lines
    length = y.shape[0]
    #number of columns, if row vector, only one element in the shape tuple
    #width = x.shape[1]
    width = x.shape[0]

    x = np.matlib.repmat(x, length, 1)
    y = np.matlib.repmat(y, 1, width)
    z = c[0] * np.ones((length,width), dtype=np.float64)

    index = 0
    for i in range(1, order+1):
        for j in range(i+1):
            index += 1
            z += c[index] * x**(i-j)*y**(j)

    return z


def read_param_for_checking_overlap(leader_file, image_file):
    import os
    import isce
    from isceobj.Sensor import xmlPrefix
    import isceobj.Sensor.CEOS as CEOS


    #read from leader file
    fsampConst = { 104: 1.047915957140240E+08,
                   52: 5.239579785701190E+07,
                   34: 3.493053190467460E+07,
                   17: 1.746526595233730E+07 }

    fp = open(leader_file,'rb')
    leaderFDR = CEOS.CEOSDB(xml=os.path.join(xmlPrefix,'alos2_slc/leader_file.xml'),dataFile=fp)
    leaderFDR.parse()
    fp.seek(leaderFDR.getEndOfRecordPosition())
    sceneHeaderRecord = CEOS.CEOSDB(xml=os.path.join(xmlPrefix,'alos2_slc/scene_record.xml'),dataFile=fp)
    sceneHeaderRecord.parse()
    fp.seek(sceneHeaderRecord.getEndOfRecordPosition())

    fsamplookup = int(sceneHeaderRecord.metadata['Range sampling rate in MHz'])
    rangeSamplingRate = fsampConst[fsamplookup]
    fp.close()
    #print('{}'.format(rangeSamplingRate))


    #read from image file
    fp = open(image_file, 'rb')
    imageFDR = CEOS.CEOSDB(xml=os.path.join(xmlPrefix,'alos2_slc/image_file.xml'), dataFile=fp)
    imageFDR.parse()
    fp.seek(imageFDR.getEndOfRecordPosition())
    imageData = CEOS.CEOSDB(xml=os.path.join(xmlPrefix,'alos2_slc/image_record.xml'), dataFile=fp)
    imageData.parseFast()

    width = imageFDR.metadata['Number of pixels per line per SAR channel']
    near_range = imageData.metadata['Slant range to 1st data sample']
    fp.close()
    #print('{}'.format(width))
    #print('{}'.format(near_range))


    return [rangeSamplingRate, width, near_range]


def check_overlap(ldr_m, img_m, ldr_s, img_s):
    import isce
    from isceobj.Constants import SPEED_OF_LIGHT
    #0                    1       2
    #[rangeSamplingRate, width, near_range]

    mparam = read_param_for_checking_overlap(ldr_m, img_m)
    sparam = read_param_for_checking_overlap(ldr_s, img_s)
    mcenter = mparam[2] + (mparam[1] - 1) / 2.0 * 0.5 * SPEED_OF_LIGHT / mparam[0]
    mwidth = (mparam[1] - 1) * 0.5 * SPEED_OF_LIGHT / mparam[0]
    scenter = sparam[2] + (sparam[1] - 1) / 2.0 * 0.5 * SPEED_OF_LIGHT / sparam[0]
    swidth = (sparam[1] - 1) * 0.5 * SPEED_OF_LIGHT / sparam[0]

    #check swath overlap
    overlap = overlapFrequency(mcenter, mwidth, scenter, swidth)
    if overlap == None:
        #no overlap
        overlap_ratio = 0.0
    else:
        overlap_ratio = (overlap[1] - overlap[0]) / mwidth

    #print('++++++++++++++++++++++++')
    #print('{}'.format(overlap_ratio))

    return overlap_ratio


def read_insar_arg(xmlfile):
    import os
    import glob

    #################################
    #define processing settings
    #################################
    insar_arg = Dummy()

    #mandatory
    #insar_arg.sensor = None #str
    insar_arg.masterDir = None # str
    insar_arg.masterFrames = None #str
    insar_arg.masterPolarization = 'HH'
    insar_arg.slaveDir = None # str
    insar_arg.slaveFrames =  None #str
    insar_arg.slavePolarization = 'HH'
    insar_arg.dem = None #str DEM dir and its name, this DEM is used for geometrical calculations, except geocoding
    insar_arg.demGeo = None #str DEM dir and its name, this DEM is used for geocoding

    #semi optional
    insar_arg.startingSubswath = 1  #int
    insar_arg.endingSubswath = 5  #int
    insar_arg.burstOverlapThreshhold = 85.0 #float

    #optional
    insar_arg.nrlks0 = 1
    insar_arg.nalks0 = 14
    insar_arg.nrlks = 5  #int
    insar_arg.nalks = 2  #int
    insar_arg.nrlksMatch = insar_arg.nrlks  #int
    insar_arg.nalksMatch = insar_arg.nalks  #int
    insar_arg.nrlksIon = insar_arg.nrlks * 8  #int
    insar_arg.nalksIon = insar_arg.nalks * 8  #int
    insar_arg.filterStrength = 0.5 #float
    insar_arg.bbox = 'useMasterDefaultBbox' #str
    insar_arg.deleteFlag = 'False' #str

    #internal
    insar_arg.sensor = 'ALOS2'
    insar_arg.combination = 0


    #################################
    #read processing settings
    #################################
    #mandatory
    insar_arg.masterDir = get_content(xmlfile, 'component.property.value', 'master.data directory.')
    if insar_arg.masterDir == None:
        raise Exception('master data directory not set!')
    else:
        insar_arg.masterDir = os.path.abspath(insar_arg.masterDir)

    insar_arg.slaveDir = get_content(xmlfile, 'component.property.value', 'slave.data directory.')
    if insar_arg.slaveDir == None:
        raise Exception('slave data directory not set!')
    else:
        insar_arg.slaveDir = os.path.abspath(insar_arg.slaveDir)

    optionalParameter = get_content(xmlfile, 'component.property.value', 'master.polarization.')
    if optionalParameter != None:
        insar_arg.masterPolarization = optionalParameter

    optionalParameter = get_content(xmlfile, 'component.property.value', 'slave.polarization.')
    if optionalParameter != None:
        insar_arg.slavePolarization = optionalParameter

##############################################################################################################
    leader_files_m = sorted(glob.glob(insar_arg.masterDir + "/" + "LED-ALOS2*-*-*"))
    leader_files_s = sorted(glob.glob(insar_arg.slaveDir + "/" + "LED-ALOS2*-*-*"))
    first_frame_m = leader_files_m[0].split('-')[-3][-4:]
    first_frame_s = leader_files_s[0].split('-')[-3][-4:]
    first_images_m = sorted(glob.glob(insar_arg.masterDir + "/" + "IMG-{}-ALOS2*{}-*-*".format(insar_arg.masterPolarization, first_frame_m)))
    first_images_s = sorted(glob.glob(insar_arg.slaveDir + "/" + "IMG-{}-ALOS2*{}-*-*".format(insar_arg.slavePolarization, first_frame_s)))

    #internal, determine mode combination here
    mode_m = leader_files_m[0].split('-')[-1][0:3]
    mode_s = leader_files_s[0].split('-')[-1][0:3]
    stripmap_modes = ['UBS', 'UBD', 'HBS', 'HBD', 'HBQ', 'FBS', 'FBD', 'FBQ']
    scansar_modes = ['WBS', 'WBD', 'WWS', 'WWD', 'VBS', 'VBD']
    if mode_m in stripmap_modes and mode_s in stripmap_modes:
        insar_arg.combination = 0
    elif mode_m in scansar_modes and mode_s in stripmap_modes:
        insar_arg.combination = 1
    elif mode_m in scansar_modes and mode_s in scansar_modes:
        insar_arg.combination = 2
    elif mode_m in stripmap_modes and mode_s in scansar_modes:
        raise Exception('For ScanSAR-stripmap inteferometry, please use ScanSAR as master')
    else:
        raise Exception('unknown acquisition modes')
##############################################################################################################

    optionalParameter = get_content(xmlfile, 'component.property.value', 'master.frames.')
    if optionalParameter != None:
        insar_arg.masterFrames = optionalParameter.strip().strip('[').strip(']').split(',')
        for i, framex in enumerate(insar_arg.masterFrames):
            insar_arg.masterFrames[i] = framex.strip()
    else:
        if len(leader_files_m) > 1:
            raise Exception('There are more than one frames in {}, but you did not specify the frames'.format(insar_arg.masterDir))
        elif len(leader_files_m) == 0:
            raise Exception('No data file found')
        else:
            insar_arg.masterFrames = [leader_files_m[0].split('-')[-3][-4:]]

    optionalParameter = get_content(xmlfile, 'component.property.value', 'slave.frames.')
    if optionalParameter != None:
        insar_arg.slaveFrames = optionalParameter.strip().strip('[').strip(']').split(',')
        for i, framex in enumerate(insar_arg.slaveFrames):
            insar_arg.slaveFrames[i] = framex.strip()
    else:
        if len(leader_files_s) > 1:
            raise Exception('There are more than one frames in {}, but you did not specify the frames'.format(insar_arg.slaveDir))
        elif len(leader_files_s) == 0:
            raise Exception('No data file found')
        else:
            insar_arg.slaveFrames = [leader_files_s[0].split('-')[-3][-4:]]

    if len(insar_arg.masterFrames) != len(insar_arg.slaveFrames):
        raise Exception('number of master frames != number of slave frames')


################################################################################################
# there may be same frames, we add sequence number to frame. 01-MAY-2017
# only three files are changed due to this: funcs.py, alos2app.py, alos2app_burst.py
# frames are used to do the following two things:
# 1. find original data files: if insar_arg.masterFrames and insar_arg.slaveFrames are changed
#                              the associated things need to be changed accordingly
# 2. handle frame directories: if insar_arg.masterFrames and insar_arg.slaveFrames are changed
#                              the associated things do not need to be changed, as frame
#                              directory can use any name.
    for i in range(len(insar_arg.masterFrames)):
        insar_arg.masterFrames[i] = "{}_".format(i+1) + insar_arg.masterFrames[i]
        insar_arg.slaveFrames[i] = "{}_".format(i+1) + insar_arg.slaveFrames[i]
################################################################################################


    insar_arg.dem = get_content(xmlfile, 'component.property.value', 'dem.geometrical coregistration.')
    if insar_arg.dem == None:
        raise Exception('dem for geometrical coregistration not set!')
    else:
        insar_arg.dem = os.path.abspath(insar_arg.dem)

    insar_arg.demGeo = get_content(xmlfile, 'component.property.value', 'dem.geocoding.')
    if insar_arg.demGeo == None:
        raise Exception('dem for geocoding!')
    else:
        insar_arg.demGeo = os.path.abspath(insar_arg.demGeo)


    #semi optional
    optionalParameter = get_content(xmlfile, 'component.property.value', 'scansar.starting subswath.')
    if optionalParameter != None:
        insar_arg.startingSubswath = int(optionalParameter)

    optionalParameter = get_content(xmlfile, 'component.property.value', 'scansar.ending subswath.')
    if optionalParameter != None:
        insar_arg.endingSubswath = int(optionalParameter)

    #for ScanSAR-stripmap interferometry, update starting subswath and ending subswath here
    if insar_arg.combination == 1:
        overlap_subswaths = []
        for i, images_x in enumerate(first_images_m):
            overlap_ratio = check_overlap(leader_files_m[0], images_x, leader_files_s[0], first_images_s[0])
            if overlap_ratio > 1.0 / 4.0:
                overlap_subswaths.append(i+1)
        if overlap_subswaths == []:
            raise Exception('There is no overlap area between the pair')
        insar_arg.startingSubswath = int(overlap_subswaths[0])
        insar_arg.endingSubswath = int(overlap_subswaths[-1])

    optionalParameter = get_content(xmlfile, 'component.property.value', 'scansar.burst overlap threshhold for mbf filtering (%).')
    if optionalParameter != None:
        insar_arg.burstOverlapThreshhold = float(optionalParameter)


    #optional
    optionalParameter = get_content(xmlfile, 'component.property.value', 'insar.number of range looks when forming interferogram.')
    if optionalParameter != None:
        insar_arg.nrlks0 = int(optionalParameter)

    optionalParameter = get_content(xmlfile, 'component.property.value', 'insar.number of azimuth looks when forming interferogram.')
    if optionalParameter != None:
        insar_arg.nalks0 = int(optionalParameter)

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

    optionalParameter = get_content(xmlfile, 'component.property.value', 'insar.number of range looks for ionospheric correction.')
    if optionalParameter != None:
        insar_arg.nrlksIon = int(optionalParameter)

    optionalParameter = get_content(xmlfile, 'component.property.value', 'insar.number of azimuth looks for ionospheric correction.')
    if optionalParameter != None:
        insar_arg.nalksIon = int(optionalParameter)

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
    for attr, value in sorted(vars(insar_arg).items()):
        print("{}: {}".format(attr, value))
    print('=====================================\n')


    return insar_arg


def set_filename(insar_arg):
    import glob

    fn = Dummy()

    #basic names
    leader_files = sorted(glob.glob(insar_arg.masterDir + "/" + "LED-ALOS2*-*-*"))
    fn.master = leader_files[0].split('-')[-2]
    leader_files = sorted(glob.glob(insar_arg.slaveDir + "/" + "LED-ALOS2*-*-*"))
    fn.slave = leader_files[0].split('-')[-2]
    ms = fn.master + '-' + fn.slave
    ml = '_{}rlks_{}alks'.format(insar_arg.nrlks, insar_arg.nalks) #do not use str(), which will change the variable type.

    #names from settings
    fn.dem = insar_arg.dem
    fn.demGeo = insar_arg.demGeo

    #output file names
    fn.masterSlc = fn.master + '.slc'
    fn.slaveSlc = fn.slave + '.slc'
    #fn.resampleSlaveSlc = fn.slave + '_resamp.slc'
    fn.masterFrame = fn.master + '.slc.pck'
    fn.slaveFrame = fn.slave + '.slc.pck'
    fn.masterParameter = fn.master + '.slc.par.xml'
    fn.slaveParameter = fn.slave + '.slc.par.xml'

    fn.basline = ms + '_baseline' + '.xml'
    fn.interferogram = ms + '.int'
    fn.amplitude = ms + '.amp'
    #fn.coherence = ms + '.cor' #actually not used, as we will calculate coherence after taking looks
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
    for attr, value in sorted(vars(fn).items()):
        print("{}: {}".format(attr, value))
    print('=====================================\n')


    return fn


