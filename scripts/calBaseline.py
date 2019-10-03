#!/usr/bin/env python3

import os
import sys
import pickle
import argparse

import isce
import isceobj
import mroipac
from mroipac.baseline.Baseline import Baseline


def calBaseline(masterFrame, slaveFrame, baselineName):

    #setting up catalog name
    masterYear = str(masterFrame.sensingStart.year)
    masterMonth = '%02d' % masterFrame.sensingStart.month
    masterDay = '%02d' % masterFrame.sensingStart.day

    slaveYear = str(slaveFrame.sensingStart.year)
    slaveMonth = '%02d' % slaveFrame.sensingStart.month
    slaveDay = '%02d' % slaveFrame.sensingStart.day

    catalogName = masterYear + masterMonth + masterDay + '-' + slaveYear + slaveMonth + slaveDay + '_baseline'
    if baselineName != None:
        if baselineName.endswith('.xml'):
            catalogName = os.path.splitext(baselineName)[0]
        else:
            catalogName = baselineName


    catalog = isceobj.Catalog.createCatalog(catalogName)

    #calculate baseline
    optlist = ['all', 'top', 'middle', 'bottom']
    success=False
    baseLocation = None

    for option in optlist:
        baseObj = Baseline()
        baseObj.configure()
        baseObj.baselineLocation = option
        baseObj.wireInputPort(name='masterFrame',object=masterFrame)
        baseObj.wireInputPort(name='slaveFrame',object=slaveFrame)
        try:
            baseObj.baseline()
            success=True
            baseLocation=option
        except:
            print('Baseline computation with option {0} Failed'.format(option))
            pass
        
        if success:
            break

    if not success:
        raise Exception('Baseline computation failed with all possible options. Images may not overlap.')

    print('horizontal_baseline_top: {}'.format(baseObj.hBaselineTop))
    print('horizontal_baseline_rate: {}'.format(baseObj.hBaselineRate))
    print('horizontal_baseline_acc: {}'.format(baseObj.hBaselineAcc))
    print('vertical_baseline_top: {}'.format(baseObj.vBaselineTop))
    print('vertical_baseline_rate: {}'.format(baseObj.vBaselineRate))
    print('vertical_baseline_acc: {}'.format(baseObj.vBaselineAcc))
    print('perp_baseline_top: {}'.format(baseObj.pBaselineTop))
    print('perp_baseline_bottom: {}'.format(baseObj.pBaselineBottom))
    print('baseline_location: {}'.format(baseLocation))

    catalog.addItem('horizontal_baseline_top', baseObj.hBaselineTop, 'baseline')
    catalog.addItem('horizontal_baseline_rate', baseObj.hBaselineRate, 'baseline')
    catalog.addItem('horizontal_baseline_acc', baseObj.hBaselineAcc, 'baseline')
    catalog.addItem('vertical_baseline_top', baseObj.vBaselineTop, 'baseline')
    catalog.addItem('vertical_baseline_rate', baseObj.vBaselineRate, 'baseline')
    catalog.addItem('vertical_baseline_acc', baseObj.vBaselineAcc, 'baseline')
    catalog.addItem('perp_baseline_top', baseObj.pBaselineTop, 'baseline')
    catalog.addItem('perp_baseline_bottom', baseObj.pBaselineBottom, 'baseline')
    catalog.addItem('baseline_location', baseLocation, 'baseline') 

    catalog.renderXml()


def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser( description='Calculate Baseline')
    parser.add_argument('-m', '--master', dest='master', type=str, required=True,
            help = 'master frame')
    parser.add_argument('-s', '--slave', dest='slave', type=str, required=True,
            help = 'slave frame')
    parser.add_argument('-o','--output', dest='output', type=str, default=None,
            help = '(output) baseline file')

    if len(sys.argv) <= 1:
        print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':

    inps = cmdLineParse()

    with open(inps.master, 'rb') as mf:
        masterFrame = pickle.load(mf)
    with open(inps.slave, 'rb') as sf:
        slaveFrame = pickle.load(sf)

    calBaseline(masterFrame, slaveFrame, inps.output)


#example:
#./calBaseline.py -m 20130927.slc.pck -s 20141211.slc.pck



