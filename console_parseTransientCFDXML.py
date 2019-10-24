# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 15:54:21 2019
@author: spaceparticle
For Transient CFD simulation with Fluent, one may have multiple injections of particles. However, to use this code,
one should make an .xml file (exported by Fluent) for each particle size. An xml file containing data for particles with different
sizes is not supported by this code at this moment (and the author is not planning to go further). Please define an
injection in Fluent for each particle size you are concerned.
"""

from parseFluentXML import *
import xml.dom.minidom
import numpy as np
import pandas as pd
from pandas import DataFrame
import logging
import os
from datetime import datetime

# --- for logging purpose, storing important info into a file. ---
logFileName = 'log_' + datetime.strftime(datetime.now(), '%Y%M%d_%H%M%S') + '.log'
logger = logging.getLogger("main")
logger.handlers = []
logging.basicConfig(level=logging.INFO)
filehandler = logging.FileHandler(logFileName)
filehandler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
filehandler.setFormatter(formatter)
logger.addHandler(filehandler)
#
# console = logging.StreamHandler() # print on screen # # console.setLevel(logging.DEBUG)
# logger.info("this is info")#logger.debug("this is debug")#logger.warning("this is warning")#logger.error("this error")#logger.critical("this is critical")
#

def parseTransientCFDXMLFile(XML_FILE='pars.xml', zones_key=None, zones_xbound=None, flu_niu=None, ldbg=False):
    '''
    This function is not general-purpose as those in parseFluentXML.py. It contains problem-specific codes and one can use
    :param XML_FILE: the name of the source xml file
    :param zones_key: keys of zones_xbound
    :param zones_xbound: like: {'0_5D':[0.0, 5.0], '5_10D':[5.0, 10.0], '10_15D':[10.0, 15.0], '15_20D':[15.0, 20.0], '20_30D':[20.0, 30.0]}
    :param flu_niu: niu of fluid
    :param ldbg: print debug info or not
    :return: No explicit return.
    '''
    # %% -------------- read the XML file and parse information of elements of "Items" and "Section" ----------------------
    print '========= parsing file: {} ============='.format(XML_FILE)
    logger.info('========= parsing file: {} ============='.format(XML_FILE))
    DT = xml.dom.minidom.parse(XML_FILE)
    Con = DT.documentElement
    #
    ptItems = Con.getElementsByTagName("Item")
    ptSections = Con.getElementsByTagName("Section")
    #
    ptItemClmns = ['id', 'name', 'type', 'units']
    ptSecClmns = ['id', 'length']
    #
    DF_Items = getElmntInfoInDF(ptItems, ptItemClmns)
    DF_Items['id'] = DF_Items['id'].values.astype(np.int)
    DF_Items.set_index('id', drop=True, inplace=True)
    if ldbg: print 'DF_Items:\n', DF_Items, '\n', DF_Items.index
    #
    DF_Sections = getElmntInfoInDF(ptSections, ptSecClmns)
    DF_Sections.set_index('id', drop=True, inplace=True)
    if ldbg: print 'DF_Sections:\n', DF_Sections, '\n', DF_Sections.index
    # %% ---------------- read, parse and store data for each particle-------------------------
    # --------- problem-wise --------
    cylD_default = 0.02 # the author was modelling flow past a cylinder in Fluent. cylD is the diameter of the cylinder.
    cylD = float(raw_input('Please input the value of the cylinder diameter, unit is m, default is {}:  '.format(
        cylD_default)) or cylD_default)
    print 'cylD is set to: {} (m)'.format(cylD)
    logger.info('cylD is set to: {} (m)'.format(cylD))
    # ---- end ---
    #
    print 'collecting information of time moments (partimes) and particle sizes (parsizes)...'
    parsizes = np.array([], dtype=np.float)
    partimes = np.array([], dtype=np.float)
    for i, ptSec in enumerate(ptSections): # each ptSec corresponds to data at a single time moment (partime)
        # print 'Now Sec #:', i
        datalength = int(ptSec.getAttribute('length'))
        secDatas = ptSec.getElementsByTagName('Data')
        # DF_secDatas = getElmntInfoInDF(secDatas, ptSecDataClmns)
        _, partime = getParDataByParamName(secDatas, DF_Items, 'Particle Time', datalength, unit='s')
        _, parsize = getParDataByParamName(secDatas, DF_Items, 'Particle Diameter', datalength, unit='m')
        _, zpars = getParDataByParamName(secDatas, DF_Items, 'Particle Z Position', datalength, unit='m')
        # print partime, len(partime)
        # print parsize, len(parsize)
        partimes = np.append(partimes, partime)
        parsizes = np.append(parsizes, parsize)
        assert len(zpars) == 1, 'Error, the program currently only support 2D problems, i.e. len(zpars) should be exactly 1'
    # DO NOT modify (like using 'sort') the orders of 'partimes'
    print '---'
    # print partimes, len(partimes)
    # print parsizes, len(parsizes)
    assert len(np.unique(parsizes)) == 1, 'Error, particles having more than 1 size---Not supported by this program'
    # ^ as mentioned at the head of this file, we do not support multi-sized particles in one xml file.
    partimes_unique = np.unique(partimes)
    time_round_decimals = 5
    parsize_round_decimals = 2
    parsizes = [round(x * 1.e6, parsize_round_decimals) / 1.e6 for x in parsizes]
    par_stats = DataFrame(index=pd.MultiIndex.from_product([zones_xbound.keys(), [str(round(s, time_round_decimals)) \
                                                                                  for s in partimes_unique]],
                                                           names=['zone', 'time']), \
                          columns=['NPars', 'Res_mean', 'Res_min', 'Res_max', 'xmin', 'xmax', 'xmean', 'ymin', 'ymax',
                                   'ymean', 'velmin', 'velmax', 'velmean']) # stats --- statistics
    for ptime_unq in partimes_unique:
        print 'working on t = {} s'.format(ptime_unq)
        idx_Secs = list(np.where(partimes == ptime_unq)[0])
        # print idx
        parids_s = np.array([], dtype=np.float)
        xpars_s = np.array([], dtype=np.float)
        ypars_s = np.array([], dtype=np.float)
        pardps_s = np.array([], dtype=np.float)
        parvels_s = np.array([], dtype=np.float)
        parRes_s = np.array([], dtype=np.float)
        for iSec in idx_Secs:
            if ldbg: print 'Now Sec #:', iSec
            ptSec = ptSections[iSec]
            datalength = int(ptSec.getAttribute('length'))
            secDatas = ptSec.getElementsByTagName('Data')
            _, parids = getParDataByParamName(secDatas, DF_Items, 'Particle ID', datalength)
            _, partime = getParDataByParamName(secDatas, DF_Items, 'Particle Time', datalength, unit='s')
            _, xpars = getParDataByParamName(secDatas, DF_Items, 'Particle X Position', datalength, unit='m')
            _, ypars = getParDataByParamName(secDatas, DF_Items, 'Particle Y Position', datalength, unit='m')
            _, zpars = getParDataByParamName(secDatas, DF_Items, 'Particle Z Position', datalength, unit='m')
            _, pardps = getParDataByParamName(secDatas, DF_Items, 'Particle Diameter', datalength, unit='m')
            _, parvels = getParDataByParamName(secDatas, DF_Items, 'Particle Velocity Magnitude', datalength, unit='m s^-1')
            _, parRes = getParDataByParamName(secDatas, DF_Items, 'Particle Reynolds Number', datalength)
            #
            assert len(partime) == 1 and partime[0] == ptime_unq, 'length != 1 or partime inconsistent'
            assert len(parids) == len(np.unique(parids)), 'parids Not unique'
            parids_s = np.append(parids_s, parids)
            xpars_s = np.append(xpars_s, xpars)
            ypars_s = np.append(ypars_s, ypars)
            pardps_s = np.append(pardps_s, pardps)
            parvels_s = np.append(parvels_s, parvels)
            parRes_s = np.append(parRes_s, parRes)
        pt = str(round(ptime_unq, time_round_decimals))
        for zone_key in zones_xbound.keys():
            ids = (xpars_s >= zones_xbound[zone_key][0] * cylD) & (xpars_s < zones_xbound[zone_key][1] * cylD)
            par_stats.loc[(zone_key, pt)]['NPars'] = len(parids_s[ids])
            par_stats.loc[(zone_key, pt)]['Res_mean', 'Res_min', 'Res_max'] = [np.mean(parRes_s[ids]),
                                                                               np.min(parRes_s[ids]),
                                                                               np.max(parRes_s[ids])]
            par_stats.loc[(zone_key, pt)]['xmin', 'xmax', 'xmean'] = [np.min(xpars_s[ids]), np.max(xpars_s[ids]),
                                                                      np.mean(xpars_s[ids])]
            par_stats.loc[(zone_key, pt)]['ymin', 'ymax', 'ymean'] = [np.min(ypars_s[ids]), np.max(ypars_s[ids]),
                                                                      np.mean(ypars_s[ids])]
            par_stats.loc[(zone_key, pt)]['velmin', 'velmax', 'velmean'] = [np.min(parvels_s[ids]),
                                                                            np.max(parvels_s[ids]),
                                                                            np.mean(parvels_s[ids])]
            # ^ Re --- particles Reynolds number. vel --- particle velocity.
    #
    logger.info('par_stats: \n {}'.format(par_stats))
    assert len(np.unique(pardps_s)) == 1, 'Error, particles diameter is not unique ---Not supported by this program'
    par_diam = pardps_s[0]
    #
    print 'making statistics of particle Re...'
    # making average of each zone over different time moments
    Res_stats = DataFrame(index=zones_xbound.keys(), columns=['mean', 'rel_std'])
    for zone_key in zones_xbound.keys():
        Res_stats.loc[zone_key]['mean', 'rel_std'] = [np.mean(par_stats.loc[zone_key]['Res_mean']),
                                                      np.std(par_stats.loc[zone_key]['Res_mean']) / np.mean(
                                                          par_stats.loc[zone_key]['Res_mean'])]
    #
    Res_stats = Res_stats.reindex(zones_key)
    if np.max(Res_stats['rel_std']) > 0.03:
        print 'Warning: relative std of Res over different time moments is too big as {}'.format(
            np.max(Res_stats['rel_std']))
        # ^ If so, one had better add more time-moments to the output file of each particle size. Time span and number of
        #  time-moments should be big enough such that they are sufficiently representative of the simulated flow.
    print '--> Res_stats: \n', Res_stats
    logger.info('Res_stats: \n {}'.format(Res_stats))
    rel_vel = Res_stats['mean'] * flu_niu / par_diam
    rel_vel = rel_vel.reindex(zones_key)
    print '\n--> relative velocity are: \n', rel_vel
    logger.info('rel_vel: \n {}'.format(rel_vel))
    #
    print '> Done of {}'.format(XML_FILE)


if __name__ == '__main__':
    zones_key = ['0_5D', '5_10D', '10_15D', '15_20D', '20_30D']
    zones_xbound = {'0_5D': [0.0, 5.0], '5_10D': [5.0, 10.0], '10_15D': [10.0, 15.0], '15_20D': [15.0, 20.0],
                    '20_30D': [20.0, 30.0]}
    assert sorted(zones_key) == sorted(zones_xbound.keys()), 'Error, keys of zones_xbound do not match zones_key'
    XML_File_Dir = './2D_x24D_Standard'
    if XML_File_Dir:
        XML_Files = []
        for fl in os.listdir(XML_File_Dir):
            if fl and fl[-4:] == '.xml':
                XML_Files.append(XML_File_Dir + '/' + fl)
    else:
        XML_Files = ['par_1.67um.xml', ]  # , par_5um.xml', 'par_50um.xml')
    logger.info('>>> pwd is: {0}, files to parse: {1}'.format(os.path.dirname(os.path.realpath(__file__)), XML_Files))
    for xml_file in XML_Files:
        parseTransientCFDXMLFile(XML_FILE=xml_file, zones_key=zones_key, zones_xbound=zones_xbound, flu_niu=1.5e-5)
    #
