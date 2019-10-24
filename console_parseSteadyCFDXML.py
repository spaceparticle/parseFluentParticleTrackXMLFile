# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 15:54:21 2019
@author: spaceparticle
"""
from parseFluentXML import *
import xml.dom.minidom
import numpy as np
from pandas import DataFrame
from pandas import Series

if __name__ == '__main__':
    # %% ------------- read the XML file and parse information of elements of "Items" and "Section" ----------------------
    XML_FILE = 'en1.xml'  # 'pars_a.xml' #'en1.xml'
    #
    lTransient = False  # False, DO NOT Modify this
    '''
    In Transient simulations of Fluent, particles' data are output (in the src .xml file) at discrete "time moments", 
    while in Steady cases, particles' data are evolved over time.
    '''
    DT = xml.dom.minidom.parse(XML_FILE)
    Con = DT.documentElement
    #
    # The xml file is composed of mainly Two parts, 1st is the definition of each 'item', item refers to physical variable,
    # ... e.g. "Injection","Particle ID","Particle Time" etc. The 2nd part contains values of items (variables).
    ptItems = Con.getElementsByTagName("Item") # 'pt'---point (particle)
    ptSections = Con.getElementsByTagName("Section")
    #
    ptItemClmns = ['id', 'name', 'type', 'units'] # 'columns'(keys) of each Item record, keep consistent with xml file
    ptSecClmns = ['id', 'length'] # 'columns'(keys) of each Section record, keep consistent with xml file
    # ---- Items part ---
    DF_Items = getElmntInfoInDF(ptItems, ptItemClmns) # 'Items' part are parsed and stored as a DataFrame
    DF_Items['id'] = DF_Items['id'].values.astype(np.int) # type conversion
    DF_Items.set_index('id', drop=True, inplace=True) # set index on 'id'
    print 'DF_Items:\n', DF_Items, '\n', DF_Items.index
    ''' e.g.
                               name       type    units
    id                                                 
    0                     Injection     OPTION         
    1                   Particle ID  INTEGER32         
    2                        Region     OPTION         
    3                 Periodic Side  INTEGER32         
    4                 Particle Time      FLOAT        s
    5           Particle X Position      FLOAT        m
    ...
    '''
    # ---- Sections part ---
    DF_Sections = getElmntInfoInDF(ptSections, ptSecClmns)
    DF_Sections.set_index('id', drop=True, inplace=True)
    print 'DF_Sections:\n', DF_Sections, '\n', DF_Sections.index
    ''' e.g.
       length
    id       
    0   10000
    1   10000
    ...
    '''
    # Note that DF_Sections is not used elsewhere.
    # %% ----------------------- read, parse and store data for each particle ------------------------------
    """
    DF structure for record of a specific particle is like:
               xpos ypos zpos vel parRe
    time 1
    time 2
    time ...
    ....
    The overall DF for all particles should be Three-dimensional. BUT as DF is 2D, 
    a Series values of which are list of DFs is used, with particles' id as index.
    """
    pars = Series([])  # The Series
    pardiams = Series([])
    ptSecDataClmns = ['item', 'dataFormat'] # keep consistent with xml file
    count_Sections = len(ptSections)
    for i, ptSec in enumerate(ptSections):
        print '... Parsing Sec {0} of total {1} (0-{2}):'.format(i, count_Sections, count_Sections-1)
        datalength = int(ptSec.getAttribute('length'))
        secDatas = ptSec.getElementsByTagName('Data')
        DF_secDatas = getElmntInfoInDF(secDatas, ptSecDataClmns)
        # 'ParamName' (e.g. 'Particle ID') in following lines is set according to the 'Items' part of the xml file.
        _, parids = getParDataByParamName(secDatas, DF_Items, 'Particle ID', datalength)
        _, partimes = getParDataByParamName(secDatas, DF_Items, 'Particle Time', datalength, unit='s')
        _, xpars = getParDataByParamName(secDatas, DF_Items, 'Particle X Position', datalength, unit='m')
        _, ypars = getParDataByParamName(secDatas, DF_Items, 'Particle Y Position', datalength, unit='m')
        _, zpars = getParDataByParamName(secDatas, DF_Items, 'Particle Z Position', datalength, unit='m')
        _, pardps = getParDataByParamName(secDatas, DF_Items, 'Particle Diameter', datalength, unit='m')
        _, parvels = getParDataByParamName(secDatas, DF_Items, 'Particle Velocity Magnitude', datalength, unit='m s^-1')
        _, parRes = getParDataByParamName(secDatas, DF_Items, 'Particle Reynolds Number', datalength)
        #
        assert len(zpars) == 1, 'Error, the program currently only support 2D problems, i.e. len(zpars) should be exactly 1'
        if len(pardps) == 1:
            pardps = np.ones((datalength,)) * pardps[0]
        if lTransient and len(partimes) == 1:
            partimes = np.ones((datalength,)) * partimes[0]
        if len(parids) == 1:
            parids = np.ones((datalength,), dtype=np.int) * parids[0]
        #
        for ipar in np.unique(parids):
            print 'ipar: ', ipar
            idx = parids == ipar
            pardiams[ipar] = pardps[idx][0]
            #        print 'partimes[idx]\n', partimes[idx]
            data = zip(partimes[idx], xpars[idx], ypars[idx], parvels[idx], parRes[idx], pardps[idx])
            parclmns = ['partime', 'xpar', 'ypar', 'parvel', 'parRe', 'pardp'] # specify self-defined names for each column of <data>
            if not ipar in pars.index:
                pars[ipar] = DataFrame(columns=parclmns)
            newdf = DataFrame(data=data, columns=parclmns)
            pars[ipar] = pars[ipar].append(newdf, ignore_index=True)
    #
    # %% -----------------------group particle diameters-------------------
    dps = np.unique(pardiams.values)
    print 'category of particle diameters: \n', dps
    # dps may look like this '[  9.99999997e-07   1.00000000e-06   2.99999992e-05   3.00000000e-05]' while we expect [1e-7, 1e-6, 3e-5...]
    # so, we do round operation:
    dps_rounded = np.unique([1. / round((1. / dp), 0) for dp in dps])
    print 'category of particle diameters rounded to certain precision: \n', dps_rounded
    #
    # %% ------------------- do some post-calculation or statistics------------
    parRe_mean_groupbypardp = Series([])
    err_rel = 0.001 # due to limited precision of float numbers, particle diameters within certain tolerance (e.g. +-0.1%) are considered having the same size
    for dp in dps_rounded:
        ipars = list(pardiams[(pardiams.values >= dp * (1 - err_rel)) & (pardiams.values <= dp * (1 + err_rel))].index)
        # print 'ipars:\n', ipars
        reSum = 0.
        for ipar in ipars:
            time_ = pars[ipar]['partime'].values[-1] - pars[ipar]['partime'].values[0]
            reSum += np.sum(pars[ipar]['parRe'][:-1] * np.diff(pars[ipar]['partime'])) / time_ # weighted by 'time span'
        reMean = reSum / len(ipars)
        parRe_mean_groupbypardp[dp] = reMean
    print 'parRe_mean_groupbypardp \n', parRe_mean_groupbypardp
    #
    flu_niu = 1.5e-5
    rel_vel = parRe_mean_groupbypardp * flu_niu / parRe_mean_groupbypardp.index
    print 'relative velocity for each particle size are: \n', rel_vel
