# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 15:54:21 2019

@author: adming
"""
import sys
import xml.dom.minidom
import numpy as np
import pandas as pd
from pandas import DataFrame 
from pandas import Series
import base64
import struct

def convBase64LEToFloat(encoded):
#    encoded = 'pC/TNplSEDlD0II5+nyPOR0doDnZ6tU5FX38OUawDTrQ5iM6o+kvOg6qNzofo1E6OKRtOvnfcjq133s6TtCJOlY8kDos9Zk6MI9VNj/lgzimGhc5ayKDOfvhkjm7H545wlHZOQ==' # gliu: particle time
    # decode the string
    data = base64.standard_b64decode(encoded)
    # ensure that there's enough data for 32-bit floats
    assert len(data) % 4 == 0
    # determine how many floats there are
    count = len(data) // 4
    # unpack the data as floats
    result = struct.unpack('<{0}f'.format(count), # one big structure of `count` floats
                           data)                  # results returned as a tuple
#    print result, '\n count: ', len(result)
    return result

def getElmntData(elmnt, tagname):
    return elmnt.getElementsByTagName('tagname')[0].childNodes[0].data


def getElmntInfoInDF(elmnt, attrs):
    DF_Elmnts = DataFrame(columns=attrs)
    for i, elmnt_i in enumerate(elmnt):
#        print 'id: %s, name: %s' % (ptItem.getAttribute('id'), ptItem.getAttribute('name'))
        DF_Elmnts.loc[i] = {tag : elmnt_i.getAttribute(tag) for tag in attrs} 
    return DF_Elmnts

def getItemNumByItemName(DF_Items, itemname):
#    print 'got Item Num: ', np.argmax(DF_Items['name']==itemname)
    return np.argmax(DF_Items['name']==itemname) # 'Particle X Position'

def getElmntByItemNum(elmnt, itemnum):
    for i, elmnt_i in enumerate(elmnt):
        if int(elmnt_i.getAttribute('item')) == itemnum:
#            print 'elmnt_i, elmnt_i.childNodes[0].data \n', elmnt_i, '\n', elmnt_i.childNodes[0].data
            return elmnt_i, elmnt_i.childNodes[0].data
    return

def getParDataByParamName(secDatas, DF_secDatas, paramname, datalength, unit=''):
    param_item_num = getItemNumByItemName(DF_secDatas, paramname)
    assert param_item_num is not None, 'Error! param_item_num is None, paramname is: ' + paramname
    if unit != '':
        assert DF_secDatas.loc[param_item_num]['units'] == unit, 'Error, unit does not match, ' + unit + ',' + DF_secDatas.loc[param_item_num]['units']
    elmnt_data, data_vals = getElmntByItemNum(secDatas, param_item_num)
    assert data_vals is not None, 'Error! data_vals is None, paramname and param_item_num are ' + paramname + ' ' + str(param_item_num)
    data_vals = data_vals.split('\n')
    assert len(data_vals) == 3, 'Error! "Abnormal" structure of raw data of <data_vals>'# raw data is like this: u'\n0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 \n      '
    data_vals = data_vals[1].strip()
    data_conv = None
    if elmnt_data.getAttribute('dataFormat')=='ASCII':
        if DF_secDatas.loc[param_item_num]['type'] == 'INTEGER32' :
            data_conv = np.asarray(data_vals.split(' '), dtype=np.int)
        elif DF_secDatas.loc[param_item_num]['type'] == 'FLOAT':
            data_conv = np.asarray(data_vals.split(' '), dtype=np.float)
    elif elmnt_data.getAttribute('dataFormat')== r'Base64/LE':
#        print 'Base64/LE here'
        if DF_secDatas.loc[param_item_num]['type'] == 'FLOAT':
            data_conv = np.asarray(convBase64LEToFloat(data_vals), dtype=np.float)
        else:
            print 'CAUTION, this case is unexpected'
            data_conv = np.asarray(convBase64LEToFloat(data_vals))
    else:
        data_conv = np.asarray(data_vals.split(' '))
    assert len(data_conv) == datalength or len(data_conv) == 1, 'Error!, length does not match ' + paramname + ',' + str(len(data_conv)) + ',' + str(datalength)
#    print '\n', paramname, '\n', data_conv
    return elmnt_data, data_conv
    
#%% -------------- read the XML file and parse information of elements of "Items" and "Section" ----------------------
XML_FILE = 'en1.xml' #'pars_a.xml' #'en1.xml'
lTransient = False # in Transient simulations, particles' data are output (in the src .xml file) at discrete "time moments", while in Steady cases, particles' data are evolved over time.
# ... 
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
print 'DF_Items:\n', DF_Items, '\n', DF_Items.index

DF_Sections = getElmntInfoInDF(ptSections, ptSecClmns)  
DF_Sections.set_index('id', drop=True, inplace=True)
print 'DF_Sections:\n', DF_Sections, '\n', DF_Sections.index
#%% ---------------- read, parse and store data for each particle-------------------------
"""
DF structure for record of a specific particle:
        xpos ypos zpos vel parRe
time

The overall DF for all particles should be Three-dimensional. 
BUT as DF is 2D, I decide to use a Series, values of which are list of DFs, 
indexes of which are particles' id.
"""
pars = Series([]) # DO NOT use this: 'pars = Series()'
pardiams = Series([])
ptSecDataClmns = ['item', 'dataFormat']
for i, ptSec in enumerate(ptSections):
    print 'Now Sec #:', i
    datalength = int(ptSec.getAttribute('length'))
    secDatas = ptSec.getElementsByTagName('Data')
    DF_secDatas = getElmntInfoInDF(secDatas, ptSecDataClmns)
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
        partimes= np.ones((datalength,)) * partimes[0]
    if len(parids) == 1:
        parids = np.ones((datalength,), dtype=np.int) * parids[0]
#    
    for ipar in np.unique(parids):
        print 'ipar: ', ipar
        idx = parids==ipar
        pardiams[ipar] = pardps[idx][0]
#        print 'partimes[idx]\n', partimes[idx]
        data = zip(partimes[idx], xpars[idx], ypars[idx], parvels[idx], parRes[idx], pardps[idx])
        parclmns = ['partime', 'xpar', 'ypar', 'parvel', 'parRe', 'pardp']        
        if not ipar in pars.index:
            pars[ipar] = DataFrame(columns=parclmns)  # this DOES NOT work
        newdf = DataFrame(data=data, columns=parclmns)
        pars[ipar] = pars[ipar].append(newdf, ignore_index=True) # this will work
#    if i >= 3: break
#        
#%% -----------------------group particles diameters-------------------
dps = np.unique(pardiams.values)
print 'category of particle diameters: \n', dps
# dps may look like this '[  9.99999997e-07   1.00000000e-06   2.99999992e-05   3.00000000e-05]'
dps_rounded = np.unique([1./round((1./dp),0) for dp in dps])
print 'category of particle diameters rounded to certain precision: \n', dps_rounded
#
parRe_mean_groupbypardp = Series([])
err_rel = 0.001
for dp in dps_rounded:
    ipars = list(pardiams[(pardiams.values >= dp * (1-err_rel)) & (pardiams.values <= dp * (1+err_rel))].index)
    print 'ipars:\n',ipars
    reSum = 0.
    for ipar in ipars:
#        reSum += np.mean(pars[ipar]['parRe'])
        time_ = pars[ipar]['partime'].values[-1] - pars[ipar]['partime'].values[0]
        reSum += np.sum(pars[ipar]['parRe'][:-1] * np.diff(pars[ipar]['partime'])) / time_
    reMean = reSum / len(ipars)
    parRe_mean_groupbypardp[dp] = reMean
print 'parRe_mean_groupbypardp \n', parRe_mean_groupbypardp
#%% ------------------- do some post-calculation or statistics------------
flu_niu = 1.5e-5
rel_vel = parRe_mean_groupbypardp * flu_niu / parRe_mean_groupbypardp.index
print 'relative velocity are: \n', rel_vel        
    
    
    
        