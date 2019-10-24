# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 15:54:21 2019
@author: spaceparticle
"""
import xml.dom.minidom
import numpy as np
from pandas import DataFrame
from pandas import Series
import base64
import struct


def convBase64LEToFloat(encoded):
    # sample: # encoded = 'pC/TNplSEDlD0II5+nyPOR0doDnZ6tU5FX38OUawDTrQ5iM6o+kvOg6qNzofo1E6OKRtOvnfcjq133s6TtCJOlY8kDos9Zk6MI9VNj/lgzimGhc5ayKDOfvhkjm7H545wlHZOQ==' # gliu: particle time
    # decode the string
    data = base64.standard_b64decode(encoded)
    # ensure that there's enough data for 32-bit floats
    assert len(data) % 4 == 0
    # determine how many floats there are
    count = len(data) // 4
    # unpack the data as floats
    result = struct.unpack('<{0}f'.format(count),  # one big structure of `count` floats
                           data)  # results returned as a tuple
    # print result, '\n count: ', len(result)
    return result


# def getElmntData(elmnt, tagname):
#     return elmnt.getElementsByTagName('tagname')[0].childNodes[0].data


def getElmntInfoInDF(elmnt, attrs):
    DF_Elmnts = DataFrame(columns=attrs)
    for i, elmnt_i in enumerate(elmnt):
        # print 'id: %s, name: %s' % (ptItem.getAttribute('id'), ptItem.getAttribute('name'))
        DF_Elmnts.loc[i] = {tag: elmnt_i.getAttribute(tag) for tag in attrs}
    return DF_Elmnts


def getItemNumByItemName(DF_Items, itemname):
    # print 'got Item Num: ', np.argmax(DF_Items['name']==itemname)
    if itemname in list(DF_Items['name']):
        return np.argmax(DF_Items['name'] == itemname)  # 'Particle X Position'
    else:
        return None


def getElmntByItemNum(elmnt, itemnum):
    for i, elmnt_i in enumerate(elmnt):
        if int(elmnt_i.getAttribute('item')) == itemnum:
            # print 'elmnt_i, elmnt_i.childNodes[0].data \n', elmnt_i, '\n', elmnt_i.childNodes[0].data
            return elmnt_i, elmnt_i.childNodes[0].data
    return


def getParDataByParamName(secDatas, DF_secDatas, paramname, datalength, unit=''):
    param_item_num = getItemNumByItemName(DF_secDatas, paramname)
    assert param_item_num is not None, 'Error! param_item_num is None for paramname {0}'.format(paramname) + \
                                        '. One Possible Cause: This param does not exist in source XML file, pls check the file .'
    if unit != '':
        assert DF_secDatas.loc[param_item_num]['units'] == unit, 'Error, unit does not match for param {0}, '.format(paramname) \
                                                            + unit + ',' + DF_secDatas.loc[param_item_num]['units']
    elmnt_data, data_vals = getElmntByItemNum(secDatas, param_item_num)
    assert data_vals is not None, 'Error! data_vals is None, paramname and param_item_num are ' + paramname + ' ' + str(
        param_item_num)
    data_vals = data_vals.split('\n')
    assert len(data_vals) == 3, 'Error! "Abnormal" structure of raw data of <data_vals>'  
    # ^ raw data is like: u'\n0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 \n      '
    data_vals = data_vals[1].strip()
    data_conv = None
    if elmnt_data.getAttribute('dataFormat') == 'ASCII':
        if DF_secDatas.loc[param_item_num]['type'] == 'INTEGER32':
            data_conv = np.asarray(data_vals.split(' '), dtype=np.int)
        elif DF_secDatas.loc[param_item_num]['type'] == 'FLOAT':
            data_conv = np.asarray(data_vals.split(' '), dtype=np.float)
    elif elmnt_data.getAttribute('dataFormat') == r'Base64/LE':
        #        print 'Base64/LE here'
        if DF_secDatas.loc[param_item_num]['type'] == 'FLOAT':
            data_conv = np.asarray(convBase64LEToFloat(data_vals), dtype=np.float)
        else:
            print 'CAUTION, this case is unexpected'
            data_conv = np.asarray(convBase64LEToFloat(data_vals))
    else:
        data_conv = np.asarray(data_vals.split(' '))
    assert len(data_conv) == datalength or len(
        data_conv) == 1, 'Error!, length does not match ' + paramname + ',' + str(len(data_conv)) + ',' + str(
        datalength)
    #    print '\n', paramname, '\n', data_conv
    return elmnt_data, data_conv