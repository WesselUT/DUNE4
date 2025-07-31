#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 17:17:33 2024

@author: wesselut
"""
import os
import numpy as np

def reshape_OF(ngridx,ngridy,vec):
    #architecture of OF vectors: 
        # [(xleft, ylowest), (xleft+1, ylowest)] etc.
    mat = np.flipud(np.reshape(vec, shape=(ngridy,ngridx)))
    return mat 

def magicfloat(var):
    temp = list(var)
    temp = [char for char in temp if char.isdigit() or char == '.']
    var = "".join(temp)
    return float(var)

def new_folder_number(path):
    folders = os.listdir(path)   
    folder_numbers = []
    for it_folders in folders:
        folder_numbers.append(magicfloat(it_folders))    
    newfolder = str(int(max(folder_numbers)) + 1)
    return newfolder

def readMomSource(fn):
    momSourceFile = open(fn)
    content = momSourceFile.readlines()
    momSource = magicfloat(content[18])
    momSourceFile.close()
    return momSource

def readFlowRate(fn):
    flowRateFile = open(fn)
    content = flowRateFile.readlines()
    flowRate = magicfloat(content[24])
    flowRateFile.close()  
    return flowRate

def get_forcingTerms(path): 
    # this function retrieves the forcing terms over time
    
    # get time directories
    dir_list = os.listdir(path)
    time_list = []

    for directory in dir_list:
        try:
            float(directory)
            time_list.append(directory)
        except:
            pass
    time_list.sort(key=float)
    time_list=np.array(time_list)
    
    time_series_momSource = np.empty(0)
    for timename in time_list[1:len(time_list)]:
         momSource = readMomSource('./'+timename+'/uniform/momentumSourceProperties')
         time_series_momSource = np.append(time_series_momSource, momSource)
    
    return time_series_momSource, time_list