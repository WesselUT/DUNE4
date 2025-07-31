#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 13:54:24 2025

@author: wesselut
"""
import pickle
import numpy as np

def save(filename, *args):
    # Get global dictionary
    glob = globals()
    d = {}
    for v in args:
        # Copy over desired values
        d[v] = glob[v]
    with open(filename, 'wb') as f:
        # Put them in the file 
        pickle.dump(d, f)

def load(filename):
    # Get global dictionary
    glob = globals()
    with open(filename, 'rb') as f:
        for k, v in pickle.load(f).items():
            # Set each global variable to the value from the file
            glob[k] = v
            
def save_to_file(fn, var):
    file = open(fn, 'a') # 'a' to append to file
    np.savetxt(file, var, newline=" ")
    file.write("\n")
    file.close()