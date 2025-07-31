#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""
import subprocess

def restore_fromLast (LastTime_str, sea, ebb):
    if (sea):
        if (ebb):
            LastTime_str = "ebb/" + LastTime_str
        else: # flood
            LastTime_str = "flood/" + LastTime_str
            
    copy_from = LastTime_str + "/epsilon " + LastTime_str + "/k " + LastTime_str + "/nut " + LastTime_str + "/omega " + LastTime_str + "/p " + LastTime_str + "/U "
    copy_to = "temp_IC" # folder to temporarily store final flow results to be used as primer for next timestep
    makedir_copyTo_command = "mkdir -p " + copy_to;
    subprocess.run(makedir_copyTo_command, shell=True)
    copy_command = "cp " + copy_from + copy_to + " 2>/dev/null || :"
    subprocess.run(copy_command, shell=True)
    
    if (sea):
        # remove flood or ebb folder
        remove_oldFolder = "rm -r " + LastTime_str
        subprocess.run(remove_oldFolder, shell=True)
        # restore acceleration terms. Done for all forcing types - fvOptions is leading in this regard
        if (ebb):
            subprocess.run(['cp'], ['ebb/acceleration-terms.dat'], ['constant/'])
        elif (not ebb): 
            subprocess.run(['cp'], ['flood/acceleration-terms.dat'], ['constant/'])
        
    subprocess.run("./Allclean")
    copy_from = copy_to + "/epsilon " + copy_to + "/k " + copy_to + "/nut " + copy_to + "/omega " + copy_to + "/p " + copy_to + "/U "
    copy_to = "0"
    makedir_0_command = "mkdir -p " + copy_to
    subprocess.run(makedir_0_command, shell=True)
    copy_command = "cp " + copy_from + copy_to + " 2>/dev/null || :"
    subprocess.run(copy_command, shell=True)

