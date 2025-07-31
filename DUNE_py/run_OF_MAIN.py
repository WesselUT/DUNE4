#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: wesselut
"""
import sys
import subprocess
sys.path.append("../")
from FOAM_funcs_wesselut_py.blockMeshDict_maker import make_blockMeshDict
from FOAM_funcs_wesselut_py.fvOptions_maker import make_initial_fvOptions
from FOAM_funcs_wesselut_py.fvOptions_maker import make_fvOptions_afterInit
from FOAM_funcs_wesselut_py.restoreFromLast_funcs import restore_fromLast

import numpy as np
import openfoamparser as ofp
import fluidfoam
import FOAM_funcs_wesselut_py.retrieve_dat as retrieve_dat
from FOAM_funcs_wesselut_py.save_load import save_to_file as save_to_file

import define_vars
from FOAM_funcs_wesselut_py import numFuncs 
from FOAM_funcs_wesselut_py.controlDict_maker import make_controlDict 

# before start, run 'openfoam2206' and 'conda activate foamfun' and finally 'spyder' in command prompt

# file management
path = './'

# get variables in struct v
v = define_vars.define_vars_child()
v = define_vars.initial_topography(v)

path_morphDat_stored = '../morphDat_stored'
path_hydDat_stored = '../hydDat_stored'

subprocess.run(['mkdir','-p', path_morphDat_stored+'/morphDat0']) # ensure there is one dummy folder 

if v.save_hydDat:
    subprocess.run(['mkdir','-p', path_hydDat_stored+'/hydDat0']) # ensure there is one dummy folder 
    new_folder_number_hydDat = retrieve_dat.new_folder_number(path_hydDat_stored)

subprocess.run(['rm','-r','morphDat', 'hydDat']) # delete morphDat and hydDat directories if they exist
subprocess.run(['mkdir', 'morphDat', 'hydDat']) # ... and make new ones

save_to_file('morphDat/xgrid.txt', v.x_half)

if (v.sea):
    subprocess.run(['rm','-r', 'flood', 'ebb'])
    subprocess.run(['mkdir', 'flood', 'ebb'])

x_grid, y_grid, Ux, Uy, p, k, t, momForcing, h_stored, dhdt_stored, tau, depth_stored, morph_time_vec = ([] for i in range(13))

for i in range(0,v.Nt): # run morphodynamic simulation
    print('Starting loop nr ' + str(i+1) + '\n')
    if (v.sea):
        v.ebb = (i % 2 == 1)
        if (v.ebb):
            print('starting ebb flow calculation')
        elif (not v.ebb):
            print('starting flood flow calculation')
            
    if (v.sea and v.ebb) or (v.river):
        v.lambda_sign_base = 1 
    elif (v.sea and not v.ebb):
        v.lambda_sign_base = -1
    
    make_blockMeshDict(v) 
    
    # initialisation
    if (v.river and i==0) or (v.sea and i<=1): #initialisation is first flow solution (river) or first two flow solutions (in case of sea; for both ebb and flood)
        print('Initialising...\n')
        make_controlDict(v.T_end_hyd_init)
        make_initial_fvOptions(v)
        
        # initial flow
        init_flow = "init_flow/" + "rest" # only initial flow at rest is implemented as of yet
        copy_init_flow_command = ("cp "+init_flow +"/k "+init_flow+
        "/nut "+init_flow+"/epsilon "+init_flow+"/omega "+init_flow+
        "/p "+init_flow+"/U 0.orig")
        subprocess.run(copy_init_flow_command, shell=True)
        
        # run OpenFOAM for flow solution at first timestep
        subprocess.run("./Allclean", shell=True)
        subprocess.run("./Allrun_first", shell=True);		
        
        latest_time_str = str(fluidfoam.readof._find_latesttime(path))
        if (v.sea):
            v.t_star_ebb, v.t_star_flood = numFuncs.t_star_calc(v)
            if (v.ebb):
                latest_time_str_ebb = latest_time_str
            else: #flood
                latest_time_str_flood = latest_time_str                
        
        if (v.forcing_type == "momForcing"):
            momForcing_temp, t_temp = retrieve_dat.get_forcingTerms(path)
            momSource = momForcing_temp[-1]            
            if (v.sea and v.ebb) or (v.river):
                make_fvOptions_afterInit(-momSource)
            else:
                make_fvOptions_afterInit(momSource) # for some reason, OpenFOAM doesn't output the momentumsource properties with a consistent sign     
            if (v.sea): # if sea, save the fvOptions file in either ebb or flood directory. Hydro dat is saved in same directory
                src = "constant/acceleration-terms.dat"
                if (v.ebb):
                    dst = "ebb/"
                else: 
                    dst = "flood/"			
                subprocess.run(['mv', src, dst]); 
            
            print('momentumsource: '+str(momSource))
    else: # after initialisation       
        make_controlDict(v.T_end_hyd_afterInit)
        restore_fromLast (latest_time_str, v.sea, v.ebb) # copy from last time to 0, this script includes ./Allclean
        if (v.forcing_type == "velocityForcing"):
            make_initial_fvOptions(v) # keeps velocity constant over simulation by using fvOptions --> velocityForcing
        subprocess.run("./Allrun_fromSecond");	# run OpenFOAM without restoring from 0dir (Allrun_fromSecond)	
        latest_time_str = str(fluidfoam.readof._find_latesttime(path)) # find name of folder with latest 'timestep' of simulation
        if (v.forcing_type == "velocityForcing"):        
            momForcing_temp, t_temp = retrieve_dat.get_forcingTerms(path)
            momSource = momForcing_temp[-1]          
    
    print('Number of hydro timesteps = ' + latest_time_str)
    
    if (v.save_hydDat):
        if (i % v.save_hydDat_interval==0): 
            newFolder_hydDat = path_hydDat_stored + '/hydDat' + new_folder_number_hydDat+'/'+str(i) +'/'
            subprocess.run(['mkdir', '-p', newFolder_hydDat])
            move_hydDat_command = 'cp -r ' + latest_time_str+'/* ' + newFolder_hydDat
            subprocess.run(move_hydDat_command, shell=True)
            
    
    # obtain shear stress
    tau_dict = ofp.parse_boundary_field(path + latest_time_str + '/wallShearStress') 
    tau_temp = tau_dict[b'bed'][b'value']
    
    # obtain flow rate and domain-average velocity
    fn = path + latest_time_str + "/uniform/functionObjects/functionObjectProperties"
    flowRate = retrieve_dat.readFlowRate(fn)
    v.u_avg = flowRate / v.domain_width / v.depth
    if (v.river):
        v = numFuncs.adapt_depth_river(v, i)

    # postprocessing tau
    tau_rho = tau_temp * v.rho # from volumetric to actual shear stress

    tau_x = -tau_rho[:,0]
    tau_z = -tau_rho[:,1]

    tau_x_parallel = tau_x / np.sqrt(1 + v.dhdx_half**2)
    tau_z_parallel = tau_z * v.dhdx_half**2 / np.sqrt(1 + v.dhdx_half**2)

    tau_x = tau_x_parallel + tau_z_parallel

    tau.append(tau_x) # store bed-shear stress
    
    # calculate bed evolution    
    v = numFuncs.dhdt_calc(v,tau_x)
    
    if (i==0):
        morph_time_vec.append(v.dt_morph)
    else: 
        if (v.sea and (not v.ebb)):
            v.dt_morph = 0
        morph_time_vec.append(morph_time_vec[i-1] + v.dt_morph) 
        print("time = " + str(morph_time_vec[i]/3600/24) + "d")
        
    # store variables
    if (v.sea and v.ebb) or (v.river):
        save_to_file('morphDat/h.txt', v.h_half)
        save_to_file('morphDat/dhdt.txt', v.dhdt)
        save_to_file('morphDat/dhdx.txt', v.dhdx_half)
        save_to_file('morphDat/qb.txt', v.qb_flow)
        save_to_file('morphDat/morphTime.txt', np.array([morph_time_vec[i]]))
        save_to_file('morphDat/U.txt', np.array([v.u_avg]))
        save_to_file('morphDat/depth.txt', np.array([v.depth]))
        save_to_file('morphDat/forcing.txt', np.array([momSource]))
    if v.river:
        save_to_file('morphDat/tau_x.txt', tau_x)
    elif v.sea and v.ebb:
        save_to_file('morphDat/tau_x_ebb.txt', tau_x)
    elif v.sea and not v.ebb:
        save_to_file('morphDat/tau_x_flood.txt', tau_x)
if v.save_hydDat:
    print('hydDat saved as hydDat' + new_folder_number_hydDat)

if v.save_morphDat: # move morphDat folder to ../morphDat_stored
    new_folder_number_morphDat = retrieve_dat.new_folder_number(path_morphDat_stored)
    subprocess.run(['mv', 'morphDat', path_morphDat_stored+'/morphDat'+new_folder_number_morphDat])
    print('morphDat stored as morphDat'+new_folder_number_morphDat)

             
