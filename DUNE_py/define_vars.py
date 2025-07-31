#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: wesselut
"""
import numpy as np
import random 

class struct():
    pass

def define_vars_child():
    v = struct()
    
    # define case type
    v.river = True
    v.sea = False
    v.ebb = True # the value of this variable only matters when v.sea = True
    v.initPerturbedBed = True # as of yet there is no other option to define the initial topography
    v.forcing_type = "momForcing"
    # forcing options: 
    #   momForcing (initially choosing velocity, keeping momentum forcing constant)
    #   init_momForcing (initially choosing momentum forcing, keeping it constant)
    #   velocityForcing (keeping velocity constant)
    v.Urep_calc_bool = True
    # if Urep_calc_bool is false, use the following representative velocities:
    v.Urep_ebb_given = 1
    v.Urep_fl_given = -1
    
    v.tau_crit_bool = False
    
    v.save_morphDat = True
    v.save_hydDat = False
    v.save_hydDat_interval = 1 # interval at which hydDat is saved: as this entail quite some data it is recommended to choose this appropriately
    
    # define numerical and environmental variables
    v.tolerance = 0.0015 # heuristically determined timestep control
    
    v.nx = 25
    v.nz = 70
    v.Nt = 20
    v.dt_morph = 0.0 # necessary to initialize
    v.L = 50
    v.x = np.linspace(0,v.L,v.nx+1)
    v.dx = v.x[1]-v.x[0]
    v.x_half = np.linspace(v.dx/2,v.L-v.dx/2,v.nx)
    v.T_end_hyd_afterInit = 500
    v.T_end_hyd_init = 15000
    
    v.rho = 1000
    v.xp = 2*np.pi*v.x/v.L
    v.xp_half = 2*np.pi*v.x_half/v.L
    
    v.depth = 7
    v.Ures = 0.8
    v.UM2 = 1.0 # only relevant in tidal case
    v.domain_width = 0.1
    
    v.vertGrading = 50
    v.h_amp = 0.3 # initial bed amplitude
    v.nConstituents = 1 # highest mode in initially perturbed bed (26 is standard for nonlinear runs)
    v.nstep = 1
    
    v.p_or = 0.4 
    v.a_b = 1.56e-5
    v.beta = 1.5
    v.a_bs = 2.5
    
    v.theta_rad = 30 * np.pi / 180    
    v.lambda_asymptotic = v.a_bs * (1- np.tan(v.theta_rad))
    v.toosteep = 0.466
    
    rho_s = 2650
    nu_w = 1e-6
    g = 9.81
    D50 = 0.35e-3
	
    D_star = (D50 * ((rho_s/v.rho - 1) * g / (nu_w**2.0))**(1.0/3.0) )
    if (D_star < 4):
        v.theta_crit = 0.24 * (D_star**(-1.0))
    elif (D_star < 10):
        v.theta_crit = 0.14 * (D_star**(-0.64))	
    elif (D_star < 20):
        v.theta_crit = 0.04 * (D_star**(-0.1))	
    elif (D_star < 150):    
        v.theta_crit = 0.013 * (D_star**0.29);	
    elif (D_star >= 150):
        v.theta_crit = 0.055
	
    v.tau_crit = v.theta_crit * D50 * v.rho * g * (rho_s/v.rho - 1)
    v.a_b = 0.5 * D50 / (v.tau_crit * (v.rho**(1.0/2.0))) * (D_star**-0.3)
    
    if (not v.tau_crit_bool):
        v.tau_crit = 0
			
    
    v.lambda_Taylor_coefs = [-1.732050807568877, 7.0, -36.373066958946424, 249.0, -2156.403255423252, 22455.0, -272252.4061877140, 3770865.0,
                             -58781967.91634893, 1018232775.0, -19399939902.55034, 403210355625.0, -9078950687045.4, 
                             220154258899575.0, -5719775428751049.0] #oddly enough, the (un?)even coefficients of the Taylor expansion are integers in this case
    return v

def initial_topography(v):
    if (v.initPerturbedBed):
        v.h = np.zeros(shape=(v.nx+1,))
        v.dhdx = v.h
        v.dhdx_half = np.zeros(shape=(v.nx,))
        v.h_half = np.zeros(shape=(v.nx,))
        random.seed(a=0, version=2)
        
        range_forLoop = range(1,v.nConstituents+1, v.nstep)
        length_forloop = len(range_forLoop)
    
        for _ in range_forLoop:
            phase = random.random() * 2 * np.pi
            #phase = 0 ## EDIT
            v.h = v.h + v.h_amp/length_forloop * np.cos(v.xp*_ - phase)
            v.dhdx = v.dhdx - v.h_amp/length_forloop * 2 *np.pi*_/v.L * np.sin(v.xp*_ - phase)
            v.dhdx_half = v.dhdx_half - v.h_amp/length_forloop * 2 *np.pi*_/v.L * np.sin(v.xp_half*_ - phase)
            v.h_half = v.h_half + v.h_amp/length_forloop * np.cos(v.xp_half*_ - phase)
        
        return v
    #elif could add different options here for initial topography, such as loading from file or some other mathematical definition
 

if __name__ == "__main__":
    v = define_vars_child()
    print(v.lambda_Taylor_coefs)       
    #dhdx_test = (np.roll(v.h,-1) - np.roll(v.h,1))/(2*v.dx)

    #plt.plot(v.dhdx)
    #plt.plot(dhdx_test)