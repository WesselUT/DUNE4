#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  7 13:57:29 2025

@author: wesselut
"""

import numpy as np
import math
import define_vars
import scipy
from scipy.interpolate import CubicSpline

def Urep_calc(v):
    if (v.ebb):
        plusminus = +1
    else: #flood
        plusminus = -1
    xi = v.Ures/v.UM2
    Urep = plusminus*( (11*xi**2 + 4) / (3*(np.pi + plusminus*2*np.asin(xi))) * np.sqrt(1 - xi**2) + plusminus * xi*(xi**2 + 3/2) )**(1/3)
    return Urep

def t_star_calc(v):
    xi = v.Ures/v.UM2
    t_star_ebb = 1/(2*np.pi) * (np.pi + 2*np.asin(xi))
    t_star_flood = 1/(2*np.pi) * (np.pi - 2*np.asin(xi))
    return t_star_ebb, t_star_flood

def half_h(h,dx):
    h_half = 0.5*(h + np.roll(h,1))
    dhdx_half = (np.roll(h_half,-1) - h_half) / dx
    return h_half, dhdx_half

def roll_matrix(I,shift):
    I_shift = np.zeros(np.shape(I))
    for i in range(0,np.shape(I)[0]):
        I_shift[i,:] = np.roll(I[i,:],shift)
    return I_shift

def derivative_upWind(y,dx):
    dy = (y-np.roll(y,1)) / dx
    return dy

def derivative_downWind(y,dx):
    dy = (np.roll(y,-1) - y) / dx
    return dy

def lambda_calc(v, dhdx_upwind, dhdx_downwind):
    lambda_downwind = np.zeros(v.nx); lambda_upwind = np.zeros(v.nx)
    for i in range(1,np.size(v.lambda_Taylor_coefs)+1):
        lambda_upwind = lambda_upwind + v.a_bs * 1/math.factorial(i) * v.lambda_Taylor_coefs[i-1] * dhdx_upwind**(i-1)   * v.lambda_sign_base**(i+1)
        lambda_downwind = lambda_downwind + v.a_bs * 1/math.factorial(i) * v.lambda_Taylor_coefs[i-1] * dhdx_downwind**(i-1)   * v.lambda_sign_base**(i+1)
    
    # correct for behavior at large stoss slopes
    lambda_upwind[np.where(dhdx_upwind > v.toosteep)] = - v.lambda_sign_base * v.lambda_asymptotic / dhdx_upwind[np.where(dhdx_upwind > v.toosteep)]
    lambda_downwind[np.where(dhdx_downwind > v.toosteep)] = - v.lambda_sign_base * v.lambda_asymptotic / dhdx_downwind[np.where(dhdx_downwind > v.toosteep)]
    lambda_upwind = -lambda_upwind; lambda_downwind = -lambda_downwind # keeping this cumbersome implementation the same as in previous C++ version
    return lambda_upwind, lambda_downwind

def phi_minMod(r):
    phi = np.max([0, np.min([1,r])])
    return phi

def dhdt_calc(v,tau):
    
    tau_abs_pow = np.abs(tau)**v.beta    
    tau_abs_pow_plusOne = np.roll(tau_abs_pow,-1)
    tau_abs_pow_minOne = np.roll(tau_abs_pow,1)
    
    Heaviside = np.zeros(np.shape(tau))
    Heaviside[np.abs(tau) > v.tau_crit] = 1
    tau_abs_pow_crit = np.abs(tau)**(v.beta-1) * (np.abs(tau)  - v.tau_crit) * Heaviside
    
    h_plusOne = np.roll(v.h_half,-1)
    h_minOne = np.roll(v.h_half,1)
    
    I = np.eye(np.size(v.h_half))
    I_minOne = roll_matrix(I,-1)
    I_plusOne = roll_matrix(I,+1)
    
    dhdx_upwind = derivative_upWind(v.h_half, v.dx)
    dhdx_downwind = derivative_downWind(v.h_half, v.dx)
    
    lambda_upwind, lambda_downwind = lambda_calc(v, dhdx_upwind, dhdx_downwind)
    
    v.lambda_upwind = lambda_upwind
    v.lambda_downwind = lambda_downwind
    
    tau_abs_pow_plusHalf = (tau_abs_pow_plusOne + tau_abs_pow) / 2
    tau_abs_pow_minHalf = (tau_abs_pow_minOne + tau_abs_pow) / 2
    
    if v.river:
        tau_abs_pow_lambda_minHalf = tau_abs_pow_minHalf * lambda_upwind
        tau_abs_pow_lambda_plusHalf = tau_abs_pow_plusHalf * lambda_downwind
    elif v.sea: 
        if v.ebb:
            tau_abs_pow_lambda_minHalf_ebb = tau_abs_pow_minHalf * lambda_upwind
            tau_abs_pow_lambda_plusHalf_ebb =  tau_abs_pow_plusHalf * lambda_downwind
				
            tau_abs_pow_lambda_minHalf  = v.t_star_flood * v.tau_abs_pow_lambda_minHalf_flood  + v.t_star_ebb * tau_abs_pow_lambda_minHalf_ebb
            tau_abs_pow_lambda_plusHalf = v.t_star_flood * v.tau_abs_pow_lambda_plusHalf_flood + v.t_star_ebb * tau_abs_pow_lambda_plusHalf_ebb
        
            q_flow_ebb = v.a_b * tau_abs_pow_crit * (tau/np.abs(tau))
        if not v.ebb:
            v.tau_abs_pow_lambda_minHalf_flood = tau_abs_pow_minHalf * lambda_upwind
            v.tau_abs_pow_lambda_plusHalf_flood =  tau_abs_pow_plusHalf * lambda_downwind
            
            v.q_flow_flood = v.a_b * tau_abs_pow_crit * (tau/np.abs(tau))
            
    if (v.river) or (v.sea and v.ebb): # update bed for each timestep (river) or after each ebb (tidal)
        # determine timestep
        diffusion_speed = v.a_b * np.abs(tau_abs_pow_lambda_plusHalf - tau_abs_pow_lambda_minHalf) / v.dx
        mean_diffusion_speed = np.mean(diffusion_speed)
        v.dt_morph = v.tolerance * v.dx / mean_diffusion_speed    
        print('dt_morph = ', str(v.dt_morph/3600), ' h')
    
        # calculate dhdt
        C1_2 = (2 * v.dx**2 * (1 - v.p_or)) / (v.a_b * v.dt_morph)
        # semi-implicit slope correction: 
        LHS_matrix_1 = C1_2 * I
        tau_abs_pow_lambda_minHalf_column = tau_abs_pow_lambda_minHalf[...,None] # promote to column vector
        tau_abs_pow_lambda_plusHalf_column = tau_abs_pow_lambda_plusHalf[...,None] # promote to column vector
        LHS_matrix_2 = (tau_abs_pow_lambda_plusHalf_column + tau_abs_pow_lambda_minHalf_column) * I
        LHS_matrix_3 = - tau_abs_pow_lambda_minHalf_column * I_minOne
        LHS_matrix_4 = - tau_abs_pow_lambda_plusHalf_column * I_plusOne
        LHS_matrix = LHS_matrix_1 + LHS_matrix_2 + LHS_matrix_3 + LHS_matrix_4
        
        RHS_vector_slope = (C1_2 * v.h_half - v.h_half * (tau_abs_pow_lambda_plusHalf + tau_abs_pow_lambda_minHalf) + 
                            h_minOne * tau_abs_pow_lambda_minHalf + h_plusOne * tau_abs_pow_lambda_plusHalf)

        # calculate dhdt_flow
        if (v.sea and v.ebb):
            qb_flow = v.t_star_ebb * q_flow_ebb + v.t_star_flood * v.q_flow_flood
        if (v.river):
            qb_flow = v.a_b * tau_abs_pow_crit * (tau/np.abs(tau))
        
        v.qb_flow = qb_flow
        
        r_h = (v.h_half - np.roll(v.h_half,1)) / (np.roll(v.h_half,-1) - v.h_half)
        r_q = (qb_flow - np.roll(qb_flow,1)) / (np.roll(qb_flow,-1) - qb_flow)
        r_upw = np.zeros(np.shape(r_h))
        for i in range(0,np.size(r_h)):
            r_upw[i] = np.min([r_h[i], r_q[i]])
                
        qb_flow_minOne = np.roll(qb_flow,1)
        qb_flow_plusOne = np.roll(qb_flow,-1)
        qb_flow_plusTwo = np.roll(qb_flow,-2)
        r_minOne = np.roll(r_upw, 1)
        r_plusOne = np.roll(r_upw, -1)    
        q_star = np.zeros(np.shape(qb_flow))
        
        for i in range(0,np.size(qb_flow)):
            q_R = qb_flow_plusOne[i] - phi_minMod(r_plusOne[i]) * (qb_flow_plusTwo[i] - qb_flow_plusOne[i]) / 2
            q_L = qb_flow[i] + phi_minMod(r_upw[i]) * (qb_flow_plusOne[i] - qb_flow[i]) / 2
            S = 0.5 * (q_L + q_R)
            if (q_L > q_R):
                if (S > 0):
                    q_star[i] = q_L
                else:
                    q_star[i] = q_R
            else:
                if (q_L > 0):
                    q_star[i] = q_L
                elif ((q_L < 0) and (q_R>0)):
                    q_star[i] = 0
                else:
                    q_star[i] = q_R
         
        q_iplusHalf = q_star
        q_iminHalf = np.roll(q_star,1)
        
        q_flow_RHS = - (2 * v.dx / v.a_b) * (q_iplusHalf - q_iminHalf)
        
        RHS_vector = RHS_vector_slope + q_flow_RHS
        
        h_plusOne = np.linalg.solve(LHS_matrix, RHS_vector)
        
        v.dhdt = (h_plusOne - v.h_half) / v.dt_morph 
        v.h_half = h_plusOne # update bed 
        v.h_spline = CubicSpline(np.append(v.x_half, v.x_half[-1]+v.dx), np.append(v.h_half, v.h_half[0]), bc_type='periodic')           
        v.dhdx_half = v.h_spline(v.x_half, 1)
        v.h = v.h_spline(v.x)
        v.dhdx = v.h_spline(v.x,1)
        
        
    return v

def adapt_depth_river(v, i):
    if (i==0):
        v.q_disch = v.u_avg * v.depth
    else: 
        newdepth = v.q_disch / v.u_avg
        print("New depth = " + str(newdepth) + "m")
        v.depth = newdepth   
    return v
    
if __name__ == "__main__":
    initPerturbedBed = True
    v = define_vars.define_vars_child()
    v.lambda_sign_base = 1
    v = define_vars.initial_topography(initPerturbedBed, v)
    tau = np.sin(v.x_half)
    dhdt_calc(v,tau)


