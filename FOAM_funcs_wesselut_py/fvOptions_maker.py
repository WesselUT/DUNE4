#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: wesselut
"""

from FOAM_funcs_wesselut_py.numFuncs import Urep_calc
import subprocess

def make_initial_fvOptions(v):
# make initial fvOptions, either for: velocityForcing or momForcing (NOT init_momForcing) OR init_momForcing
    if (v.forcing_type != "init_momForcing"):
        if (v.sea):
            if (v.Urep_calc_bool):
                Urep = Urep_calc(v)
            else: 
                if(v.ebb):
                    Urep = v.Urep_ebb_given
                elif (not v.ebb):
                    Urep = v.Urep_fl_given
    
            make_fvOptions_init (   Urep )                 	
            print("Urep = " , Urep )	
			
        elif (v.river):
        	make_fvOptions_init (v.Ures)
    
    else: # v.forcing_type == "init_momForcing"
        if (v.sea):
        	if (v.ebb):
        		make_fvOptions_afterInit(v.momentumSource_ebb)       	
        	elif (not v.ebb):
        		make_fvOptions_afterInit(v.momentumSource_flood);        	
        elif (v.river):
        	make_fvOptions_afterInit(v.momentumSource_river);
 
def make_fvOptions_init(Ux):

    bmd_file = open("constant/fvOptions", "w") # overwrite existing fvOptions
    
    foam_preamble = r"""/*--------------------------------*- C++ -*----------------------------------*\ 
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2206                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      fvOptions;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""
    
    bmd_file.write(foam_preamble)
    
    fvOptions_1 = """
momentumSource
{
    type        meanVelocityForce;
    active	    yes;
	
    meanVelocityForceCoeffs
    {
    	    selectionMode   all;
"""

    bmd_file.write(fvOptions_1)

    fvOptions_velocity = """
	    Ubar            (""" + str(Ux) + """ 0 0);
"""
    bmd_file.write(fvOptions_velocity)
    
    fvOptions_2 = """
	    fields          (U);
	    relaxation	    1.0;
    }
}

// ************************************************************************* //
"""

    bmd_file.write(fvOptions_2)

def make_fvOptions_afterInit(momSource):
    # first copy the momentum source fvOptions file to the constant folder
    path_to_src = '../FOAM_funcs_wesselut_py/dicts/fvOptions_notInit'
    path_to_destination = 'constant/fvOptions'    
    command_copy_fvOptions = ['cp', path_to_src, path_to_destination]
    subprocess.call(command_copy_fvOptions)
   	
    # now create acceleration terms file (called by fvOptions)
    bmd_file = open("constant/acceleration-terms.dat", "w") # overwrite existing acceleration-terms.dat file   	
   	
    text_acceleration_dat = """
//     t      acc    omega   dOmegadt  
(
(0.0       ((""" + str(momSource) + """ 0 0) (0 0 0) (0 0 0)))
(100000.0  ((""" + str(momSource) + """ 0 0) (0 0 0) (0 0 0)))
); 
// here: t = time [s] acc = acceleration vector [m/s2] omega = rotation vector [rad/s] dOmegadt = angular acceleration [rad/s2]

""" 
    bmd_file.write(text_acceleration_dat)

if __name__ == "__main__":
    # did not manage to write a proper test for this module, as it must be run inside the OF folder
    test = 0
    
   
