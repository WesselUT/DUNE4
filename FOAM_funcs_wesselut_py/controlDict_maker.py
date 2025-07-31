#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: wesselut
"""

import subprocess

def make_controlDict(T_end_hyd):
# make controlDict i
    bmd_file = open("system/controlDict", "w") # overwrite existing fvOptions
    
    controlDict_1 = r"""/*--------------------------------*- C++ -*----------------------------------*\ 
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
    object      controlDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

application     simpleFoam;
// application     pimpleFoam;

startFrom       startTime;

startTime       0;

stopAt          endTime;
"""
     
    bmd_file.write(controlDict_1)
    controlDict_T_end = """
endTime         """ + str(T_end_hyd) + """; """ 

    bmd_file.write(controlDict_T_end)

    controlDict_2 = """
    
deltaT  1;   

//writeControl    adjustable;
writeControl	timeStep;

writeInterval   1;

purgeWrite	1;

writeFormat     ascii;

writePrecision  6;

writeCompression off;

timeFormat      general;

timePrecision   6;

runTimeModifiable true;

adjustTimeStep  on;

maxCo           0.8;

tolerance	1e-8;

functions
{
	/*
	shearStress1
	{
		type	shearStress;
		libs	(fieldFunctionObjects);
		   // Optional (inherited) entries
		    writePrecision  6;
		    writeToFile     true;
		    writeControl    timeStep;
	    	    writeInterval   1;
		    purgeWrite	    1;
	}
	*/
	
	#includeFunc flowRatePatch(name=outlet)
	wallShearStress1
	{
	    // Mandatory entries (unmodifiable)
	    type            wallShearStress;
    	    libs            (fieldFunctionObjects);
    	    
	    // Optional (inherited) entries
	    writePrecision  6;
	    writeToFile     true;
	    writeControl    timeStep;
    	    writeInterval   1;
	    purgeWrite	    1;
	    patches	    ("fixedWalls");
	}
}


// ************************************************************************* //

"""
    
    bmd_file.write(controlDict_2)
    

if __name__ == "__main__":
    # did not manage to write a proper test for this module, as it must be run inside the OF folder
    test = 0
    
   
