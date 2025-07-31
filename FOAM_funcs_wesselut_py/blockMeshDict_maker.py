#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: wesselut
"""
import numpy as np

def make_blockMeshDict(v):    
    bmd_file = open("system/blockMeshDict", "w") # overwrite existing blockMeshDict
    
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
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scale   1.000000;
"""
    
    bmd_file.write(foam_preamble)
    width = v.domain_width
    vertices_mat = np.array([[0, v.h[0], 0], 
                             [v.L, v.h[0], 0], # using here that v.h[0] = v.h[-1]
                            [v.L, v.depth, 0],
                            [0, v.depth, 0],
                            [0, v.h[0], width],
                            [v.L, v.h[0], width],
                            [v.L, v.depth, width],
                            [0, v.depth, width]])
    foam_vertices_blocks = """
vertices 
(
""" 
    for _ in range(0,np.shape(vertices_mat)[0]):
        foam_vertices_blocks = foam_vertices_blocks + "(" + str(vertices_mat[_,0]) + "   " + str(vertices_mat[_,1]) + "   " + str(vertices_mat[_,2]) + ")\n"

    foam_vertices_blocks = foam_vertices_blocks + """);

blocks
(
    hex (0 1 2 3 4 5 6 7) (""" + str(v.nx) + " " + str(v.nz) + " " + """1) simpleGrading (1 """ + str(v.vertGrading) +""" 1)
);
"""
    
    bmd_file.write(foam_vertices_blocks)
    
    edges = """
edges 
(
     polyLine 0 1 ("""
    for _ in range(0,np.size(v.h_half)):
        edges = edges + "(" + str(v.x_half[_]) + " " + str(v.h_half[_]) + " 0) "
    edges = edges + """)
    polyLine 4 5 ("""
    
    for _ in range(0,np.size(v.h_half)):
        edges = edges + "(" + str(v.x_half[_]) + " " + str(v.h_half[_]) + " " + str(width) + ") "
    edges = edges + ")\n); \n"
    
    bmd_file.write(edges)
    
    boundaries = """
boundary
(

    rigidLid
    {
        type wall;
        faces
        (
            (3 7 6 2)
        );
    }
    
    bed
    {
        type wall;
        faces
        (
            (1 5 4 0)
        );
    }
    inlet
    {
     type   cyclic;
     neighbourPatch outlet;
     faces ((0 4 7 3));
    }
    outlet
    {
     type   cyclic;
     neighbourPatch  inlet;
     faces  ((1 2 6 5));
    }            
    frontAndBack
    {
        type empty;
        faces
        (
            (0 3 2 1)
            (4 5 6 7)
        );
    }
);

mergePatchPairs
(
);

// ************************************************************************* //
"""
    bmd_file.write(boundaries)

    
    
    
    
"""
nx = 50
nz = 70
L = 50
x = np.linspace(0, L, nx)
h = 1 * np.cos(2*np.pi*x/L)
depth = 7
vertGrading = 50
make_blockMeshDict(nx, nz, x, h, L, depth, vertGrading)
"""
    