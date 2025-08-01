# DUNE4
Morphodynamic model of underwater dunes in river and tidal environments

------------------------------------------DUNE RECIPE------------------------------------------
1. Install OpenFOAM (version 2206 was used in developing this model) and learn how to use it
2. Install the conda environment 'foamfun' (specified by the file foamfun-env.yaml) - conda installation required
3. Copy the folders DUNE_py and FOAM_funcs_wesselut_py somewhere OpenFOAM can find it
4. Enter 'openfoam2206' (or whichever version you're using) and 'conda activate foamfun' on the command line
5. Go to the folder DUNE_py and enter 'spyder' (which is included in the foamfun environment) on the command line, open 'run_OF_MAIN.py'
6. Press run
7. To adjust parameters of the problem, edit the file 'define_vars.py'
8. Initial conditions and boundary conditions are specified in the files located in 'init_flow/rest', hence must be changed here in case you wish to use different conditions
9. The files fvSchemes and fvSolution (located in 'system/') are not altered by the program and thus must be changed manually if needed (same as for the initial/boundary conditions)

For more information on the model's rationale, please see https://doi.org/10.1016/j.geomorph.2025.109649 (river dunes) and https://doi.org/10.3990/1.9789036560320 (PhD thesis, chapter 5 contains information on the implementation of the tidal setting). 
