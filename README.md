# FEM Parallel Fluid Analisis
## Overview
This code includes CFD solvers based on Message Passing Interface. <br>
The equations are discretized in orthogonal grid using Finite Element Method. <br>
## Dependencies
・metis: Domain partitioning (these sub domains are owned by each MPI process) <br>
・PETSc: MPI based library for soling matrix system <br>
・TextParser: Enables us to handle text parameter
## Usage
    * sh build.sh
    * cd /<example_dir>
    * mpirun -n <process> ./<solver_dir>/<solver_name> <tp_name>.tp petsc_options.dat
## Features
・SteadyStokesSolver: Steady Stokes Equation Solver <br>
・SteadyNavierStokesSolver: Steady Navier-Stokes Equation Solver <br>
・UnsteadyNavierStokesSolver: Unsteady Navier-Stokes Equation Solver <br>
## References
https://github.com/oubiomechlab/voxelFEMfluid.git <br>
https://github.com/chennachaos/cfdsfemmpi.git