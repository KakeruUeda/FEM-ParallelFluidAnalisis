# Fluid Data Assimilation
## Overview
This code includes CFD and Optimization solvers based on Message Passing Interface. <br>
The equations are discretized in orthogonal grid using Finite Element Method. <br>
To express discontinuity on the interface, eXtended Finite Element Method is also used.
## Dependencies
・metis: Domain partitioning (these sub domains are owned by each MPI process) <br>
・PETSc: MPI based library for soling matrix system <br>
・TextParser: Enables us to handle text parameter
## Usage
    * sh build.sh
    * cd /<example_dir>
    * mpirun -n <process> ./<solver_dir>/<solver_name> test.tp petsc_options.dat
## Features
・SteadyStokesSolver: Solve Stokes equation <br>
・SteadyNavierStokesSolver: Solve Steady Navier Stokes equation <br>
・UnsteadyNavierStokesSolver: Solve Unsteady Navier Stokes equation <br>
・BCsOptimization: Not yet implemented <br>
## Reference
## Author
Kakeru Ueda
 
