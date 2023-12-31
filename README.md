# FEM Parallel Fluid Analisis
## Overview
This code includes CFD and Optimization solvers based on Message Passing Interface. <br>
The equations are discretized in orthogonal grid using Finite Element Method. <br>
Users can apply either eXtended Finite Element Method or Darcy's law to boundary interface.
## Dependencies
・metis: Domain partitioning (these sub domains are owned by each MPI process) <br>
・PETSc: MPI based library for soling matrix system <br>
・TextParser: Enables us to handle text parameter
## Usage
    * sh build.sh
    * cd /<example_dir>
    * mpirun -n <process> ./<solver_dir>/<solver_name> test.tp petsc_options.dat
## Features
・SteadyStokesSolver: Steady Stokes Equation <br>
・SteadyNavierStokesSolver: Steady Navier-Stokes Equation <br>
・UnsteadyNavierStokesSolver: Unsteady Navier-Stokes Equation <br>
・3DVar: 3D Varialional Data Assimilation (Not yet implemented) <br>
・4DVar: 4D Varialional Data Assimilation (Not yet implemented) <br>
## Reference
## Author
Kakeru Ueda
 
