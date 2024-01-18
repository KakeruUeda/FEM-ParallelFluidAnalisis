# FEM Parallel Fluid Analisis
## Overview
This code includes CFD solvers and fluid optimizer based on Message Passing Interface. <br>
The equations are discretized in orthogonal grid using Finite Element Method. <br>
Users can apply either eXtended Finite Element Method or Darcy's law to wall boundary elements.
## Dependencies
・metis: Domain partitioning (these sub domains are owned by each MPI process) <br>
・PETSc: MPI based library for soling matrix system <br>
・TextParser: Enables us to handle text parameter
## Usage
    * sh build.sh
    * cd /<example_dir>
    * mpirun -n <process> ./<solver_dir>/<solver_name> <tp_name>.tp petsc_options.dat
## Features
・SteadyStokesSolver: Solver for Steady Stokes Equation <br>
・SteadyNavierStokesSolver: Solver for Steady Navier-Stokes Equation <br>
・UnsteadyNavierStokesSolver: Solver for Unsteady Navier-Stokes Equation <br>
・3DVar: 3D Variational Data Assimilation (Not yet implemented) <br>
・4DVar: 4D Variational Data Assimilation (Not yet implemented) <br>
## Reference
## Author
Kakeru Ueda
