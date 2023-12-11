## --- Fluid Data Assimilation ---
## Overview
This code optimizes MRI velociy field using adjoint variable method. <br>
・FEM <br>
・XFEM (interface) <br>
・Petsc <br>
・MPI
## Usage
    * sh build.sh
    * cd /<example_dir>
    * mpirun -n <process> ./<solver_dir>/<solver> test.tp
## Features
・SteadyStokesSolver: Solve stokes equation <br>
・SteadyNavierStokesSolver: Solve navier stokes equation <br>
・BCsOptimization: Not yet implemented <br>
## Reference
## Author
Kakeru Ueda
 
