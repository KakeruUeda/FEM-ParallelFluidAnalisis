# Fluid Data Assimilation
## Overview
This code optimizes MRI velociy field using adjoint variable method. <br>
FEM &nbsp; &nbsp; &nbsp;・・・ Discretization in orthogonal grid <br>
XFEM &nbsp; &nbsp;・・・ Interface discretization <br>
Petsc &nbsp; &nbsp;・・・ Solver for matrix system <br>
MPICH &nbsp;・・・ Parallel algorithm depending on petsc
## Usage
    * sh build.sh
    * cd /<example_dir>
    * mpirun -n <process> ./<solver_dir>/<solver_name> test.tp petsc_options.dat
## Features
・SteadyStokesSolver: Solve Stokes equation <br>
・SteadyNavierStokesSolver: Solve Navier Stokes equation <br>
・BCsOptimization: Not yet implemented <br>
## Reference
## Author
Kakeru Ueda
 
