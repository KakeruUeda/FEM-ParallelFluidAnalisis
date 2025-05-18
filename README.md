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
## examples

・Flow through a stenosed tube

<p>
  <img src="images/stenosis0_x.png" alt="Stenosis 0" width="500" style="margin-bottom: 10px;"><br>
  <img src="images/stenosis1_x.png" alt="Stenosis 1" width="500" style="margin-bottom: 10px;"><br>
  <img src="images/stenosis2_x.png" alt="Stenosis 2" width="500">
</p>

・MPI performance

<p>
  <img src="images/mpi_performance.png" alt="MPI Performance" width="400">
</p>