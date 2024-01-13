#include "FEM.h"

int main(int argc, char* argv[]){

  MPI_Init(NULL, NULL);
  string  petscfile   = argv[2];
  PetscInitialize(NULL, NULL, petscfile.c_str(), NULL);
  
  FEM NAVIER;
  NAVIER.solver = SOLVER::UNSTEADY_NAVIERSTOKES;
  
  //input and check tp file
  string input_file = argv[1];
  int ierror;
  if ((ierror = NAVIER.tp.read(input_file)) != TP_NO_ERROR) {
   printf("\t Error at reading '%s' file\n", input_file.c_str());
    return 1;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  NAVIER.initialize();
  MPI_Barrier(MPI_COMM_WORLD);

  NAVIER.nu = NAVIER.mu / NAVIER.rho;
  NAVIER.Re = 1e0 / NAVIER.nu;

  omp_set_num_threads(NAVIER.numOfOMP);
  MPI_Barrier(MPI_COMM_WORLD);

  NAVIER.UnsteadyNavierStokes();
  MPI_Barrier(MPI_COMM_WORLD);

  //NAVIER.postCaluculation();
  //MPI_Barrier(MPI_COMM_WORLD);
    
  NAVIER.solverDeallocate();
  MPI_Barrier(MPI_COMM_WORLD);

  PetscFinalize(); MPI_Finalize();

  return 0;
}