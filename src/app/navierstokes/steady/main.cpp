#include "FEM.h"
#include <omp.h>
int main(int argc, char* argv[]){

  MPI_Init(NULL, NULL);
  string  petscfile   = argv[2];
  PetscInitialize(NULL, NULL, petscfile.c_str(), NULL);
  
  FEM NAVIER;
  NAVIER.solver = SOLVER::STEADY_NAVIERSTOKES;
  
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

  NAVIER.nu = NAVIER.mu/NAVIER.rho;
  NAVIER.Re = (NAVIER.U*NAVIER.D)/NAVIER.nu;
  
  if(NAVIER.myId == 1){
    string vtiFile;
    vtiFile = NAVIER.outputDir + "/domain.vti";
    NAVIER.export_vti_domain(vtiFile);
  }

  omp_set_num_threads(NAVIER.numOfOMP);

  NAVIER.SteadyNavierStokes();
  MPI_Barrier(MPI_COMM_WORLD);

  NAVIER.postCaluculation();
  MPI_Barrier(MPI_COMM_WORLD);
    
  NAVIER.deallocate();
  MPI_Barrier(MPI_COMM_WORLD);

  PetscFinalize();
  MPI_Finalize();

  return 0;
}