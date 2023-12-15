#include "FEM.h"
#include <omp.h>
int main(int argc, char* argv[]){

  MPI_Init(NULL, NULL);
  string  petscfile   = argv[2];
  PetscInitialize(NULL, NULL, petscfile.c_str(), NULL);
 
  FEM STOKES;
  STOKES.solver = SOLVER::STEADY_STOKES;
  
  //input and check tp file
  string input_file = argv[1];
  int ierror;
  if ((ierror = STOKES.tp.read(input_file)) != TP_NO_ERROR) {
   printf("\t Error at reading '%s' file\n", input_file.c_str());
    return 1;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  STOKES.initialize();
  MPI_Barrier(MPI_COMM_WORLD);
  
  if(STOKES.myId == 1){
    string vtiFile;
    vtiFile = STOKES.outputDir + "/domain.vti";
    STOKES.export_vti_domain(vtiFile);
  }

  omp_set_num_threads(STOKES.numOfOMP);

  STOKES.Stokes();
  MPI_Barrier(MPI_COMM_WORLD);

  STOKES.postCaluculation();
  MPI_Barrier(MPI_COMM_WORLD);
    
  STOKES.deallocate();
  MPI_Barrier(MPI_COMM_WORLD);

  PetscFinalize();
  MPI_Finalize();
  
  return 0;
  
}