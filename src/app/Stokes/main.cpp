#include "FEM.h"
#include <omp.h>
int main(int argc, char* argv[]){

  MPI_Init(NULL, NULL);
  string  petscfile   = argv[2];
  PetscInitialize(NULL, NULL, petscfile.c_str(), NULL);
  FEM fem;
  
  //input and check tp file
  string input_file = argv[1];
  int ierror;
  if ((ierror = fem.tp.read(input_file)) != TP_NO_ERROR) {
   printf("\t Error at reading '%s' file\n", input_file.c_str());
    return 1;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  fem.initialize();
  MPI_Barrier(MPI_COMM_WORLD);
  
  if(fem.myId == 1){
    string vtiFile;
    vtiFile = "domain.vti";
    fem.export_vti_domain(vtiFile);
  }

  omp_set_num_threads(fem.numOfOMP);

  fem.stokes();
  MPI_Barrier(MPI_COMM_WORLD);

  fem.postCaluculation();
  MPI_Barrier(MPI_COMM_WORLD);
    
  fem.deallocate();
  MPI_Barrier(MPI_COMM_WORLD);

  PetscFinalize();
  MPI_Finalize();
  
  return 0;
  
}