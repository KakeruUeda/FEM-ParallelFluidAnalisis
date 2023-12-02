#include "FEM.h"
using namespace std;
 
int main(int argc, char* argv[])
{    
  MPI_Init(NULL, NULL);
  
  string  petscfile   = argv[2];
  PetscInitialize(NULL, NULL, petscfile.c_str(), NULL);
  
  FEM Fem;
  
  //input and check tp file
  string input_file = argv[1];
  int ierror;
  if ((ierror = Fem.tp.read(input_file)) != TP_NO_ERROR) {
   printf("\tError at reading '%s' file\n", input_file.c_str());
    return 1;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  Fem.initialize();
  MPI_Barrier(MPI_COMM_WORLD);
  Fem.Stokes();
  MPI_Barrier(MPI_COMM_WORLD);


  PetscFinalize();
  return 0;
}