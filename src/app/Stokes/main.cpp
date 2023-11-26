#include "FEM.h"
using namespace std;
 
int main(int argc, char* argv[])
{
  PetscInitialize(NULL, NULL,(char *)0, NULL);

  FEM Fem;

  //input argument
  if(argc!=2){
    printf("Invalid input. Please set tp file\n");
    return 1;
  }
  
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