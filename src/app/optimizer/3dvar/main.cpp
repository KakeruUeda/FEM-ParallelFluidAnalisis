#include "VariationalDataAssimilation.h"

int main(int argc,char *argv[])
{ 
  MPI_Init(NULL, NULL);
  string  petscfile   = argv[2];
  PetscInitialize(NULL, NULL, petscfile.c_str(), NULL);
  
  VarDA var;
  var.solver = SOLVER::STEADY_NAVIERSTOKES;
  
  //input and check tp file
  string input_file = argv[1];
  int ierror;
  if ((ierror = var.tp.read(input_file)) != TP_NO_ERROR) {
   printf("\t Error at reading '%s' file\n", input_file.c_str());
    return 1;
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  var.initializeVarDA();

  var.nu = var.mu / var.rho;
  var.Re = 1e0 / var.nu;

  for(int opt_itr=0; opt_itr<var.optMaxItr; opt_itr++)
  {
    if(opt_itr != 0){
      MPI_Barrier(MPI_COMM_WORLD); 
      var.SteadyNavierStokes();
    }

    MPI_Barrier(MPI_COMM_WORLD); 
    var.calcCostFunction();
  }

 
  MPI_Barrier(MPI_COMM_WORLD);
  PetscFinalize(); MPI_Finalize();

  return 0;

}