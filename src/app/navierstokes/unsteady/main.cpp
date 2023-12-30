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
  NAVIER.Re = (NAVIER.U * NAVIER.D) / NAVIER.nu;
  
  if(NAVIER.myId == 0){
    string vtiFile;
    if(NAVIER.bd == BOUNDARY::XFEM)
    {
      vtiFile = NAVIER.outputDir + "/domain.vti";
      NAVIER.export_vti_domain(vtiFile);
    }
    else if(NAVIER.bd == BOUNDARY::DARCY)
    {
      vtiFile = NAVIER.outputDir + "/sdf_node.vti";
      NAVIER.export_vti_node(vtiFile,NAVIER.sdf_node);
      vtiFile = NAVIER.outputDir + "/phiVOF.vti";
      NAVIER.export_vti_elm(vtiFile,NAVIER.phiVOF);
    }
  }

  omp_set_num_threads(NAVIER.numOfOMP);

  NAVIER.UnsteadyNavierStokes();
  MPI_Barrier(MPI_COMM_WORLD);

  NAVIER.postCaluculation();
  MPI_Barrier(MPI_COMM_WORLD);
    
  NAVIER.deallocate();
  MPI_Barrier(MPI_COMM_WORLD);

  PetscFinalize();
  MPI_Finalize();

  return 0;
}