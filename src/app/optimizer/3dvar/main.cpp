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

  MPI_Barrier(MPI_COMM_WORLD);
  var.visualizeDomain();

  omp_set_num_threads(var.numOfOMP);

  var.nu = var.mu / var.rho;
  var.Re = 1e0 / var.nu;
  
  string costFunction = var.outputDirDA + "/costFunction.dat";
  FILE *fp;
  string vtiFile;
  
  MPI_Barrier(MPI_COMM_WORLD);
  for(int opt_itr=0; opt_itr<var.optMaxItr; opt_itr++)
  {
    if(opt_itr != 0){
      MPI_Barrier(MPI_COMM_WORLD); 
      var.SteadyNavierStokes();
    }

    MPI_Barrier(MPI_COMM_WORLD); 
    var.calcCostFunction();

    MPI_Barrier(MPI_COMM_WORLD); 
    var.calcFeedbackForce();

    VDOUBLE4D fvel(var.nz+1, VDOUBLE3D(var.ny+1, VDOUBLE2D(var.nx+1, VDOUBLE1D(3, 0e0))));

    int tmp2 = 0;
    for(int k=0; k<var.nz+1; k++){
      for(int j=0; j<var.ny+1; j++){
        for(int i=0; i<var.nx+1; i++){
          fvel[k][j][i][0] = var.feedbackForce[tmp2][0];
          fvel[k][j][i][1] = var.feedbackForce[tmp2][1];
          fvel[k][j][i][2] = var.feedbackForce[tmp2][2];
          tmp2++;
        }
      }
    }

    int num = opt_itr/var.output_itr;
    vtiFile  = var.outputDirDA + "/force_" + to_string(num)+".vti";
    var.export_file.export_vti_node_3dir(vtiFile, var.feedbackForce, var.nx, var.ny, var.nz, var.dx, var.dy, var.dz);

    vtiFile  = var.outputDirDA + "/force_2D_" + to_string(num)+".vti";
    var.export_file.export_vti_node_2dir_2D(vtiFile, var.feedbackForce, var.nx, var.ny, var.nz, var.dx, var.dy, var.dz);
    
    //var.export_file.export_vti_velocity_node(vtiFile, fvel, var.obs.nx, var.obs.ny, var.obs.nz, var.obs.dx, var.obs.dy, var.obs.dz);
    vtiFile  = var.outputDirDA + "/error_" + to_string(num)+".vti";
    var.export_file.export_vti_elm_xyz_3dir(vtiFile, var.ue, var.ve, var.we, var.obs.nx, var.obs.ny, var.obs.nz, var.obs.dx, var.obs.dy, var.obs.dz);
    PetscPrintf(MPI_COMM_WORLD, "itr=%d costFunction = %e\n", opt_itr, var.costFunction.total);
    if((fp = fopen(costFunction.c_str(), "a")) == NULL){
      printf("file open error\n");
      exit(1);
    }
    
    MPI_Barrier(MPI_COMM_WORLD); 
    exit(1);
  }


  MPI_Barrier(MPI_COMM_WORLD);
  PetscFinalize(); MPI_Finalize();

  return 0;

}