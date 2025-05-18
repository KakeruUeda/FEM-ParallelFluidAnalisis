#include "FEM.h"
#include "PostProcess.h"

int main(int argc, char* argv[])
{
  MPI_Init(NULL, NULL);
  string petscfile = argv[2];
  PetscInitialize(NULL, NULL, petscfile.c_str(), NULL);
 
  FEM stokes;
  PostProcess post;
  stokes.solver = SOLVER::STEADY_STOKES;
  
  string input_file = argv[1];
  int ierror;
  if ((ierror = stokes.tp.read(input_file)) != TP_NO_ERROR) {
   printf("\t Error at reading '%s' file\n", input_file.c_str());
    return 1;
  }
  string input_file_post = argv[1];
  int ierror_post;
  if ((ierror_post = post.tp_post.read(input_file_post)) != TP_NO_ERROR) {
   printf("\t Error at reading '%s' file\n", input_file_post.c_str());
    return 1;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  stokes.initialize();

  MPI_Barrier(MPI_COMM_WORLD);  
  post.readDAParam();

  int tmp = 0;
  MPI_Barrier(MPI_COMM_WORLD);  
  for(auto it = post.numOfObsVoxels.begin(); it != post.numOfObsVoxels.end(); it++){
    post.validateDADomain(stokes, post.numOfObsVoxels[tmp][0], 
    post.numOfObsVoxels[tmp][1], post.numOfObsVoxels[tmp][2]);
    tmp++;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  stokes.visualizeDomain();
  
  MPI_Barrier(MPI_COMM_WORLD);
  stokes.Stokes();
  
  MPI_Barrier(MPI_COMM_WORLD);
  stokes.physicalVariables();
  
  tmp = 0;
  MPI_Barrier(MPI_COMM_WORLD);  
  for(auto it = post.numOfObsVoxels.begin(); it != post.numOfObsVoxels.end(); it++){
    post.prepareForDataAssimilation(stokes, post.numOfObsVoxels[tmp][0], 
    post.numOfObsVoxels[tmp][1], post.numOfObsVoxels[tmp][2]);
    tmp++;
  }
  
  MPI_Barrier(MPI_COMM_WORLD);  
  post.extractDomain(stokes);

  MPI_Barrier(MPI_COMM_WORLD); 
  stokes.solverDeallocate();
  
  MPI_Barrier(MPI_COMM_WORLD);
  PetscFinalize(); MPI_Finalize();
  
  return 0;
  
}