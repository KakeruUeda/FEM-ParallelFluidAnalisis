#include "FEM.h"

FEM::FEM()
{
  MPI_Comm_size(MPI_COMM_WORLD, &numOfId);
  MPI_Comm_rank(MPI_COMM_WORLD, &myId);

  //numOfDofsLocal = numOfDofsGlobal = 0;
  /*
  elm = nullptr;
  solverPetsc = nullptr;
  */

  /// FLUID ONLY ///
  numOfDofsLocalFluid = numOfDofsGlobalFluid = 0;
  
  elmFluid = nullptr;
  solverPetscFluid = nullptr;
}


FEM::~FEM()
{
  /*
  if(elm != nullptr)
  {
    for(unsigned int ii=0;ii<numOfElmGlobal;++ii) delete elm[ii];

    delete [] elm;
    elm = nullptr;
  }
  */
}

int FEM::deallocate()
{
  /*
  if(solverPetsc != nullptr){
    solverPetsc->free();
  }
  */

  return 0;
}
