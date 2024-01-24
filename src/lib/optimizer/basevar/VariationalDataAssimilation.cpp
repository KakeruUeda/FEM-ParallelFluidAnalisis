#include "VariationalDataAssimilation.h"


VarDA::~VarDA()
{
}

void VarDA::setAdjointMatAndVecZero()
{
  int ii, ic, size;
  int Japan, Oita;  

  Japan = numOfNodeInElm * numOfDofsNode;
  Oita = Japan * Japan;
  
  PetscScalar  Flocal_tmp[Japan];  
  PetscScalar  Klocal_tmp[Oita];   
  VINT1D  vecIntTemp;
  
  for(ic=0; ic<numOfElmGlobalFluid; ic++)
  { 
    for(ii=0; ii<Japan; ii++)  Flocal_tmp[ii] = 0.0;
    for(ii=0; ii<Oita; ii++)  Klocal_tmp[ii] = 0.0;
    if(elmFluid[ic]->getSubdomainId() == myId)
    {
      size = elmFluid[ic]->nodeForAssyBCsFluid.size();
      vecIntTemp = elmFluid[ic]->nodeForAssyFluid;
      MatSetValues(adjointSolverPetscFluid->mtx, size, &vecIntTemp[0], size, &vecIntTemp[0], Klocal_tmp, INSERT_VALUES);
      VecSetValues(adjointSolverPetscFluid->rhsVec, size, &vecIntTemp[0], Flocal_tmp, INSERT_VALUES);
    }
  }
}