#include "FEM.h"


void FEM::assignBoundaryConditions(){

  int ii, ic, n1, n2, n3, size;
  int Japan, Oita;
  DirichletBCsFluid.resize(numOfDofsNode*numOfNodeGlobalFluid);
  for(ii=0; ii<numOfBdNodeFluid; ii++){
    n1 = int(DirichletBCsFluid_tmp[ii][0]);
    n2 = int(DirichletBCsFluid_tmp[ii][1]);
    n3 = n1*numOfDofsNode+n2;
    DirichletBCsFluid[n3] = DirichletBCsFluid_tmp[ii][2];
  }
  

  Japan = numOfNodeInElm*numOfDofsNode;
  Oita = Japan*Japan;
  
  PetscScalar  Flocal_tmp[Japan];  
  PetscScalar  Klocal_tmp[Oita];  

  vector<int>  vecIntTemp;
  for(ic=0; ic<numOfElmGlobalFluid; ic++)
  {  
    if(elmFluid[ic]->getSubdomainId() == myId)
    {
      for(ii=0; ii<Japan; ii++){
        Flocal_tmp[ii] = 0.0;
      }  
      for(ii=0; ii<Oita; ii++){
        Klocal_tmp[ii] = 0.0;
      }  
      
      size = elmFluid[ic]->nodeForAssyBCsFluid.size();
      vecIntTemp = elmFluid[ic]->nodeForAssyFluid;
      
      for(ii=0; ii<size; ii++){
        if(elmFluid[ic]->nodeForAssyBCsFluid[ii] == -1){
          Klocal_tmp[ii*numOfNodeInElm*numOfDofsNode+ii] = 1;
          Flocal_tmp[ii] = DirichletBCsFluid[elmFluid[ic]->globalDOFnumsFluid[ii]];    
        }
      }
      VecSetValues(solverPetscFluid->rhsVec, size, &vecIntTemp[0], Flocal_tmp, INSERT_VALUES);
      MatSetValues(solverPetscFluid->mtx,    size, &vecIntTemp[0], size, &vecIntTemp[0], Klocal_tmp, INSERT_VALUES);
    }
  }
  errpetsc = MPI_Barrier(MPI_COMM_WORLD);
 
  MatAssemblyBegin(solverPetscFluid->mtx,MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(solverPetscFluid->mtx, MAT_FLUSH_ASSEMBLY);
  VecAssemblyBegin(solverPetscFluid->rhsVec);
  VecAssemblyEnd(solverPetscFluid->rhsVec);

  return;
}
