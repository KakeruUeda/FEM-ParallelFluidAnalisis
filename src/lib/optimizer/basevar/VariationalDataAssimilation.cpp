#include "VariationalDataAssimilation.h"


VarDA::~VarDA()
{
}

void VarDA::setAdjointMatAndVecZero()
{
  int ii, ic, size;
  int Japan, Oita;  
  int Japan2, Oita2; 
  
  Japan = numOfNodeInElm * numOfDofsNode;
  Oita = Japan * Japan;

  Japan2 = numOfNodeInElm * numOfDofsNode + 4 * 3;
  Oita2 = Japan2 * Japan2;
  
  PetscScalar  Flocal_tmp[Japan];  
  PetscScalar  Klocal_tmp[Oita];   

  PetscScalar  FlocalBd_tmp[Japan2];  
  PetscScalar  KlocalBd_tmp[Oita2];   
  
  VINT1D  vecIntTemp;
  
  for(ic=0; ic<numOfElmGlobalFluid; ic++)
  { 
    for(ii=0; ii<Japan; ii++)  Flocal_tmp[ii] = 0e0;
    for(ii=0; ii<Oita; ii++)  Klocal_tmp[ii] = 0e0;
    for(ii=0; ii<Japan2; ii++)  FlocalBd_tmp[ii] = 0e0;
    for(ii=0; ii<Oita2; ii++)  KlocalBd_tmp[ii] = 0e0;
    
    if(elmFluid[ic]->getSubdomainId() == myId)
    {
      size = elmFluid[ic]->nodeForAssyBCsAdjointFluid.size();
      vecIntTemp = elmFluid[ic]->nodeForAssyAdjointFluid;
      
      if(!control_elm_prev_type[ic]){
        MatSetValues(adjointSolverPetscFluid->mtx, size, &vecIntTemp[0], size, &vecIntTemp[0], Klocal_tmp, INSERT_VALUES);
        VecSetValues(adjointSolverPetscFluid->rhsVec, size, &vecIntTemp[0], Flocal_tmp, INSERT_VALUES);
      }else{
        MatSetValues(adjointSolverPetscFluid->mtx, size, &vecIntTemp[0], size, &vecIntTemp[0], KlocalBd_tmp, INSERT_VALUES);
        VecSetValues(adjointSolverPetscFluid->rhsVec, size, &vecIntTemp[0], FlocalBd_tmp, INSERT_VALUES);
      }
    }
  }
}

void VarDA::assignBCsAdjoint()
{
  int ii, n1, n2, n3;
  DirichletBCsAdjointFluid.resize(numOfDofsAdjointGlobalFluid);

  for(ii=0; ii<numOfBdNodeFluid; ii++)
  {
    n1 = int(DirichletBCsFluid_tmp[ii][0]);
    n2 = int(DirichletBCsFluid_tmp[ii][1]);
    n3 = 0;
    for(int jj=0; jj<n1; jj++){
      n3 += numOfDofsNodePrevAdjointFluid[jj];
    }
    n3 += n2;
    DirichletBCsAdjointFluid[n3] = DirichletBCsFluid_tmp[ii][2];
  }

  return;
}

void VarDA::applyBCsAdjoint()
{
  int ii, ic, size;
  int Japan, Oita;  
  int Japan2, Oita2; 
  
  Japan = numOfNodeInElm * numOfDofsNode;
  Oita = Japan * Japan;

  Japan2 = numOfNodeInElm * numOfDofsNode + 4 * 3;
  Oita2 = Japan2 * Japan2;
  
  PetscScalar  Flocal_tmp[Japan];  
  PetscScalar  Klocal_tmp[Oita];   

  PetscScalar  FlocalBd_tmp[Japan2];  
  PetscScalar  KlocalBd_tmp[Oita2];   

  VINT1D  vecIntTemp;

  for(ic=0; ic<numOfElmGlobalFluid; ic++)
  {
    if(elmFluid[ic]->getSubdomainId() == myId)
    {
      for(ii=0; ii<Japan;  ii++)  Flocal_tmp[ii] = 0e0;
      for(ii=0; ii<Oita;   ii++)  Klocal_tmp[ii] = 0e0;
      for(ii=0; ii<Japan2; ii++)  FlocalBd_tmp[ii] = 0e0;
      for(ii=0; ii<Oita2;  ii++)  KlocalBd_tmp[ii] = 0e0;
      
      size = elmFluid[ic]->nodeForAssyBCsAdjointFluid.size();
      vecIntTemp = elmFluid[ic]->nodeForAssyAdjointFluid;
      
      if(!control_elm_prev_type[ic]){
        for(ii=0; ii<size; ii++){
          if(elmFluid[ic]->nodeForAssyBCsAdjointFluid[ii] == -1){
            Klocal_tmp[ii * numOfNodeInElm * numOfDofsNode + ii] = 1;
            Flocal_tmp[ii] = 0e0; // weight function would be zero at boundary(?) //DirichletBCsAdjointFluid[elmFluid[ic]->globalDOFnumsAdjointFluid[ii]];    
          }
        }
        MatSetValues(solverPetscFluid->mtx,    size, &vecIntTemp[0], size, &vecIntTemp[0], Klocal_tmp, INSERT_VALUES);
        VecSetValues(solverPetscFluid->rhsVec, size, &vecIntTemp[0], Flocal_tmp, INSERT_VALUES);
      }else{
       for(ii=0; ii<size; ii++){
          if(elmFluid[ic]->nodeForAssyBCsAdjointFluid[ii] == -1){
            KlocalBd_tmp[ii * (numOfNodeInElm * numOfDofsNode + 4 * 3) + ii] = 1;
            FlocalBd_tmp[ii] = 0e0; // weight function would be zero at boundary(?) //DirichletBCsAdjointFluid[elmFluid[ic]->globalDOFnumsAdjointFluid[ii]];    
          }
        }
        MatSetValues(solverPetscFluid->mtx,    size, &vecIntTemp[0], size, &vecIntTemp[0], KlocalBd_tmp, INSERT_VALUES);
        VecSetValues(solverPetscFluid->rhsVec, size, &vecIntTemp[0], FlocalBd_tmp, INSERT_VALUES);
      }

    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  
  MatAssemblyBegin(solverPetscFluid->mtx, MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(solverPetscFluid->mtx, MAT_FLUSH_ASSEMBLY);

  VecAssemblyBegin(solverPetscFluid->rhsVec);
  VecAssemblyEnd(solverPetscFluid->rhsVec);

  return;
}