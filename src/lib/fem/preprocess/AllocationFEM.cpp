#include "FEM.h"

void FEM::allocateObj()
{
  elmFluid = new ElementBaseFEM* [numOfElmGlobalFluid];

  for(int ee=0; ee<numOfElmGlobalFluid; ee++)
  {
    elmFluid[ee] = new ElementBaseFEM();
  }

  if(solverPetscFluid != NULL) delete solverPetscFluid;
  solverPetscFluid = NULL;

  solverPetscFluid = (PetscSolver*) new PetscSolver;

  SolnDataFluid.initialize(numOfNodeGlobalFluid * numOfDofsNode);
}


void FEM::resizeVariables()
{
  if(solver == SOLVER::STEADY_STOKES)
  {
    uFluid.resize(numOfNodeGlobalFluid, 0e0); u.resize(numOfNodeGlobal, 0e0);
    vFluid.resize(numOfNodeGlobalFluid, 0e0); v.resize(numOfNodeGlobal, 0e0);
    wFluid.resize(numOfNodeGlobalFluid, 0e0); w.resize(numOfNodeGlobal, 0e0);
    pFluid.resize(numOfNodeGlobalFluid, 0e0); p.resize(numOfNodeGlobal, 0e0);
  }
  else if(solver == SOLVER::STEADY_NAVIERSTOKES)
  {
    uFluid.resize(numOfNodeGlobalFluid, 0e0); u.resize(numOfNodeGlobal, 0e0);
    vFluid.resize(numOfNodeGlobalFluid, 0e0); v.resize(numOfNodeGlobal, 0e0);
    wFluid.resize(numOfNodeGlobalFluid, 0e0); w.resize(numOfNodeGlobal, 0e0);
    pFluid.resize(numOfNodeGlobalFluid, 0e0); p.resize(numOfNodeGlobal, 0e0);  
  }
  else if(solver == SOLVER::UNSTEADY_NAVIERSTOKES)
  {
    uf.resize(2, VDOUBLE1D(numOfNodeGlobalFluid, 0e0)); u.resize(numOfNodeGlobal, 0e0);
    vf.resize(2, VDOUBLE1D(numOfNodeGlobalFluid, 0e0)); v.resize(numOfNodeGlobal, 0e0);
    wf.resize(2, VDOUBLE1D(numOfNodeGlobalFluid, 0e0)); w.resize(numOfNodeGlobal, 0e0);
    pf.resize(2, VDOUBLE1D(numOfNodeGlobalFluid, 0e0)); p.resize(numOfNodeGlobal, 0e0);
  }
  else{
    cout << "resizeVariables solver is not set" << endl;
    exit(0);
  }  
}