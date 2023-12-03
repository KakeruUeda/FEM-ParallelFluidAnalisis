#include "FEM.h"
using namespace std;

void FEM::initialize()
{
  readInput(); 
  setDomain();
  setBoundary();
  //setImage();
  prepareMatrix();
}


void FEM::setDomain()
{

  x.resize(numOfNodeGlobal,vector<double>(3, 0));
  element.resize(numOfElmGlobal);

  for(int ic=0;ic<numOfElmGlobal;ic++){
    element[ic].resize(8);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  int tmp2=0;
  for(int k=0; k<nz+1; k++){
    for(int i=0; i<ny+1; i++){
      for(int j=0; j<nx+1; j++){
        x[tmp2][0]=j*dx;
        x[tmp2][1]=i*dy;
        x[tmp2][2]=k*dz; 
        tmp2++;
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  tmp2=0;
  for(int k=0;k<nz;k++){
    for(int i=0;i<ny;i++){
      for(int j=0;j<nx;j++){
          element[tmp2][0]=j   +i*(nx+1) +k*(nx+1)*(ny+1);
          element[tmp2][1]=j+1 +i*(nx+1) +k*(nx+1)*(ny+1);
          element[tmp2][2]=j+1 +(i+1)*(nx+1) +k*(nx+1)*(ny+1);
          element[tmp2][3]=j   +(i+1)*(nx+1) +k*(nx+1)*(ny+1);
          element[tmp2][4]=j   +i*(nx+1) +(k+1)*(nx+1)*(ny+1);
          element[tmp2][5]=j+1 +i*(nx+1) +(k+1)*(nx+1)*(ny+1);
          element[tmp2][6]=j+1 +(i+1)*(nx+1) +(k+1)*(nx+1)*(ny+1);
          element[tmp2][7]=j   +(i+1)*(nx+1) +(k+1)*(nx+1)*(ny+1);
          tmp2++;
      }
    }
  }
  prepare();
}

void FEM::prepare()
{
  elm = new ElementBaseFEM* [numOfElmGlobal];

  for(int ee=0;ee<numOfElmGlobal;++ee){
    elm[ee] = new ElementBaseFEM();
    elm[ee]->SolnData = &(SolnData);
  }

  if(solverPetsc != NULL) delete solverPetsc;
  solverPetsc = NULL;

  solverPetsc = (PetscSolver*) new PetscSolver;

  SolnData.initialise(numOfNodeGlobal*numOfDofsNode);
}



void FEM::setBoundary(){

  int count = 0;
  vector<double> vecDbTmp(3,0);

  for(int i=0;i<numOfNodeGlobal;i++){
    for(int k=0; k<3; k++){
      if(bd_iu[i][k] == 0){
        vecDbTmp[0] = i;
        vecDbTmp[1] = k;
        vecDbTmp[2] = bd_u[i][k];
        DirichletBCs_tmp.push_back(vecDbTmp);
        count++;
      }
    }
  }
  
  ///pressure///
  vecDbTmp[0] = 0;
  vecDbTmp[1] = 3;
  vecDbTmp[2] = 1e0;
  DirichletBCs_tmp.push_back(vecDbTmp);
  count++;
  numOfBdNode = count;

}

