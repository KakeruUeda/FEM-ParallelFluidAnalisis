#include "FEM.h"

void FEM::initialize()
{ 
  readInput(); 
  setDomain();
  setBoundary();
  setFluidDomain();
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
  /*
  elm = new ElementBaseFEM* [numOfElmGlobal];

  for(int ee=0;ee<numOfElmGlobal;++ee){
    elm[ee] = new ElementBaseFEM();
    elm[ee]->SolnData = &(SolnData);
  }

  if(solverPetsc != NULL) delete solverPetsc;
  solverPetsc = NULL;

  solverPetsc = (PetscSolver*) new PetscSolver;

  SolnData.initialise(numOfNodeGlobal*numOfDofsNode);
  */

}


void FEM::setBoundary(){
  
  for(int ic=0;ic<numOfElmGlobal;ic++){
    
    if(phi[ic]>1e-10) continue;

    for(int i=0;i<numOfNodeInElm;i++){
      for(int j=0;j<3;j++){
        bd_iu[element[ic][i]][j] = 0;
        bd_u[element[ic][i]][j] = 0;
        //bd_ip[element[ic][i]] = 0;
        //bd_p[element[ic][i]] = 0;
      }
    }
  }
  
  /*
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
      if(bd_ip[i] == 0){
        vecDbTmp[0] = i;
        vecDbTmp[1] = 3;
        vecDbTmp[2] = bd_p[i];
        DirichletBCs_tmp.push_back(vecDbTmp);
        count++;
      }
    }
  }

  ///pressure///

  numOfBdNode = count;
  */

}


void FEM::setFluidDomain()
{
  numOfElmGlobalFluid=0;
  
  for(int ic=0; ic<numOfElmGlobal; ic++){
    if(phi[ic]<1e-12){
      continue;
    }
    sortElm.push_back(ic);
    elementFluidPrev.push_back(element[ic]);
    numOfElmGlobalFluid++;
  }


  sortNode.resize(numOfElmGlobalFluid*numOfNodeInElm);
  
  int tmp2 = 0;
  for(int ic=0; ic<numOfElmGlobalFluid; ic++){
    for(int ii=0; ii<numOfNodeInElm; ii++){
      sortNode[tmp2] = elementFluidPrev[ic][ii];
      tmp2++;
    }
  }
  sort(sortNode.begin(), sortNode.end());
  sortNode.erase(unique(sortNode.begin(), sortNode.end()), sortNode.end());
  
  vector<int> sortNodeNew(numOfNodeGlobal,0); 
  

  ///// sortNodeNew: start with 0 /////
  numOfNodeGlobalFluid = 0;
  for(auto it = sortNode.begin(); it != sortNode.end(); it++){
    sortNodeNew[*it] = numOfNodeGlobalFluid;
    numOfNodeGlobalFluid++;
  }

  int n1, n2;
  elementFluid.resize(numOfElmGlobalFluid,vector<int>(numOfNodeInElm,0));
  for(int ii=0; ii<numOfElmGlobalFluid; ii++){
    for(int jj=0; jj<numOfNodeInElm; jj++){
      n1 = elementFluidPrev[ii][jj];
      elementFluid[ii][jj] = sortNodeNew[n1];
    }
  }

  /*
  if(myId == 0){
    for(int i=0; i<numOfElmFluidGlobal; i++){
      for(int j=0; j<numOfNodeInElm; j++){
        //cout << "i = " << i << " j= " << j << " elementFluidPrev node = " << elementFluidPrev[i][j] << endl;
      }
    }
    for(int i=0; i<numOfElmFluidGlobal; i++){
      for(int j=0; j<numOfNodeInElm; j++){
        //cout << "i = " << i << " j= " << j << " elementFluid node = " << elementFluid[i][j] << endl;
      }
    }
    //exit(1);
  }
  */

  bd_iu_fluid.resize(numOfNodeGlobalFluid,vector<int>(3,0));
  bd_u_fluid.resize(numOfNodeGlobalFluid,vector<double>(3,0e0));
  bd_ip_fluid.resize(numOfNodeGlobalFluid,0);
  bd_p_fluid.resize(numOfNodeGlobalFluid,0e0);

  phiFluid.resize(numOfElmGlobalFluid,0e0);
  phiEXFluid.resize(numOfElmGlobalFluid,0e0);
  sdfFluid.resize(numOfNodeGlobalFluid,0e0);
  
  for(int ii=0; ii<numOfNodeGlobalFluid; ii++){
    for(int kk=0; kk<3; kk++){
      bd_iu_fluid[ii][kk] = bd_iu[sortNode[ii]][kk];
      bd_u_fluid[ii][kk] = bd_u[sortNode[ii]][kk];
    }
    sdfFluid[ii] = sdf[sortNode[ii]];
    bd_ip_fluid[ii] = bd_ip[sortNode[ii]];
    bd_p_fluid[ii] = bd_p[sortNode[ii]];
  }

  for(int ii=0; ii<numOfElmGlobalFluid; ii++){
    phiFluid[ii] = phi[sortElm[ii]];
    phiEXFluid[ii] = phiEX[sortElm[ii]];
  }

  xFluid.resize(numOfNodeGlobalFluid,vector<double>(3, 0));
  for(int ii=0; ii<numOfNodeGlobalFluid; ii++){
    for(int kk=0; kk<3; kk++){
      xFluid[ii][kk] = x[sortNode[ii]][kk];
    }
  }

  /*
  if(myId == 0){
    for(int i=0; i<numOfNodeGlobalFluid; i++){
      cout << "i = " << i << " " ;
      for(int j=0; j<3; j++){
        cout << bd_iu_fluid[i][j] << " " ;
      }
      cout << endl;
    }
    for(int i=0; i<numOfNodeGlobalFluid; i++){
      cout << "i = " << i << " " ;
      for(int j=0; j<3; j++){
        cout << bd_u_fluid[i][j] << " " ;
      }
      cout << endl;
    }
    //exit(1);
  }
  */
  
  

  int count = 0;
  vector<double> vecDbTmp(3,0);

  for(int ii=0;ii<numOfNodeGlobalFluid;ii++){
    for(int kk=0; kk<3; kk++){
      if(bd_iu_fluid[ii][kk] == 0){
        vecDbTmp[0] = ii;
        vecDbTmp[1] = kk;
        vecDbTmp[2] = bd_u_fluid[ii][kk];
        DirichletBCsFluid_tmp.push_back(vecDbTmp);
        count++;
      }
    }
  }
  for(int ii=0;ii<numOfNodeGlobalFluid;ii++){
    if(bd_ip_fluid[ii] == 0){
      vecDbTmp[0] = ii;
      vecDbTmp[1] = 3;
      vecDbTmp[2] = bd_p_fluid[ii];
      DirichletBCsFluid_tmp.push_back(vecDbTmp);
      count++;
    }
  }

  /*
  vecDbTmp[0] = 0;
  vecDbTmp[1] = 3;
  vecDbTmp[2] = 0e0;
  DirichletBCsFluid_tmp.push_back(vecDbTmp);
  count++;
  */
  numOfBdNodeFluid = count;


  //// FLUID ONLY PREPARE ////
  elmFluid = new ElementBaseFEM* [numOfElmGlobalFluid];

  for(int ee=0;ee<numOfElmGlobalFluid;ee++){
    elmFluid[ee] = new ElementBaseFEM();
    elmFluid[ee]->SolnDataFluid = &(SolnDataFluid);
  }

  if(solverPetscFluid != NULL) delete solverPetscFluid;
  solverPetscFluid = NULL;

  solverPetscFluid = (PetscSolver*) new PetscSolver;

  SolnDataFluid.initialise(numOfNodeGlobalFluid*numOfDofsNode);

}


