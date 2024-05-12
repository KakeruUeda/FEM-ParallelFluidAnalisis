#include "FEM.h"

void FEM::initialize()
{ 
  readInput(); 
  setDomain();
  setBoundary();
  setFluidDomain();
  resizeVariables();

  allocateObj();
  prepareMatrix();

  //octreeSubDivision();
  //interfacePartition();
}


void FEM::setDomain()
{
  x.resize(numOfNodeGlobal,VDOUBLE1D(3, 0));
  element.resize(numOfElmGlobal);

  for(int ic=0;ic<numOfElmGlobal;ic++){
    element[ic].resize(8);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  int tmp2=0;
  for(int k=0; k<nz+1; k++){
    for(int j=0; j<ny+1; j++){
      for(int i=0; i<nx+1; i++){
        x[tmp2][0]=i*dx;
        x[tmp2][1]=j*dy;
        x[tmp2][2]=k*dz; 
        tmp2++;
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  tmp2=0;
  for(int k=0; k<nz; k++){
    for(int j=0; j<ny; j++){
      for(int i=0; i<nx; i++){
        element[tmp2][0]= i   + j*(nx+1)     + k*(nx+1)*(ny+1);
        element[tmp2][1]= i+1 + j*(nx+1)     + k*(nx+1)*(ny+1);
        element[tmp2][2]= i+1 + (j+1)*(nx+1) + k*(nx+1)*(ny+1);
        element[tmp2][3]= i   + (j+1)*(nx+1) + k*(nx+1)*(ny+1);
        element[tmp2][4]= i   + j*(nx+1)     + (k+1)*(nx+1)*(ny+1);
        element[tmp2][5]= i+1 + j*(nx+1)     + (k+1)*(nx+1)*(ny+1);
        element[tmp2][6]= i+1 + (j+1)*(nx+1) + (k+1)*(nx+1)*(ny+1);
        element[tmp2][7]= i   + (j+1)*(nx+1) + (k+1)*(nx+1)*(ny+1);
        tmp2++;
      }
    }
  }

  return;
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
  VDOUBLE1D vecDbTmp(3,0);

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
  
  VINT1D sortNodeNew(numOfNodeGlobal,0); 

  // sortNodeNew: start with 0
  numOfNodeGlobalFluid = 0;
  for(auto it = sortNode.begin(); it != sortNode.end(); it++){
    sortNodeNew[*it] = numOfNodeGlobalFluid;
    numOfNodeGlobalFluid++;
  }

  int n1, n2;
  elementFluid.resize(numOfElmGlobalFluid,VINT1D(numOfNodeInElm,0));
  for(int ii=0; ii<numOfElmGlobalFluid; ii++){
    for(int jj=0; jj<numOfNodeInElm; jj++){
      n1 = elementFluidPrev[ii][jj];
      elementFluid[ii][jj] = sortNodeNew[n1];
    }
  }

  bd_iu_fluid.resize(numOfNodeGlobalFluid,VINT1D(3,0));
  bd_u_fluid.resize(numOfNodeGlobalFluid,VDOUBLE1D(3,0e0));
  bd_ip_fluid.resize(numOfNodeGlobalFluid,0);
  bd_p_fluid.resize(numOfNodeGlobalFluid,0e0);

  phiFluid.resize(numOfElmGlobalFluid,0e0);
  phiEXFluid.resize(numOfElmGlobalFluid,0e0);
  phiVOFFluid.resize(numOfElmGlobalFluid,0e0);
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
    phiVOFFluid[ii] = phiVOF[sortElm[ii]];
  }

  xFluid.resize(numOfNodeGlobalFluid,VDOUBLE1D(3, 0));
  for(int ii=0; ii<numOfNodeGlobalFluid; ii++){
    for(int kk=0; kk<3; kk++){
      xFluid[ii][kk] = x[sortNode[ii]][kk];
    }
  }

  int count = 0;
  VDOUBLE1D vecDbTmp(3,0);

  for(int ii=0;ii<numOfNodeGlobalFluid;ii++){
    for(int kk=0; kk<3; kk++){
      if(bd_iu_fluid[ii][kk] == 0){
        vecDbTmp[0] = ii;
        vecDbTmp[1] = kk;
        vecDbTmp[2] = bd_u_fluid[ii][kk];
        DirichletBCsFluidPrev_tmp.push_back(vecDbTmp);
        count++;
      }
    }
  }
  for(int ii=0;ii<numOfNodeGlobalFluid;ii++){
    if(bd_ip_fluid[ii] == 0){
      vecDbTmp[0] = ii;
      vecDbTmp[1] = 3;
      vecDbTmp[2] = bd_p_fluid[ii];
      DirichletBCsFluidPrev_tmp.push_back(vecDbTmp);
      count++;
    }
  }


  numOfBdNodeFluid = count;

  //VINT2D().swap(bd_iu);
  //VDOUBLE2D().swap(bd_u);
  //VINT1D().swap(bd_ip);
  //VDOUBLE1D().swap(bd_p);
  //VDOUBLE2D().swap(x);
  
}


void FEM::visualizeDomain()
{
  if(myId > 0) return;
  string vtiFile;

  vtiFile = outputDirMain + "/meshPartition.vti";
  export_file.export_vti_metis(vtiFile, nodeId, elmId, nx, ny, nz, dx, dy, dz);
    
  if(bd == BOUNDARY::XFEM)
  {
    vtiFile = outputDirMain + "/domain.vti";
    export_file.export_vti_domain(vtiFile, sdf, phi, phiEX, nx, ny, nz, dx, dy, dz);
  }
  else if(bd == BOUNDARY::DARCY)
  {    
    vtiFile = outputDirMain + "/phiVOF.vti";
    export_file.export_vti_elm(vtiFile, phiVOF, nx, ny, nz, dx, dy, dz);
  }

  return;
}


