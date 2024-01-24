#include "VariationalDataAssimilation.h"

void VarDA::initializeVarDA()
{
  mainInitialize();
  adjointInitialize();

  string str, base_label, label;
  string mapFile;
  
  base_label = "/VarDA";
  label = base_label + "/velocityMAP";

  if ( !tp.getInspectedValue(label, mapFile)){
    cout << label<< " is not set" << endl;
    exit(0);
  }
  obs.inputVelocityMAP(mapFile);

  ue.resize(obs.nz, VDOUBLE2D(obs.ny, VDOUBLE1D(obs.nx, 0e0)));
  ve.resize(obs.nz, VDOUBLE2D(obs.ny, VDOUBLE1D(obs.nx, 0e0)));
  we.resize(obs.nz, VDOUBLE2D(obs.ny, VDOUBLE1D(obs.nx, 0e0)));

  feedbackForce.resize(numOfNodeGlobal, VDOUBLE1D(3,0e0));
  grad_allNode.resize(numOfNodeGlobal, VDOUBLE1D(3,0e0));

  X.resize(control_node.size(), VDOUBLE1D(2,0e0));
  grad.resize(control_node.size(), VDOUBLE1D(2,0e0));

  for(int i=0; i<control_node.size(); i++)
  {
    for(int j=0;j<2;j++)
    {
      X[i][j] = bd_u[control_node[i]][j];
    }
  }

  resizeDAVariables();

  return;
}

void VarDA::mainInitialize()
{
  FEM::initialize(); 
}

void VarDA::adjointInitialize()
{
  readInputVarDA(); 
  setControlBooundary();
  
  if(adjointSolverPetscFluid != NULL) delete adjointSolverPetscFluid;
  adjointSolverPetscFluid = NULL;
  adjointSolverPetscFluid = (PetscSolver*) new PetscSolver;
  
  setPetscAdjointSolver();
}

void VarDA::setControlBooundary()
{
  bdface_dir.resize(2, 0e0);

  // TOP
  if(controlBoundary == "top")
  {
    for(int k=0; k<nz; k++){
      for(int i=0; i<nx; i++){
        int elm = k*nx*ny + i + (nx-1)*(ny-1);
        if(phiVOF[elm] < 1e-10) continue;
        VINT1D tmp(4, 0e0);
        tmp[0] = element[elm][2];
        tmp[1] = element[elm][3];
        tmp[2] = element[elm][6];
        tmp[3] = element[elm][7];
        control_elm_node.push_back(tmp);
      }
    }
    bdface_dir[0] = 0; bdface_dir[1] = 2;
  }

  // BOTTOM
  else if(controlBoundary == "bottom")
  {
    for(int k=0; k<nz; k++){
      for(int i=0; i<nx; i++){
        int elm = k*nx*ny + i;
        if(phiVOF[elm] < 1e-10) continue;
        VINT1D tmp(4, 0e0);
        tmp[0] = element[elm][0];
        tmp[1] = element[elm][1];
        tmp[2] = element[elm][4];
        tmp[3] = element[elm][5];
        control_elm_node.push_back(tmp);
      }
    }
    bdface_dir[0] = 0; bdface_dir[1] = 2;
  }

  // LEFT
  else if(controlBoundary == "left")
  {
    for(int k=0; k<nz; k++){
      for(int j=0; j<ny; j++){
        int elm = nx*ny*k + nx*j;
        if(phiVOF[elm] < 1e-10) continue;
        VINT1D tmp(4, 0e0);
        tmp[0] = element[elm][0];
        tmp[1] = element[elm][3];
        tmp[2] = element[elm][4];
        tmp[3] = element[elm][7];
        control_elm_node.push_back(tmp);
      }
    }
    bdface_dir[0] = 1; bdface_dir[1] = 2;
  }

  // RIGHT
  else if(controlBoundary == "right")
  {
    for(int k=0; k<nz; k++){
      for(int j=0; j<ny; j++){
        int elm = nx-1 + nx*ny*k + nx*j;
        if(phiVOF[elm] < 1e-10) continue;
        VINT1D tmp(4, 0e0);
        tmp[0] = element[elm][1];
        tmp[1] = element[elm][2];
        tmp[2] = element[elm][5];
        tmp[3] = element[elm][6];
        control_elm_node.push_back(tmp);
      }
    }
    bdface_dir[0] = 0; bdface_dir[1] = 2;
  }

  // FRONT
  else if(controlBoundary == "front")
  {
    for(int j=0; j<ny; j++){
      for(int i=0; i<nx; i++){
        int elm = i + ny*i;
        if(phiVOF[elm] < 1e-10) continue;
        VINT1D tmp(4, 0e0);
        tmp[0] = element[elm][0];
        tmp[1] = element[elm][1];
        tmp[2] = element[elm][2];
        tmp[3] = element[elm][3];
        control_elm_node.push_back(tmp);
      }
    }
    bdface_dir[0] = 0; bdface_dir[1] = 1;
  }

  // BACK
  else if(controlBoundary == "back")
  {
    for(int j=0; j<ny; j++){
      for(int i=0; i<nx; i++){
        int elm = i + nx*j + nx*ny*(nz-1);
        if(phiVOF[elm]<1e-10) continue;
        VINT1D tmp(4, 0e0);
        tmp[0] = element[elm][4];
        tmp[1] = element[elm][5];
        tmp[2] = element[elm][6];
        tmp[3] = element[elm][7];
        control_elm_node.push_back(tmp);
      }
    }
    bdface_dir[0] = 0; bdface_dir[1] = 1;
  }

  for(int i=0; i<control_elm_node.size(); i++){
    for(int p=0; p<4; p++){
      control_node.push_back(control_elm_node[i][p]);
    }
  }

  sort(control_node.begin(), control_node.end());
  control_node.erase(unique(control_node.begin(), control_node.end()), control_node.end());

  return;
}

void VarDA::resizeDAVariables()
{
  cfd.u.resize(obs.nz, VDOUBLE2D(obs.ny, VDOUBLE1D(obs.nx, 0e0)));
  cfd.v.resize(obs.nz, VDOUBLE2D(obs.ny, VDOUBLE1D(obs.nx, 0e0)));
  cfd.w.resize(obs.nz, VDOUBLE2D(obs.ny, VDOUBLE1D(obs.nx, 0e0)));

  ue.resize(obs.nz, VDOUBLE2D(obs.ny, VDOUBLE1D(obs.nx, 0e0)));
  ve.resize(obs.nz, VDOUBLE2D(obs.ny, VDOUBLE1D(obs.nx, 0e0)));
  we.resize(obs.nz, VDOUBLE2D(obs.ny, VDOUBLE1D(obs.nx, 0e0)));

  ue_edge.resize(obs.nz+2, VDOUBLE2D(obs.ny+2, VDOUBLE1D(obs.nx+2, 0e0)));
  ve_edge.resize(obs.nz+2, VDOUBLE2D(obs.ny+2, VDOUBLE1D(obs.nx+2, 0e0)));
  we_edge.resize(obs.nz+2, VDOUBLE2D(obs.ny+2, VDOUBLE1D(obs.nx+2, 0e0)));
  
  return;
}

void VarDA::setPetscAdjointSolver()
{
  int jpn, nsize,size, kk;
  numOfDofsNodePrevAdjointFluid.resize(numOfNodeGlobalFluid, 0);
  numOfDofsNodeInElementPrevAdjointFluid.resize(numOfElmGlobalFluid, VINT1D(numOfNodeInElm, 0));

  for(int in=0; in<numOfNodeGlobalFluid; in++){
    numOfDofsNodePrevAdjointFluid[in] = numOfDofsNode;
    for(int ib=0; ib<control_node.size(); ib++){
      if(control_node[ib] == in){
        numOfDofsNodePrevAdjointFluid[in] = numOfDofsNode + 3;
      }
    }
  }

  for(int ie=0; ie<numOfElmGlobalFluid; ie++){
    for(int p=0; p<numOfNodeInElm; p++){
      numOfDofsNodeInElementPrevAdjointFluid[ie][p] = numOfDofsNode;
    }
    for(int ib=0; ib<control_elm_node.size(); ib++){
      for(int p=0; p<4; p++){
        if(control_elm_node[ib][p] == elementFluidPrev[ie][p]){
          numOfDofsNodeInElementPrevAdjointFluid[ie][p] = numOfDofsNode + 3;
        }
      }
    }
  }


  for(int in=0; in<numOfNodeGlobalFluid; in++){
    VBOOL1D vecBoolTempFalse(numOfDofsNodePrevAdjointFluid[in], false);
    nodeTypePrevAdjointFluid.push_back(vecBoolTempFalse);
  }
 
  for(int in=0; in<numOfNodeGlobalFluid; in++){
    VINT1D vecIntTempM1(numOfDofsNodePrevAdjointFluid[in], -1);
    nodeDofArrayBCsPrevAdjointFluid.push_back(vecIntTempM1);
    nodeDofArrayPrevAdjointFluid.push_back(vecIntTempM1);
  }

  for(int i=0; i<numOfBdNodeFluid; i++){
    nodeTypePrevAdjointFluid[DirichletBCsFluid_tmp[i][0]][DirichletBCsFluid_tmp[i][1]] = true;
  }

  for(int in=0; in<numOfNodeGlobalFluid; in++){
    for(int j=0; j<numOfDofsNodePrevAdjointFluid[in]; j++){
      nodeDofArrayPrevAdjointFluid[in][j] = numOfDofsAdjointGlobalFluid;
      nodeDofArrayBCsPrevAdjointFluid[in][j] = numOfDofsAdjointGlobalFluid++;
      if(nodeTypePrevAdjointFluid[in][j]){
        nodeDofArrayBCsPrevAdjointFluid[in][j] = -1;
      }
    }
  }

  /////////////////
  ////  SORT   ////
  /////////////////
  numOfDofsNodeAdjointFluid.resize(numOfNodeGlobalFluid, 0);
  numOfDofsNodeInElementAdjointFluid.resize(numOfElmGlobalFluid, VINT1D(numOfNodeInElm, 0));

  for(int in=0; in<numOfNodeGlobalFluid; in++){
    numOfDofsNodeAdjointFluid[in] = numOfDofsNode;
    for(int ib=0; ib<control_node.size(); ib++){
      if(control_node[ib] == nodeMapFluid[in]){
        numOfDofsNodeAdjointFluid[in] = numOfDofsNode + 3;
      }
    }
  }

  for(int ie=0; ie<numOfElmGlobalFluid; ie++){
    for(int p=0; p<numOfNodeInElm; p++){
      numOfDofsNodeInElementAdjointFluid[ie][p] = numOfDofsNode;
    }
    for(int ib=0; ib<control_elm_node.size(); ib++){
      for(int p=0; p<4; p++){
        if(control_elm_node[ib][p] == elementFluid[ie][p]){
          numOfDofsNodeInElementAdjointFluid[ie][p] = numOfDofsNode + 3;
        }
      }
    }
  }

  for(int in=0; in<numOfNodeGlobalFluid; in++){
    VBOOL1D vecBoolTempFalse(numOfDofsNodeAdjointFluid[in], false);
    nodeTypeAdjointFluid.push_back(vecBoolTempFalse);
  }

  for(int in=0; in<numOfNodeGlobalFluid; in++){
    VINT1D vecIntTempM1(numOfDofsNodeAdjointFluid[in], -1);
    nodeDofArrayBCsAdjointFluid.push_back(vecIntTempM1);
    nodeDofArrayAdjointFluid.push_back(vecIntTempM1);
  }

  for(int i=0; i<numOfBdNodeFluid; i++){ 
    int n1 = nodeMapFluid[DirichletBCsFluid_tmp[i][0]];
    DirichletBCsFluid_tmp[i][0] = n1;
    nodeTypeAdjointFluid[n1][DirichletBCsFluid_tmp[i][1]] = true;
  }

  jpn = 0;
  for(int in=0; in<numOfNodeGlobalFluid; in++){
    for(int j=0; j<numOfDofsNodeAdjointFluid[in]; j++){
      nodeDofArrayAdjointFluid[in][j] = jpn;
      nodeDofArrayBCsAdjointFluid[in][j] = jpn++;
      if(nodeTypeAdjointFluid[in][j]){
        nodeDofArrayBCsAdjointFluid[in][j] = -1;
      }
    }
  }

  if(jpn != numOfDofsAdjointGlobalFluid)
  {
    cerr << "(Adjpoint) Something wrong with NodeDofArrayNew " << endl;
    cout << "jpn = " << jpn << '\t' << numOfDofsAdjointGlobalFluid << '\t' << myId << endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  for(int i=0; i<node_start; i++){
    row_start_adjoint += numOfDofsNodeAdjointFluid[i];
  }
  for(int i=0; i<node_end+1; i++){
    row_end_adjoint += numOfDofsNodeAdjointFluid[i];
  }
  for(int i=node_start; i<=node_end; i++){
    for(int j=0; j<numOfDofsNodeAdjointFluid[i]; j++){
      numOfDofsAdjointLocalFluid++;
    }
  }
  printf(" numOfDofsLocalFluidAdjoint = %5d/%5d \t row_start_adjoint  = %5d \t row_end_adjoint  = %5d \t myId  = %5d \n", numOfDofsAdjointLocalFluid, numOfDofsAdjointGlobalFluid, row_start_adjoint, row_end_adjoint, myId);
  MPI_Barrier(MPI_COMM_WORLD);
  
  jpn=0;
  errpetsc = MPI_Allreduce(&numOfDofsAdjointLocalFluid, &jpn, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
 
  if(jpn != numOfDofsAdjointGlobalFluid)
  {
    cerr << " (Adjoint) Sum of local problem sizes is not equal to global size" << endl;
    cout << " jpn = " << jpn << '\t' << numOfDofsAdjointGlobalFluid << '\t' << myId << endl;
  }
  /////////////////
  //// SORTEND ////
  /////////////////
  MPI_Barrier(MPI_COMM_WORLD);

  for(int ie=0; ie<numOfElmGlobalFluid; ie++){
    if(elmFluid[ie]->getSubdomainId() == myId){
      nsize = 0;
      for(int p=0; p<numOfNodeInElm; p++){
        nsize += numOfDofsNodeInElementAdjointFluid[ie][p];
      }
      //elmFluid[ie]->prepareElemData(nsize);
      //numOfNodeInElm = elmFluid[ie]->nodeNumsFluid.size();

      elmFluid[ie]->nodeForAssyBCsAdjointFluid.resize(nsize);
      elmFluid[ie]->nodeForAssyAdjointFluid.resize(nsize);
      for(int p=0; p<numOfNodeInElm; p++){
        jpn = numOfDofsNodeInElementAdjointFluid[ie][p] * p;
        kk  = elmFluid[ie]->nodeNumsFluid[p];
        for(int j=0; j<numOfDofsNodeInElementAdjointFluid[ie][p]; j++){
          elmFluid[ie]->nodeForAssyBCsAdjointFluid[jpn+j] = nodeDofArrayBCsAdjointFluid[kk][j];
          elmFluid[ie]->nodeForAssyAdjointFluid[jpn+j] =  nodeDofArrayAdjointFluid[kk][j];
        }
      }
    }
  }   
  MPI_Barrier(MPI_COMM_WORLD);

  assyForSolnAdjointFluid.resize(numOfDofsAdjointGlobalFluid);
  jpn = 0;
  for(int in=0; in<numOfNodeGlobalFluid; in++){
    for(int j=0; j<numOfDofsNodeAdjointFluid[in]; j++){
      assyForSolnAdjointFluid[jpn++] = in * numOfDofsNodeAdjointFluid[in] + j;
    }
  }

  vector<set<int>> forAssyMatAdjointFluid;
  set<int>::iterator it;
  int  *tt; int r;
  nsize = 0;

  forAssyMatAdjointFluid.resize(numOfDofsAdjointGlobalFluid);

  for(int ie=0; ie<numOfElmGlobalFluid; ie++)
  {
    if(elmFluid[ie]->getSubdomainId() == myId)
    {
      tt = &(elmFluid[ie]->nodeForAssyAdjointFluid[0]);
      nsize = elmFluid[ie]->nodeForAssyAdjointFluid.size();

      for(int ii=0;ii<nsize; ii++)
      {
        r = tt[ii];

        if(r != -1)
        {
          if(r >= row_start_adjoint && r <= row_end_adjoint)
          {
          for(int jj=0;jj<nsize;jj++)
            {
              if(tt[jj] != -1)
              {
                forAssyMatAdjointFluid[r].insert(tt[jj]);
              }
            }
          }
        }
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  adjointSolverPetscFluid->nnz_max_row = 0;
  for(int i=row_start_adjoint; i<=row_end_adjoint; i++){
    size = forAssyMatAdjointFluid[i].size();
    adjointSolverPetscFluid->nnz_max_row = max(adjointSolverPetscFluid->nnz_max_row, size);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  adjointSolverPetscFluid->initialize(numOfDofsAdjointLocalFluid, numOfDofsAdjointGlobalFluid);

  return;
}