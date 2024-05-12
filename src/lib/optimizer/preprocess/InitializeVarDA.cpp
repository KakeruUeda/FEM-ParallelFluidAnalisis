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

  feedbackForceFluid.resize(numOfNodeGlobalFluid, VDOUBLE1D(3,0e0));

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
  cout << "0_0_0" << endl;
  readInputVarDA(); 
  setControlBooundary();
  cout << "0_0_1" << endl;
  
  adjointSolverPetscFluid = (PetscSolver*) new PetscSolver;

  setPetscAdjointSolver();
  
  SolnDataFluid.vecsize = numOfDofsAdjointGlobalFluid;
  SolnDataFluid.solnAdjoint.resize(SolnDataFluid.vecsize);
  SolnDataFluid.solnAdjoint.setZero();
}

void VarDA::setControlBooundary()
{
  bdface_dir.resize(2, 0);
  bdnodeInElm.resize(4, 0);

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
    control_elm_prev_type.resize(numOfElmGlobal, false);
    control_node_prev_type.resize(numOfNodeGlobal, false);
    control_elm_node_prev_type.resize(numOfElmGlobal, VBOOL1D(numOfNodeInElm, false));
    
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
        control_elm_prev_type[elm] = true;
        for(int i=0; i<4; i++){
          control_node_prev_type[tmp[i]] = true;
        }
        control_elm_node_prev_type[elm][0] = true;
        control_elm_node_prev_type[elm][3] = true;
        control_elm_node_prev_type[elm][4] = true;
        control_elm_node_prev_type[elm][7] = true;
      }
    }
    bdface_dir[0] = 1; bdface_dir[1] = 2;
    bdnodeInElm[0] = 0; bdnodeInElm[1] = 2;
    bdnodeInElm[2] = 4; bdnodeInElm[3] = 6;
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
        if(phiVOF[elm] < 1e-10) continue;
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

  if(myId == 0){
    ofstream out_control_node(outputDirTest + "/control_node.dat");
    for(int i=0; i<control_node.size(); i++){
      out_control_node << control_node[i] << endl;
    }
    out_control_node.close();
  }
  MPI_Barrier(MPI_COMM_WORLD);

  control_elm_prev_type_fluid.resize(numOfElmGlobalFluid, false);
  for(int ie=0; ie<numOfElmGlobalFluid; ie++){
    control_elm_prev_type_fluid[ie] = control_elm_prev_type[sortElm[ie]];
  }

  if(myId == 0){
    ofstream out_control_elm_prev_type_fluid(outputDirTest + "/control_elm_prev_type_fluid.dat");
    for(int i=0; i<numOfElmGlobalFluid; i++){
      out_control_elm_prev_type_fluid << control_elm_prev_type_fluid[i] << endl;
    }
    out_control_elm_prev_type_fluid.close();
  }
  MPI_Barrier(MPI_COMM_WORLD);

  control_node_prev_type_fluid.resize(numOfNodeGlobalFluid, false);
  for(int in=0; in<numOfNodeGlobalFluid; in++){
    control_node_prev_type_fluid[in] = control_node_prev_type[sortNode[in]];
  }

  if(myId == 0){
    ofstream out_control_node_prev_type_fluid(outputDirTest + "/control_node_prev_type.dat");
    for(int i=0; i<numOfNodeGlobalFluid; i++){
      out_control_node_prev_type_fluid << control_node_prev_type_fluid[i] << endl;
    }
    out_control_node_prev_type_fluid.close();
  }
  MPI_Barrier(MPI_COMM_WORLD);

  control_elm_node_prev_type_fluid.resize(numOfElmGlobal, VBOOL1D(numOfNodeInElm, false));
  for(int ie=0; ie<numOfElmGlobalFluid; ie++){
    control_elm_node_prev_type_fluid[ie] = control_elm_node_prev_type[sortElm[ie]];
  }

  if(myId == 0){
    ofstream out_control_elm_node_prev_type_fluid(outputDirTest + "/control_elm_node_prev_type_fluid.dat");
    for(int i=0; i<numOfElmGlobalFluid; i++){
      for(int p=0; p<numOfNodeInElm; p++){
        out_control_elm_node_prev_type_fluid << control_elm_node_prev_type_fluid[i][p] << " ";
      }
      out_control_elm_node_prev_type_fluid << endl;
    }
    out_control_elm_node_prev_type_fluid.close();
  }
  MPI_Barrier(MPI_COMM_WORLD);

  control_node_fluid.resize(control_node.size());
  for(int in=0; in<numOfNodeGlobalFluid; in++){
    for(int ib=0; ib<control_node.size(); ib++){
      if(control_node[ib] == in){
        control_node_fluid[in] = control_node[sortNode[in]];
      }
    }
  }

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

  u_adjoint_fluid.resize(numOfNodeGlobalFluid);
  v_adjoint_fluid.resize(numOfNodeGlobalFluid);
  w_adjoint_fluid.resize(numOfNodeGlobalFluid);
  p_adjoint_fluid.resize(numOfNodeGlobalFluid);
  lambda_u_fluid.resize(numOfNodeGlobalFluid);
  lambda_v_fluid.resize(numOfNodeGlobalFluid);
  lambda_w_fluid.resize(numOfNodeGlobalFluid);

  u_adjoint.resize(numOfNodeGlobal);
  v_adjoint.resize(numOfNodeGlobal);
  w_adjoint.resize(numOfNodeGlobal);
  p_adjoint.resize(numOfNodeGlobal);
  lambda_u.resize(numOfNodeGlobal);
  lambda_v.resize(numOfNodeGlobal);
  lambda_w.resize(numOfNodeGlobal);

  return;
}

void VarDA::setPetscAdjointSolver()
{
  cout << "0_0" << endl;
  int jpn, nsize, size, kk;
  numOfDofsNodePrevAdjointFluid.resize(numOfNodeGlobalFluid, 0);
  numOfDofsNodeInElementPrevAdjointFluid.resize(numOfElmGlobalFluid, VINT1D(numOfNodeInElm, 0));

  for(int in=0; in<numOfNodeGlobalFluid; in++){
    if(control_node_prev_type_fluid[in]){
      numOfDofsNodePrevAdjointFluid[in] = numOfDofsNode + 3;
    }else{
      numOfDofsNodePrevAdjointFluid[in] = numOfDofsNode;
    }
  }

  //if(myId == 0){
  //  ofstream out_numOfDofsNodePrevAdjointFluid(outputDirTest + "/numOfDofsNodePrevAdjointFluid.dat");
  //  for(int in=0; in<numOfNodeGlobalFluid; in++){
  //    out_numOfDofsNodePrevAdjointFluid << numOfDofsNodePrevAdjointFluid[in] << endl;
  //  }
  //  out_numOfDofsNodePrevAdjointFluid.close();
  //}
  //MPI_Barrier(MPI_COMM_WORLD);

  for(int ie=0; ie<numOfElmGlobalFluid; ie++){
    for(int p=0; p<numOfNodeInElm; p++){
      if(control_elm_node_prev_type_fluid[ie][p]){
        numOfDofsNodeInElementPrevAdjointFluid[ie][p] = numOfDofsNode + 3;
      }else{
        numOfDofsNodeInElementPrevAdjointFluid[ie][p] = numOfDofsNode;
      }
    }
  }


  //if(myId == 0){
  //  ofstream out_numOfDofsNodeInElementPrevAdjointFluid(outputDirTest + "/numOfDofsNodeInElementPrevAdjointFluid.dat");
  //  for(int ie=0; ie<numOfElmGlobalFluid; ie++){
  //    for(int p=0; p<numOfNodeInElm; p++){
  //      out_numOfDofsNodeInElementPrevAdjointFluid << numOfDofsNodeInElementPrevAdjointFluid[ie][p] << " ";
  //    }
  //    out_numOfDofsNodeInElementPrevAdjointFluid << endl;
  //  }
  //  out_numOfDofsNodeInElementPrevAdjointFluid.close();
  //}
  //MPI_Barrier(MPI_COMM_WORLD);

  int dofsPrevStart = 0;
  nodeDofsMapPrevFluid.resize(numOfNodeGlobalFluid);
  for(int in=0; in<numOfNodeGlobalFluid; in++){
    nodeDofsMapPrevFluid[in] = dofsPrevStart;
    dofsPrevStart += numOfDofsNodePrevAdjointFluid[in];
  }

  //if(myId == 0){
  //  ofstream out_nodeDofsMapPrevFluid(outputDirTest + "/nodeDofsMapPrevFluid.dat");
  //  for(int in=0; in<numOfNodeGlobalFluid; in++){
  //    out_nodeDofsMapPrevFluid << nodeDofsMapPrevFluid[in] << endl;
  //  }
  //  out_nodeDofsMapPrevFluid.close();
  //}
  //MPI_Barrier(MPI_COMM_WORLD);

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
    nodeTypePrevAdjointFluid[DirichletBCsFluidPrev_tmp[i][0]][DirichletBCsFluidPrev_tmp[i][1]] = true;
  }

  //if(myId == 0){
  //  ofstream out_nodeTypePrevAdjointFluid(outputDirTest + "/nodeTypePrevAdjointFluid.dat");
  //  for(int in=0; in<numOfNodeGlobalFluid; in++){
  //    for(int p=0; p<numOfDofsNodePrevAdjointFluid[in]; p++){
  //      out_nodeTypePrevAdjointFluid << nodeTypePrevAdjointFluid[in][p] << " ";
  //    }
  //    out_nodeTypePrevAdjointFluid << endl;
  //  }
  //  out_nodeTypePrevAdjointFluid.close();
  //}
  //MPI_Barrier(MPI_COMM_WORLD);

  numOfDofsAdjointGlobalFluid = 0;
  for(int in=0; in<numOfNodeGlobalFluid; in++){
    for(int j=0; j<numOfDofsNodePrevAdjointFluid[in]; j++){
      nodeDofArrayPrevAdjointFluid[in][j] = numOfDofsAdjointGlobalFluid;
      nodeDofArrayBCsPrevAdjointFluid[in][j] = numOfDofsAdjointGlobalFluid++;
      if(nodeTypePrevAdjointFluid[in][j]){
        nodeDofArrayBCsPrevAdjointFluid[in][j] = -1;
      }
    }
  }

  //  if(myId == 0){
  //    ofstream out_nodeDofArrayBCsPrevAdjointFluid(outputDirTest + "/nodeDofsArrayBCsPrevAdjointFluid.dat");
  //    for(int in=0; in<numOfNodeGlobalFluid; in++){
  //      for(int p=0; p<numOfDofsNodePrevAdjointFluid[in]; p++){
  //        out_nodeDofArrayBCsPrevAdjointFluid << nodeDofArrayBCsPrevAdjointFluid[in][p] << " ";
  //      }
  //      out_nodeDofArrayBCsPrevAdjointFluid << endl;
  //    }
  //    out_nodeDofArrayBCsPrevAdjointFluid.close();
  //  }
  //  MPI_Barrier(MPI_COMM_WORLD);
  //
  //  if(myId == 0){
  //    ofstream out_nodeDofArrayPrevAdjointFluid(outputDirTest + "/nodeDofsArrayPrevAdjointFluid.dat");
  //    for(int in=0; in<numOfNodeGlobalFluid; in++){
  //      for(int p=0; p<numOfDofsNodePrevAdjointFluid[in]; p++){
  //        out_nodeDofArrayPrevAdjointFluid << nodeDofArrayPrevAdjointFluid[in][p] << " ";
  //      }
  //      out_nodeDofArrayPrevAdjointFluid << endl;
  //    }
  //    out_nodeDofArrayPrevAdjointFluid.close();
  //    exit(1);
  //  }
  //  MPI_Barrier(MPI_COMM_WORLD);
  cout << "0_2" << endl;

  
  for(int ie=0; ie<numOfElmGlobalFluid; ie++){
    elmFluid[ie]->dofsNumsPrevFluid.resize(numOfNodeInElm);
    for(int p=0; p<numOfNodeInElm; p++){
      elmFluid[ie]->dofsNumsPrevFluid[p] = nodeDofArrayPrevAdjointFluid[elmFluid[ie]->nodeNumsPrevFluid[p]][0];
    }
  }

  //f(myId == 0){
  // ofstream out_dofsNumsPrevFluid(outputDirTest + "/dofsNumsPrevFluid.dat");
  // for(int ie=0; ie<numOfElmGlobalFluid; ie++){
  //   for(int p=0; p<numOfNodeInElm; p++){
  //     out_dofsNumsPrevFluid << elmFluid[ie]->dofsNumsPrevFluid[p] << " ";
  //   }
  //   out_dofsNumsPrevFluid << endl;
  // }
  // out_dofsNumsPrevFluid.close();
  // exit(1);
  //
  //PI_Barrier(MPI_COMM_WORLD);



  /////////////////
  ////  SORT   ////
  /////////////////
  control_node_type_fluid.resize(numOfNodeGlobalFluid, false);
  control_elm_type_fluid.resize(numOfNodeGlobalFluid, false);
  numOfDofsNodeAdjointFluid.resize(numOfNodeGlobalFluid, 0);
  numOfDofsNodeInElementAdjointFluid.resize(numOfElmGlobalFluid, VINT1D(numOfNodeInElm, 0));

  for(int in=0; in<numOfNodeGlobalFluid; in++){
    int n1 = nodeMapPrevFluid[in];
    control_node_type_fluid[n1] = control_node_prev_type_fluid[in];
  }

  for(int in=0; in<numOfNodeGlobalFluid; in++){
    if(control_node_type_fluid[in]){
      numOfDofsNodeAdjointFluid[in] = numOfDofsNode + 3;
    }else{
      numOfDofsNodeAdjointFluid[in] = numOfDofsNode;
    }
  }
  
  for(int ie=0; ie<numOfElmGlobalFluid; ie++){
    for(int p=0; p<numOfNodeInElm; p++){
      numOfDofsNodeInElementAdjointFluid[ie][p] = numOfDofsNodeAdjointFluid[elmFluid[ie]->nodeNumsPrevFluid[p]];
    }
  }

  for(int ie=0; ie<numOfElmGlobalFluid; ie++){
    bool check = false;
    for(int p=0; p<numOfNodeInElm; p++){
      if(numOfDofsNodeInElementAdjointFluid[ie][p] > numOfDofsNode) check = true;
    }
    if(check) control_elm_type_fluid[ie] = true;
  }

  int dofsStart = 0;
  nodeDofsMapFluid.resize(numOfNodeGlobalFluid);
  for(int in=0; in<numOfNodeGlobalFluid; in++){
    nodeDofsMapFluid[in] = dofsStart;
    dofsStart += numOfDofsNodeAdjointFluid[in];
  }
 //if(myId == 0){
 //  ofstream out_nodeDofsMapFluid(outputDirTest + "/nodeDofsMapFluid.dat");
 //  for(int in=0; in<numOfNodeGlobalFluid; in++){
 //    out_nodeDofsMapFluid << nodeDofsMapFluid[in] << endl;
 //  }
 //  out_nodeDofsMapFluid.close();
 //}
 //MPI_Barrier(MPI_COMM_WORLD);
  


  for(int in=0; in<numOfNodeGlobalFluid; in++){
    VBOOL1D vecBoolTempFalse(numOfDofsNodeAdjointFluid[in], false);
    nodeTypeAdjointFluid.push_back(vecBoolTempFalse);
  }

  for(int in=0; in<numOfNodeGlobalFluid; in++){
    VINT1D vecIntTempM1(numOfDofsNodeAdjointFluid[in], -1);
    nodeDofArrayBCsAdjointFluid.push_back(vecIntTempM1);
    nodeDofArrayAdjointFluid.push_back(vecIntTempM1);
  }

  for(int ib=0; ib<numOfBdNodeFluid; ib++)
  { 
    int n1 = nodeMapFluid[DirichletBCsFluidPrev_tmp[ib][0]];
    DirichletBCsFluid_tmp[ib][0] = n1;
    nodeTypeAdjointFluid[n1][DirichletBCsFluid_tmp[ib][1]] = true;
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

  //if(myId == 0){
  // ofstream out_nodeNumsFluid(outputDirTest + "/nodeNumsFluid.dat");
  // for(int ie=0; ie<numOfElmGlobalFluid; ie++){
  //   for(int p=0; p<numOfNodeInElm; p++){
  //     out_nodeNumsFluid << elmFluid[ie]->nodeNumsFluid[p] << " ";
  //   }
  //   out_nodeNumsFluid << endl;
  // }
  // out_nodeNumsFluid.close();
  //}
  //MPI_Barrier(MPI_COMM_WORLD);

  //if(myId == 0){
  // ofstream out_nodeMapFluid(outputDirTest + "/nodeMapFluid.dat");
  // for(int in=0; in<numOfNodeGlobalFluid; in++){
  //   out_nodeMapFluid << nodeMapFluid[in] << endl;
  // }
  // out_nodeMapFluid.close();
  //}
  //MPI_Barrier(MPI_COMM_WORLD);


  //if(myId == 0){
  // ofstream out_nodeDofsArrayAdjointFluid(outputDirTest + "/nodeDofsArrayAdjointFluid.dat");
  // for(int in=0; in<numOfNodeGlobalFluid; in++){
  //   for(int p=0; p<numOfDofsNodeAdjointFluid[in]; p++){
  //     out_nodeDofsArrayAdjointFluid << nodeDofArrayAdjointFluid[in][p] << " ";
  //   }
  //   out_nodeDofsArrayAdjointFluid << endl;
  // }
  // out_nodeDofsArrayAdjointFluid.close();
  //}
  //MPI_Barrier(MPI_COMM_WORLD);


  for(int ie=0; ie<numOfElmGlobalFluid; ie++){
    elmFluid[ie]->dofsNumsFluid.resize(numOfNodeInElm);
    for(int p=0; p<numOfNodeInElm; p++){
      elmFluid[ie]->dofsNumsFluid[p] = nodeDofArrayAdjointFluid[elmFluid[ie]->nodeNumsPrevFluid[p]][0];
    }
  }

  //if(myId == 0){
  // ofstream out_dofsNumsFluid(outputDirTest + "/dofsNumsFluid.dat");
  // for(int ie=0; ie<numOfElmGlobalFluid; ie++){
  //   for(int p=0; p<numOfNodeInElm; p++){
  //     out_dofsNumsFluid << elmFluid[ie]->dofsNumsFluid[p] << " ";
  //   }
  //   out_dofsNumsFluid << endl;
  // }
  // out_dofsNumsFluid.close();
  //}
  //MPI_Barrier(MPI_COMM_WORLD);


  if(jpn != numOfDofsAdjointGlobalFluid)
  {
    cerr << "(Adjoint) Something wrong with NodeDofArrayNew " << endl;
    cout << "jpn = " << jpn << '\t' << numOfDofsAdjointGlobalFluid << '\t' << myId << endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  
  row_start_adjoint = 0;
  for(int i=0; i<node_start; i++){
    row_start_adjoint += numOfDofsNodeAdjointFluid[i];
  }

  row_end_adjoint = 0;
  for(int i=0; i<=node_end; i++){
    row_end_adjoint += numOfDofsNodeAdjointFluid[i];
  }
  row_end_adjoint = row_end_adjoint - 1;

  numOfDofsAdjointLocalFluid = 0;
  for(int i=node_start; i<=node_end; i++){
    for(int j=0; j<numOfDofsNodeAdjointFluid[i]; j++){
      numOfDofsAdjointLocalFluid++;
    }
  }
  printf(" numOfDofsLocalFluidAdjoint = %5d/%5d \t row_start_adjoint  = %5d \t row_end_adjoint  = %5d \t myId  = %5d \n",numOfDofsAdjointLocalFluid, numOfDofsAdjointGlobalFluid, row_start_adjoint, row_end_adjoint, myId);
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
      if(!control_elm_type_fluid[ie]){
        elmFluid[ie]->prepareAdjointElemData(nsize);
      }else{
        elmFluid[ie]->prepareAdjointElemDataBd(numOfDofsNodeInElementAdjointFluid, nsize, ie);
      }

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
