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