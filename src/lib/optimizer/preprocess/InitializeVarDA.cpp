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

  ue.resize(nz, VDOUBLE2D(ny, VDOUBLE1D(nx, 0e0)));
  ve.resize(nz, VDOUBLE2D(ny, VDOUBLE1D(nx, 0e0)));
  we.resize(nz, VDOUBLE2D(ny, VDOUBLE1D(nx, 0e0)));

  feedbackForce.resize(numOfNodeGlobalFluid, VDOUBLE1D(3,0e0));
  grad_allNode.resize(numOfNodeGlobalFluid, VDOUBLE1D(3,0e0));

  X.resize(control_node.size(), VDOUBLE1D(2,0e0));
  grad.resize(control_node.size(), VDOUBLE1D(2,0e0));

  for(int i=0; i<control_node.size(); i++)
  {
    for(int j=0;j<2;j++)
    {
      X[i][j] = bd_u[control_node[i]][j];
    }
  }

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

  // TOP
  if(controlBoundary == "top")
  {
    for(int i=0; i<nz+1; i++){
      for(int j=0; j<nx+1; j++){
        int elm = i*(nx+1)*(ny+1) + (nx+1)*ny;
        if(phiVOF[elm] < 1e-10) continue;
        VINT1D tmp(4, 0e0);
        tmp[0] = element[elm][2];
        tmp[1] = element[elm][3];
        tmp[2] = element[elm][6];
        tmp[3] = element[elm][7];
        control_elm_node.push_back(tmp);
      }
    }
  }

  // BOTTOM
  else if(controlBoundary == "bottom")
  {
    for(int i=0; i<nz+1; i++){
      for(int j=0; j<nx+1; j++){
        int elm = i*(nx+1)*(ny+1);
        if(phiVOF[elm] < 1e-10) continue;
        VINT1D tmp(4, 0e0);
        tmp[0] = element[elm][0];
        tmp[1] = element[elm][1];
        tmp[2] = element[elm][4];
        tmp[3] = element[elm][5];
        control_elm_node.push_back(tmp);
      }
    }
  }


  // LEFT
  else if(controlBoundary == "left")
  {
    for(int i=0; i<nz+1; i++){
      for(int j=0; j<ny+1; j++){
        int elm = (nx+1)*(ny+1)*i + (nx+1)*j;
        if(phiVOF[elm] < 1e-10) continue;
        VINT1D tmp(4, 0e0);
        tmp[0] = element[elm][0];
        tmp[1] = element[elm][3];
        tmp[2] = element[elm][4];
        tmp[3] = element[elm][7];
        control_elm_node.push_back(tmp);
      }
    }
  }


  // RIGHT
  else if(controlBoundary == "right")
  {
    for(int i=0; i<nz+1; i++){
      for(int j=0; j<ny+1; j++){
        int elm = nx + (nx+1)*(ny+1)*i + (nx+1)*j;
        if(phiVOF[elm] < 1e-10) continue;
        VINT1D tmp(4, 0e0);
        tmp[0] = element[elm][1];
        tmp[1] = element[elm][2];
        tmp[2] = element[elm][5];
        tmp[3] = element[elm][6];
        control_elm_node.push_back(tmp);
      }
    }
  }


  // FRONT
  else if(controlBoundary == "front")
  {
    for(int i=0; i<ny+1; i++){
      for(int j=0; j<nx+1; j++){
        int elm = j + (nx+1)*i;
        if(phiVOF[elm] < 1e-10) continue;
        VINT1D tmp(4, 0e0);
        tmp[0] = element[elm][0];
        tmp[1] = element[elm][1];
        tmp[2] = element[elm][2];
        tmp[3] = element[elm][3];
        control_elm_node.push_back(tmp);
      }
    }
  }
  

  // BACK
  else if(controlBoundary == "back")
  {
    for(int i=0; i<ny+1; i++){
      for(int j=0; j<nx+1; j++){
        int elm = j + (nx+1)*i + (nx+1)*(ny+1)*nz;
        if(phiVOF[elm]<1e-10) continue;
        VINT1D tmp(4, 0e0);
        tmp[0] = element[elm][4];
        tmp[1] = element[elm][5];
        tmp[2] = element[elm][6];
        tmp[3] = element[elm][7];
        control_elm_node.push_back(tmp);
      }
    }
  }


  for(int i=0; i<control_elm_node.size(); i++){
    for(int p=0; p<4; p++){
      control_node.push_back(control_elm_node[i][p]);
    }
  }

  sort(control_node.begin(), control_node.end());
  control_node.erase(unique(control_node.begin(), control_node.end()), control_node.end());

}