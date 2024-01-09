#include "FEM.h"


void FEM::assignBCs()
{
  int ii, n1, n2, n3;

  DirichletBCsFluid.resize(numOfDofsNode*numOfNodeGlobalFluid);
  
  for(ii=0; ii<numOfBdNodeFluid; ii++)
  {
    n1 = int(DirichletBCsFluid_tmp[ii][0]);
    n2 = int(DirichletBCsFluid_tmp[ii][1]);
    n3 = n1 * numOfDofsNode + n2;
    DirichletBCsFluid[n3] = DirichletBCsFluid_tmp[ii][2];
  }

  return;
}


void FEM::assignPulsatileBCs(const double t_itr)
{
  int ii, n1, n2, n3;
  double time_now = t_itr * dt;
  double pulse = sin((2e0 * PI / T) * time_now) + 1e0;
  
  for(ii=0; ii<numOfBdNodeFluid; ii++)
  {
    n1 = int(DirichletBCsFluid_tmp[ii][0]);
    n2 = int(DirichletBCsFluid_tmp[ii][1]);
    n3 = n1 * numOfDofsNode+n2;

    if(DirichletBCsFluid_tmp[ii][2] > 0){
      DirichletBCsFluid[n3] = DirichletBCsFluid_tmp[ii][2] * pulse;
    }else{
      DirichletBCsFluid[n3] = DirichletBCsFluid_tmp[ii][2];
    }
  }
  return;
}


void FEM::applyBCs()
{
  int ii, ic, size;
  int Japan, Oita;  

  Japan = numOfNodeInElm * numOfDofsNode;
  Oita = Japan * Japan;
  
  PetscScalar  Flocal_tmp[Japan];  
  PetscScalar  Klocal_tmp[Oita];  

  VINT1D  vecIntTemp;
  for(ic=0; ic<numOfElmGlobalFluid; ic++)
  {
    if(elmFluid[ic]->getSubdomainId() == myId)
    {
      for(ii=0; ii<Japan; ii++){
        Flocal_tmp[ii] = 0e0;
      }  
      for(ii=0; ii<Oita; ii++){
        Klocal_tmp[ii] = 0e0;
      }
      size = elmFluid[ic]->nodeForAssyBCsFluid.size();
      vecIntTemp = elmFluid[ic]->nodeForAssyFluid;
      
      for(ii=0; ii<size; ii++){
        if(elmFluid[ic]->nodeForAssyBCsFluid[ii] == -1){
          Klocal_tmp[ii * numOfNodeInElm * numOfDofsNode + ii] = 1;
          Flocal_tmp[ii] = DirichletBCsFluid[elmFluid[ic]->globalDOFnumsFluid[ii]];    
        }
      }
      VecSetValues(solverPetscFluid->rhsVec, size, &vecIntTemp[0], Flocal_tmp, ADD_VALUES);
      MatSetValues(solverPetscFluid->mtx,    size, &vecIntTemp[0], size, &vecIntTemp[0], Klocal_tmp, ADD_VALUES);
    }
  }

  return;
}

double FEM::calc_tau(const double (&dxdr)[3][3], const double (&vel)[3])
{
  double tau = 0e0;
  double drdx[3][3];
  BasicFunctions::calcInverseMatrix_3x3(drdx,dxdr);

  double term1 = (2e0 / dt) * (2e0 / dt);
  
  double G[3][3];
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      G[i][j] = 0e0;
      for(int k=0;k<3;k++) G[i][j] += drdx[k][i] * drdx[k][j];
    }
  }
 
  double term2 = 0e0;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      term2 += vel[i] * G[i][j] * vel[j];
    }
  }
  
  double CI = 36e0;
  
  double term3 = 0e0;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      term3 += G[i][j] * G[i][j];
    }
  }
  term3 = CI * term3 / Re;
  
  return tau = pow(term1 + term2 + term3, -5e-1);
}


double FEM::calc_tau2(const double (&vel)[3])
{
  double tau = 0e0;
  double velMag = sqrt(vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);

  double term1 = (2e0 / dt) * (2e0 / dt);

  double max_length = dx;
  
  if(max_length < dy) max_length = dy;
  if(max_length < dz) max_length = dz;
  
  double he = max_length;
  
  double term2 = (2e0 * velMag / he) * (2e0 * velMag / he);
  
  double term3 = (4e0 / (Re * he * he)) * (4e0 / (Re * he * he));
  
  return tau = pow(term1 + term2 + term3, -5e-1);
}


double FEM::calc_tau3(VDOUBLE2D &dNdx, double (&vel)[3])
{
  double tau = 0e0;
  double velMag = sqrt(vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);

  double term1 = (2e0/dt)*(2e0/dt);

  double he = 0e0;

  if(velMag == 0)
  {
    double max_length = dx;
    if(max_length < dy) max_length = dy;
    if(max_length < dz) max_length = dz;
    he = max_length;
  }
  else if(velMag != 0)
  {
    for(int i=0; i<numOfNodeInElm; i++)
    {
      double h_tmp = 0e0;
      for(int k=0; k<3; k++)
      {
        h_tmp += vel[k]/velMag * dNdx[i][k];
      }
      he += fabs(h_tmp);
    }
    he = 2e0/he;
  }


  double term2 = (2e0 * velMag / he) * (2e0 * velMag / he);
  
  double term3 = (4e0 / (Re * he * he)) * (4e0 / (Re * he * he));
  
  return tau = pow(term1+term2+term3,-5e-1);
}