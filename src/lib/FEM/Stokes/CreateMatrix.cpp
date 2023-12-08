#include "FEM.h"
using namespace std;

void FEM::calcStokesMatrix(const int ic,MatrixXd &Klocal, VectorXd &Flocal)
{
  vector<vector<double>> x_current(numOfNodeInElm,vector<double>(3,0e0));
  vector<double> N(numOfNodeInElm,0);
  vector<vector<double>> dNdr(numOfNodeInElm,vector<double>(3,0e0));

  for(int i=0;i<numOfNodeInElm;i++){
    for(int j=0;j<3;j++){
      x_current[i][j] = x[elm[ic]->nodeNumsPrev[i]][j];
    }
  }
  //#pragma omp parallel for
  Gauss g1(1),g2(2),g3(3);
  for(int i1=0; i1<2; i1++){
    for(int i2=0; i2<2; i2++){
      for(int i3=0; i3<2; i3++){
        ShapeFunction3D::C3D8_N(N,g1.point[i1],g1.point[i2],g1.point[i3]);
        ShapeFunction3D::C3D8_dNdr(dNdr,g1.point[i1],g1.point[i2],g1.point[i3]);
        DiffusionInGaussIntegral(Klocal,Flocal,dNdr,x_current,numOfNodeInElm,g1.weight[i1]*g1.weight[i2]*g1.weight[i3],ic);
        PressureInGaussIntegral(Klocal,Flocal,N,dNdr,x_current,numOfNodeInElm,g1.weight[i1]*g1.weight[i2]*g1.weight[i3],ic);
        PSPGInGaussIntegral(Klocal,Flocal,dNdr,x_current,numOfNodeInElm,g1.weight[i1]*g1.weight[i2]*g1.weight[i3],ic);
      }
    }
  }
}


void FEM::DiffusionInGaussIntegral(MatrixXd &Klocal, VectorXd &Flocal,vector<vector<double>> &dNdr,vector<vector<double>> &x_current,const int numOfNodeInElm,const double weight,const int ic)
{
  int ii, jj;
  int TI,TIp1,TIp2,TIp3;
  int TJ,TJp1,TJp2,TJp3;

  vector<vector<double>> dNdx(numOfNodeInElm,vector<double>(3,0));

  double dxdr[3][3];

  MathFEM::calc_dxdr(dxdr,dNdr,x_current,numOfNodeInElm);
  double detJ=dxdr[0][0]*dxdr[1][1]*dxdr[2][2]+dxdr[0][1]*dxdr[1][2]*dxdr[2][0]+dxdr[0][2]*dxdr[1][0]*dxdr[2][1]
             -dxdr[0][2]*dxdr[1][1]*dxdr[2][0]-dxdr[0][1]*dxdr[1][0]*dxdr[2][2]-dxdr[0][0]*dxdr[1][2]*dxdr[2][1];
  MathFEM::calc_dNdx(dNdx,dNdr,dxdr,numOfNodeInElm);

  vector<vector<double>> K(numOfNodeInElm,vector<double>(numOfNodeInElm,0));


  for(ii=0;ii<numOfNodeInElm;ii++)
  {
    TI   = 4*ii;
    TIp1 = TI+1;
    TIp2 = TI+2;
    TIp3 = TI+3;
    for(jj=0;jj<numOfNodeInElm;jj++)
    {
      TJ   = 4*jj;
      TJp1 = TJ+1;
      TJp2 = TJ+2;
      TJp3 = TJ+3;
     
     K[ii][jj] = 0e0;

      for(int k=0;k<3;k++)
      {
        K[ii][jj] += dNdx[ii][k]*dNdx[jj][k];
      }
      Klocal(TI,   TJ)   -= mu * K[ii][jj] * detJ * weight;
      Klocal(TIp1, TJp1) -= mu * K[ii][jj] * detJ * weight;
      Klocal(TIp2, TJp2) -= mu * K[ii][jj] * detJ * weight;
    }
  }
}



void FEM::PressureInGaussIntegral(MatrixXd &Klocal, VectorXd &Flocal,vector<double> &N,vector<vector<double>> &dNdr,vector<vector<double>> &x_current,const int numOfNodeInElm,const double weight,const int ic)
{ 
  int ii, jj;
  int TI,TIp1,TIp2,TIp3;
  int TJ,TJp1,TJp2,TJp3;

  vector<vector<double>> dNdx(numOfNodeInElm,vector<double>(3,0));

  double dxdr[3][3];
  MathFEM::calc_dxdr(dxdr,dNdr,x_current,numOfNodeInElm);
  double detJ=dxdr[0][0]*dxdr[1][1]*dxdr[2][2]+dxdr[0][1]*dxdr[1][2]*dxdr[2][0]+dxdr[0][2]*dxdr[1][0]*dxdr[2][1]
             -dxdr[0][2]*dxdr[1][1]*dxdr[2][0]-dxdr[0][1]*dxdr[1][0]*dxdr[2][2]-dxdr[0][0]*dxdr[1][2]*dxdr[2][1];
  
  MathFEM::calc_dNdx(dNdx,dNdr,dxdr,numOfNodeInElm);

  for(ii=0;ii<numOfNodeInElm;ii++){
    TI   = 4*ii;
    TIp1 = TI+1;
    TIp2 = TI+2;
    TIp3 = TI+3;
    for(jj=0;jj<numOfNodeInElm;jj++){
      TJ   = 4*jj;
      TJp1 = TJ+1;
      TJp2 = TJ+2;
      TJp3 = TJ+3;
      // pressure term
      Klocal(TI,   TJp3) += N[jj] * dNdx[ii][0] * detJ * weight;
      Klocal(TIp1, TJp3) += N[jj] * dNdx[ii][1] * detJ * weight;
      Klocal(TIp2, TJp3) += N[jj] * dNdx[ii][2] * detJ * weight;

      // continuity equation
      Klocal(TIp3, TJ)   += N[ii] * dNdx[jj][0] * detJ * weight;
      Klocal(TIp3, TJp1) += N[ii] * dNdx[jj][1] * detJ * weight;
      Klocal(TIp3, TJp2) += N[ii] * dNdx[jj][2] * detJ * weight;
    }
  }
}



void FEM::PSPGInGaussIntegral(MatrixXd &Klocal, VectorXd &Flocal,vector<vector<double>> &dNdr,vector<vector<double>> &x_current,const int numOfNodeInElm,const double weight,const int ic)
{  
  int ii, jj;
  int TI,TIp1,TIp2,TIp3;
  int TJ,TJp1,TJp2,TJp3;

  vector<vector<double>> dNdx(numOfNodeInElm,vector<double>(3,0));

  double dxdr[3][3];
  MathFEM::calc_dxdr(dxdr,dNdr,x_current,numOfNodeInElm);
  double detJ=dxdr[0][0]*dxdr[1][1]*dxdr[2][2]+dxdr[0][1]*dxdr[1][2]*dxdr[2][0]+dxdr[0][2]*dxdr[1][0]*dxdr[2][1]
             -dxdr[0][2]*dxdr[1][1]*dxdr[2][0]-dxdr[0][1]*dxdr[1][0]*dxdr[2][2]-dxdr[0][0]*dxdr[1][2]*dxdr[2][1];
  
  MathFEM::calc_dNdx(dNdx,dNdr,dxdr,numOfNodeInElm);

  vector<vector<double>> K(numOfNodeInElm,vector<double>(numOfNodeInElm,0));
  
  //ouble h, h_tmp, h_tmp_tmp;
  //or(int i=0; i<numOfNodeInElm; i++){
  // h_tmp_tmp = 0e0;
  // for(int k=0; k<3; k++){
  //   h_tmp_tmp += dNdx[i][k];
  // }
  // h_tmp += fabs(h_tmp_tmp);
  //
  //=2e0/h_tmp;
  //out << h << endl;
  //out << dx << endl;
  //exit(1);
  //double h = sqrt(dx*dx+dy*dy+dz*dz);
  double h = dx/2;
  double tau = h*h/mu/12e0;

  //double h = sqrt(dx*dx+dy*dy+dz*dz);
  //double tau = 9e0*pow(4e0*mu/(h*h),2e0);
  //tau = 1e0/sqrt(tau);

  //double h = sqrt(dx*dx+dy*dy+dz*dz);
  //double tau = h*h/(4e0*mu)/3e0;
    
  //double h = sqrt(dx*dx+dy*dy+dz*dz);
  //double tau = pow((4e0*mu)/(h*h),2e0);
  //tau = 1e0/sqrt(tau);


  for(ii=0;ii<numOfNodeInElm;ii++){
    TI   = 4*ii;
    TIp1 = TI+1;
    TIp2 = TI+2;
    TIp3 = TI+3;
    for(jj=0;jj<numOfNodeInElm;jj++){
      TJ   = 4*jj;
      TJp1 = TJ+1;
      TJp2 = TJ+2;
      TJp3 = TJ+3;
      K[ii][jj] = 0e0;
      for(int k=0;k<3;k++){
        K[ii][jj] += dNdx[ii][k]*dNdx[jj][k];
      }
      Klocal(TIp3, TJp3) += tau * K[ii][jj] * detJ * weight;
    }
  }
}
