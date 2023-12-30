#include "FEM.h"

void FEM::MatAssySTT(const int ic,MatrixXd &Klocal, VectorXd &Flocal)
{
  int ii, jj;
  int IU,IV,IW,IP;
  int JU,JV,JW,JP;
  
  vector<vector<double>> x_current(numOfNodeInElm,vector<double>(3,0e0));
  
  vector<double> N(numOfNodeInElm,0);
  vector<vector<double>> dNdr(numOfNodeInElm,vector<double>(3,0e0));
  vector<vector<double>> dNdx(numOfNodeInElm,vector<double>(3,0e0));

  vector<vector<double>> K(numOfNodeInElm,vector<double>(numOfNodeInElm,0e0));
  
  for(int i=0;i<numOfNodeInElm;i++){
    for(int j=0;j<3;j++){
      x_current[i][j] = xFluid[elmFluid[ic]->nodeNumsPrevFluid[i]][j];
    }
  }
  double dxdr[3][3];

  Gauss gauss(2);
  int GP = 2;

  //// stabilization parameter ////
  double h = dx/2;
  double tau = h*h/mu/12e0;

  double detJ, weight;
  
  for(int i1=0; i1<GP; i1++){
    for(int i2=0; i2<GP; i2++){
      for(int i3=0; i3<GP; i3++){
        ShapeFunction3D::C3D8_N(N,gauss.point[i1],gauss.point[i2],gauss.point[i3]);
        ShapeFunction3D::C3D8_dNdr(dNdr,gauss.point[i1],gauss.point[i2],gauss.point[i3]);
        
        MathFEM::calc_dxdr(dxdr,dNdr,x_current,numOfNodeInElm);
            
        detJ = dxdr[0][0]*dxdr[1][1]*dxdr[2][2]+dxdr[0][1]*dxdr[1][2]*dxdr[2][0]+dxdr[0][2]*dxdr[1][0]*dxdr[2][1]
                          -dxdr[0][2]*dxdr[1][1]*dxdr[2][0]-dxdr[0][1]*dxdr[1][0]*dxdr[2][2]-dxdr[0][0]*dxdr[1][2]*dxdr[2][1];
        weight = gauss.weight[i1]*gauss.weight[i2]*gauss.weight[i3];
              
        MathFEM::calc_dNdx(dNdx,dNdr,dxdr,numOfNodeInElm);
        
        for(ii=0;ii<numOfNodeInElm;ii++)
        {  
          IU = 4*ii;
          IV = IU+1;
          IW = IU+2;
          IP = IU+3;
          for(jj=0;jj<numOfNodeInElm;jj++)
          {  
            JU = 4*jj;
            JV = JU+1;
            JW = JU+2;
            JP = JU+3;

            K[ii][jj] = 0e0;
            
            for(int k=0;k<3;k++){
              K[ii][jj] += dNdx[ii][k]*dNdx[jj][k];
            }

            //// disfusion ////
            Klocal(IU, JU) -= mu * K[ii][jj] * detJ * weight;
            Klocal(IV, JV) -= mu * K[ii][jj] * detJ * weight;
            Klocal(IW, JW) -= mu * K[ii][jj] * detJ * weight;

            //// pressure ////
            Klocal(IU, JP) += N[jj] * dNdx[ii][0] * detJ * weight;
            Klocal(IV, JP) += N[jj] * dNdx[ii][1] * detJ * weight;
            Klocal(IW, JP) += N[jj] * dNdx[ii][2] * detJ * weight;
        
            //// continuity ////
            Klocal(IP, JU) += N[ii] * dNdx[jj][0] * detJ * weight;
            Klocal(IP, JV) += N[ii] * dNdx[jj][1] * detJ * weight;
            Klocal(IP, JW) += N[ii] * dNdx[jj][2] * detJ * weight;
            
            //// PSPG ////
            Klocal(IP, JP) += tau * K[ii][jj] * detJ * weight;
          }
        }
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
  
  double h = dx/2;
  double tau = h*h/(4e0*mu)/3e0;

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
