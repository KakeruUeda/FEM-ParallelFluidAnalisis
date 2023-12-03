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
  Gauss g1(1),g2(2),g3(3);
  for(int i1=0; i1<2; i1++){
    for(int i2=0; i2<2; i2++){
      for(int i3=0; i3<2; i3++){
        ShapeFunction3D::C3D8_N(N,g1.point[i1],g1.point[i2],g1.point[i3]);
        ShapeFunction3D::C3D8_dNdr(dNdr,g1.point[i1],g1.point[i2],g1.point[i3]);
        Pressure_inGaussIntegral(Klocal,Flocal,N,dNdr,x_current,numOfNodeInElm,g1.weight[i1]*g1.weight[i2]*g1.weight[i3],ic);
        Diffusion_inGaussIntegral(Klocal,Flocal,dNdr,x_current,numOfNodeInElm,g1.weight[i1]*g1.weight[i2]*g1.weight[i3],ic);
        PSPG_inGaussIntegral(Klocal,Flocal,dNdr,x_current,numOfNodeInElm,g1.weight[i1]*g1.weight[i2]*g1.weight[i3],ic);
      }
    }
  }
  /*
  ofstream outKlocal("KlocalDiffusion.dat");
  ofstream outX("x.dat");
  if(this_mpi_proc == 0){
    outKlocal << Klocal << endl;
    for(int i=0;i<numOfNodeInElm;i++){
      outX << "i = " << i << " x = " <<  x_current[i][0] <<  " y = " <<  x_current[i][1] <<  " z = " <<  x_current[i][2] << endl;
    }
    outX << endl;
    outKlocal << endl;
  }
  outKlocal.close();
  outX.close();
  */
  //exit(1);
}


void FEM::Diffusion_inGaussIntegral(MatrixXd &Klocal, VectorXd &Flocal,vector<vector<double>> &dNdr,vector<vector<double>> &x_current,const int numOfNodeInElm,const double weight,const int ic)
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
  /*
  for(int p=0;p<numOfNodeInElm;p++){
    for(int q=0;q<numOfNodeInElm;q++){
      K(p,q) = 0e0;
      for(int j=0;j<3;j++){
        K(p,q) += dNdx(p,j)*dNdx(q,j);
      }
    }
  }
  */

  /*
  for(int i=0;i<numOfNodeInElm;i++){
    for(int j=0;j<numOfNodeInElm;j++){
      LHS_uu[ic](i,j,0,0) += mu * K(i,j) * detJ * weight;
      LHS_uu[ic](i,j,1,1) += mu * K(i,j) * detJ * weight;
      LHS_uu[ic](i,j,2,2) += mu * K(i,j) * detJ * weight;
    }
  }
  */

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
      double fact;
      for(int k=0;k<3;k++){
        K[ii][jj] += dNdx[ii][k]*dNdx[jj][k];
        //Klocal(TI,   TJ)   -= mu * dNdx[ii][k]*dNdx[jj][k] * detJ * weight;
        //Klocal(TIp1, TJp1) -= mu * dNdx[ii][k]*dNdx[jj][k] * detJ * weight;
        //Klocal(TIp2, TJp2) -= mu * dNdx[ii][k]*dNdx[jj][k] * detJ * weight;
        //cout << dNdx[ii][k] << endl;
        //fact += dNdx[ii][k]*dNdx[jj][k];
      }
      //cout << fact << endl;
      //fact = 0e0;

      Klocal(TI,   TJ)   -= mu * K[ii][jj] * detJ * weight;
      Klocal(TIp1, TJp1) -= mu * K[ii][jj] * detJ * weight;
      Klocal(TIp2, TJp2) -= mu * K[ii][jj] * detJ * weight;
    }
  }
}



void FEM::Pressure_inGaussIntegral(MatrixXd &Klocal, VectorXd &Flocal,vector<double> &N,vector<vector<double>> &dNdr,vector<vector<double>> &x_current,const int numOfNodeInElm,const double weight,const int ic)
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
      Klocal(TI,   TJp3) += 1e0 * N[jj] * dNdx[ii][0] * detJ * weight;
      Klocal(TIp1, TJp3) += 1e0 * N[jj] * dNdx[ii][1] * detJ * weight;
      Klocal(TIp2, TJp3) += 1e0 * N[jj] * dNdx[ii][2] * detJ * weight;

      // continuity equation
      Klocal(TIp3, TJ)   += 1e0 * N[ii] * dNdx[jj][0] * detJ * weight;
      Klocal(TIp3, TJp1) += 1e0 * N[ii] * dNdx[jj][1] * detJ * weight;
      Klocal(TIp3, TJp2) += 1e0 * N[ii] * dNdx[jj][2] * detJ * weight;
    }
  }
}



void FEM::PSPG_inGaussIntegral(MatrixXd &Klocal, VectorXd &Flocal,vector<vector<double>> &dNdr,vector<vector<double>> &x_current,const int numOfNodeInElm,const double weight,const int ic)
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
  /*
  for(int p=0;p<numOfNodeInElm;p++){
    for(int q=0;q<8;q++){
      for(int j=0;j<3;j++){
        LHS_up[ic](p,q,j) -= 1e0 * N(q) * dNdx(p,j) * detJ * weight;
        LHS_pu[ic](p,q,j) += 1e0 * N(q) * dNdx(p,j) * detJ * weight;
      }
    }
  }
  */
  double h = sqrt(dx*dx+dy*dy+dz*dz);
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
      //K[ii][jj] = 0e0;
      for(int k=0;k<3;k++){
        K[ii][jj] += dNdx[ii][k]*dNdx[jj][k];
        //Klocal(TIp3, TJ)   += tau * dNdx[ii][k]*dNdx[jj][k] * detJ * weight;
        //Klocal(TIp3, TJp1) += tau * dNdx[ii][k]*dNdx[jj][k] * detJ * weight;
        //Klocal(TIp3, TJp2) += tau * dNdx[ii][k]*dNdx[jj][k] * detJ * weight;
        //Klocal(TIp3, TJp3) += tau * dNdx[ii][k]*dNdx[jj][k] * detJ * weight;
      }
      Klocal(TIp3, TJ)   += tau * K[ii][jj] * detJ * weight;
      Klocal(TIp3, TJp1) += tau * K[ii][jj] * detJ * weight;
      Klocal(TIp3, TJp2) += tau * K[ii][jj] * detJ * weight;
      Klocal(TIp3, TJp3) += tau * K[ii][jj] * detJ * weight;
    }
  }
}
