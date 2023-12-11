#include "FEM.h"
using namespace std;

void FEM::SteadyNavierStokesMatrix(const int ic,MatrixXd &Klocal, VectorXd &Flocal)
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

  Gauss g1(1),g2(2),g3(3);

  double dxdr[3][3];


  double detJ, weight;
  
  for(int i1=0; i1<2; i1++){
    for(int i2=0; i2<2; i2++){
      for(int i3=0; i3<2; i3++){
        ShapeFunction3D::C3D8_N(N,g1.point[i1],g1.point[i2],g1.point[i3]);
        ShapeFunction3D::C3D8_dNdr(dNdr,g1.point[i1],g1.point[i2],g1.point[i3]);
        
        MathFEM::calc_dxdr(dxdr,dNdr,x_current,numOfNodeInElm);
            
        detJ = dxdr[0][0]*dxdr[1][1]*dxdr[2][2]+dxdr[0][1]*dxdr[1][2]*dxdr[2][0]+dxdr[0][2]*dxdr[1][0]*dxdr[2][1]
                          -dxdr[0][2]*dxdr[1][1]*dxdr[2][0]-dxdr[0][1]*dxdr[1][0]*dxdr[2][2]-dxdr[0][0]*dxdr[1][2]*dxdr[2][1];
        weight = g1.weight[i1]*g1.weight[i2]*g1.weight[i3];
              
        MathFEM::calc_dNdx(dNdx,dNdr,dxdr,numOfNodeInElm);
        MathFEM::calc_dNdx(dNdx,dNdr,dxdr,numOfNodeInElm);

        double vel[3]={0e0,0e0,0e0};
        double dvdx[3][3];
        double vdvdx[3]={0e0,0e0,0e0};

        for(int p=0;p<numOfNodeInElm;p++){
          vel[0] += N[p]*uFluid[elmFluid[ic]->nodeNumsPrevFluid[p]];
          vel[1] += N[p]*vFluid[elmFluid[ic]->nodeNumsPrevFluid[p]];
          vel[2] += N[p]*pFluid[elmFluid[ic]->nodeNumsPrevFluid[p]];
        }

        for(int j=0;j<3;j++){
          dvdx[0][j] = 0e0;
          dvdx[1][j] = 0e0;
          dvdx[2][j] = 0e0;
          for(int p=0;p<numOfNodeInElm;p++){
            dvdx[0][j] += dNdx[p][j] * uFluid[elmFluid[ic]->nodeNumsPrevFluid[p]];
            dvdx[1][j] += dNdx[p][j] * vFluid[elmFluid[ic]->nodeNumsPrevFluid[p]];
            dvdx[2][j] += dNdx[p][j] * wFluid[elmFluid[ic]->nodeNumsPrevFluid[p]];
          }
        }

        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++) vdvdx[i] += vel[j]*dvdx[i][j];
        }

        double pre=0e0;
        double dpdx[3]={0e0,0e0,0e0};
        
        for(int p=0;p<numOfNodeInElm;p++){
          pre += N[p]*pFluid[elmFluid[ic]->nodeNumsPrevFluid[p]];    
        } 
        for(int j=0;j<3;j++){
          for(int p=0;p<numOfNodeInElm;p++){
            dpdx[j] += dNdx[p][j]*pFluid[elmFluid[ic]->nodeNumsPrevFluid[p]];
          }
        }

        double dudx=0e0,dvdy=0e0,dwdz=0e0;
        for(int p=0;p<numOfNodeInElm;p++){
          dudx += dNdx[p][0] * uFluid[elmFluid[ic]->nodeNumsPrevFluid[p]];
          dvdy += dNdx[p][1] * vFluid[elmFluid[ic]->nodeNumsPrevFluid[p]];
          dwdz += dNdx[p][2] * wFluid[elmFluid[ic]->nodeNumsPrevFluid[p]];
        }
        double div = dudx + dvdy + dwdz;

        vector<double> tmp(numOfNodeInElm,0e0);

        for(int p=0;p<numOfNodeInElm;p++){
          for(int kk=0;kk<3;kk++){
            tmp[p] += dNdx[p][kk] * vel[kk];
          } 
        }

        //// stabilization parameter ////
        double h = sqrt(dx*dx+dy*dy+dz*dz);
        double vMag = sqrt(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2]);
        double tau=pow(2e0*vMag/h,2e0)+9e0*pow(4e0*nu/(h*h),2e0);
        tau = 1e0/sqrt(tau);

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

            /// advection ///
            Klocal(IU, JU) += rho * N[ii] * N[jj] * dvdx[0][0] * detJ * weight;
            Klocal(IU, JV) += rho * N[ii] * N[jj] * dvdx[0][1] * detJ * weight;
            Klocal(IU, JW) += rho * N[ii] * N[jj] * dvdx[0][2] * detJ * weight;
            Klocal(IV, JU) += rho * N[ii] * N[jj] * dvdx[1][0] * detJ * weight;
            Klocal(IV, JV) += rho * N[ii] * N[jj] * dvdx[1][1] * detJ * weight;
            Klocal(IV, JW) += rho * N[ii] * N[jj] * dvdx[1][2] * detJ * weight;
            Klocal(IW, JU) += rho * N[ii] * N[jj] * dvdx[2][0] * detJ * weight;
            Klocal(IW, JV) += rho * N[ii] * N[jj] * dvdx[2][1] * detJ * weight;
            Klocal(IW, JW) += rho * N[ii] * N[jj] * dvdx[2][2] * detJ * weight;
            
            Klocal(IU, JU) += rho * N[ii]*tmp[jj] * detJ * weight;
            Klocal(IV, JV) += rho * N[ii]*tmp[jj] * detJ * weight;
            Klocal(IW, JW) += rho * N[ii]*tmp[jj] * detJ * weight;
            
            //// pressure ////
            Klocal(IU, JP) += N[jj] * dNdx[ii][0] * detJ * weight;
            Klocal(IV, JP) += N[jj] * dNdx[ii][1] * detJ * weight;
            Klocal(IW, JP) += N[jj] * dNdx[ii][2] * detJ * weight;
        
            //// continuity ////
            Klocal(IP, JU) += N[ii] * dNdx[jj][0] * detJ * weight;
            Klocal(IP, JV) += N[ii] * dNdx[jj][1] * detJ * weight;
            Klocal(IP, JW) += N[ii] * dNdx[jj][2] * detJ * weight;
            
            //// SUPG ////

            //pressure term
            Klocal(IU, JU) += tau * dNdx[ii][0] * N[jj] * dpdx[0] * detJ * weight;
            Klocal(IU, JV) += tau * dNdx[ii][1] * N[jj] * dpdx[0] * detJ * weight;
            Klocal(IU, JW) += tau * dNdx[ii][2] * N[jj] * dpdx[0] * detJ * weight;
            Klocal(IV, JU) += tau * dNdx[ii][0] * N[jj] * dpdx[1] * detJ * weight;
            Klocal(IV, JV) += tau * dNdx[ii][1] * N[jj] * dpdx[1] * detJ * weight;
            Klocal(IV, JW) += tau * dNdx[ii][2] * N[jj] * dpdx[1] * detJ * weight;
            Klocal(IW, JU) += tau * dNdx[ii][0] * N[jj] * dpdx[2] * detJ * weight;
            Klocal(IW, JV) += tau * dNdx[ii][1] * N[jj] * dpdx[2] * detJ * weight;
            Klocal(IW, JW) += tau * dNdx[ii][2] * N[jj] * dpdx[2] * detJ * weight;
            
            Klocal(IU, JP) += tau * tmp[ii] * dNdx[jj][0] * detJ * weight;
            Klocal(IV, JP) += tau * tmp[ii] * dNdx[jj][1] * detJ * weight;
            Klocal(IW, JP) += tau * tmp[ii] * dNdx[jj][2] * detJ * weight;
            
            //advection term
            for(int kk=0;kk<3;kk++){
              Klocal(IU, JU) += rho * tau * N[jj] * dNdx[ii][0] * vel[kk] * dvdx[0][kk] * detJ * weight;
              Klocal(IU, JV) += rho * tau * N[jj] * dNdx[ii][1] * vel[kk] * dvdx[0][kk] * detJ * weight;
              Klocal(IU, JW) += rho * tau * N[jj] * dNdx[ii][2] * vel[kk] * dvdx[0][kk] * detJ * weight;
              Klocal(IV, JU) += rho * tau * N[jj] * dNdx[ii][0] * vel[kk] * dvdx[1][kk] * detJ * weight;
              Klocal(IV, JV) += rho * tau * N[jj] * dNdx[ii][1] * vel[kk] * dvdx[1][kk] * detJ * weight;
              Klocal(IV, JW) += rho * tau * N[jj] * dNdx[ii][2] * vel[kk] * dvdx[1][kk] * detJ * weight;
              Klocal(IW, JU) += rho * tau * N[jj] * dNdx[ii][0] * vel[kk] * dvdx[2][kk] * detJ * weight;
              Klocal(IW, JV) += rho * tau * N[jj] * dNdx[ii][1] * vel[kk] * dvdx[2][kk] * detJ * weight;
              Klocal(IW, JW) += rho * tau * N[jj] * dNdx[ii][2] * vel[kk] * dvdx[2][kk] * detJ * weight;
            }
            Klocal(IU, JU) += rho * tau * tmp[ii] * N[jj] * dvdx[0][0] * detJ * weight;
            Klocal(IU, JV) += rho * tau * tmp[ii] * N[jj] * dvdx[0][1] * detJ * weight;
            Klocal(IU, JW) += rho * tau * tmp[ii] * N[jj] * dvdx[0][2] * detJ * weight;
            Klocal(IV, JU) += rho * tau * tmp[ii] * N[jj] * dvdx[1][0] * detJ * weight;
            Klocal(IV, JV) += rho * tau * tmp[ii] * N[jj] * dvdx[1][1] * detJ * weight;
            Klocal(IV, JW) += rho * tau * tmp[ii] * N[jj] * dvdx[1][2] * detJ * weight;
            Klocal(IW, JU) += rho * tau * tmp[ii] * N[jj] * dvdx[2][0] * detJ * weight;
            Klocal(IW, JV) += rho * tau * tmp[ii] * N[jj] * dvdx[2][1] * detJ * weight;
            Klocal(IW, JW) += rho * tau * tmp[ii] * N[jj] * dvdx[2][2] * detJ * weight;

            Klocal(IU, JU) += rho * tau * tmp[ii] * tmp[jj] * detJ * weight;
            Klocal(IV, JV) += rho * tau * tmp[ii] * tmp[jj] * detJ * weight;
            Klocal(IW, JW) += rho * tau * tmp[ii] * tmp[jj] * detJ * weight;
  

            //// PSPG ////
            Klocal(IP, JP) += (tau/rho) * K[ii][jj] * detJ * weight;
            for(int k=0;k<3;k++){
              Klocal(IP, JU) =  tau * (dNdx[ii][k]*N[jj]*dvdx[k][0]+dNdx[ii][0]*vel[k]*dNdx[jj][k]) * detJ * weight; 
              Klocal(IP, JV) =  tau * (dNdx[ii][k]*N[jj]*dvdx[k][1]+dNdx[ii][1]*vel[k]*dNdx[jj][k]) * detJ * weight; 
              Klocal(IP, JW) =  tau * (dNdx[ii][k]*N[jj]*dvdx[k][2]+dNdx[ii][2]*vel[k]*dNdx[jj][k]) * detJ * weight; 
            }
          }

          //// disfusion ////
          for(int kk=0;kk<3;kk++){
            Flocal(IU) += mu * dNdx[ii][kk] * dvdx[0][kk] * detJ * weight;
            Flocal(IV) += mu * dNdx[ii][kk] * dvdx[1][kk] * detJ * weight;
            Flocal(IW) += mu * dNdx[ii][kk] * dvdx[2][kk] * detJ * weight;
          }

          /// advection ///
          Flocal[IU] += rho * N[ii]*vdvdx[0] * detJ * weight;
          Flocal[IV] += rho * N[ii]*vdvdx[1] * detJ * weight;
          Flocal[IW] += rho * N[ii]*vdvdx[2] * detJ * weight;
            
          //// pressure ////
          Flocal[IU] -= 1e0 * pre * dNdx[ii][0] * detJ * weight;
          Flocal[IV] -= 1e0 * pre * dNdx[ii][1] * detJ * weight;
          Flocal[IW] -= 1e0 * pre * dNdx[ii][2] * detJ * weight;
  
          //// continuity ////
          Flocal[IP] += 1e0 * N[ii] * div * detJ * weight;

          //// SUPG ////
          Flocal[IU] += tau * tmp[ii] * dpdx[0] * detJ * weight; //pressure term
          Flocal[IV] += tau * tmp[ii] * dpdx[1] * detJ * weight;
          Flocal[IW] += tau * tmp[ii] * dpdx[2] * detJ * weight;
          for(int kk=0;kk<3;kk++){
            Flocal[IU] += rho* tau * tmp[ii] * vel[kk]*dvdx[0][kk] * detJ * weight; //advection
            Flocal[IV] += rho* tau * tmp[ii] * vel[kk]*dvdx[1][kk] * detJ * weight;
            Flocal[IW] += rho* tau * tmp[ii] * vel[kk]*dvdx[2][kk] * detJ * weight;
          }

          //// PSPG ////
          Flocal[IP] += (tau/rho) * dNdx[ii][0] * dpdx[0] * detJ * weight;
          Flocal[IP] += (tau/rho) * dNdx[ii][1] * dpdx[1] * detJ * weight;
          Flocal[IP] += (tau/rho) * dNdx[ii][2] * dpdx[2] * detJ * weight;
          for(int kk=0;kk<3;kk++){
            Flocal[IP] += tau * dNdx[ii][0] * vel[kk] * dvdx[0][kk] * detJ * weight; 
            Flocal[IP] += tau * dNdx[ii][1] * vel[kk] * dvdx[2][kk] * detJ * weight;           
            Flocal[IP] += tau * dNdx[ii][2] * vel[kk] * dvdx[1][kk] * detJ * weight; 
          }

        }
      }
    }
  }
}