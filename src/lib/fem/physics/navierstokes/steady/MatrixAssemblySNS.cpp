#include "FEM.h"

void FEM::MatAssySNS(const int ic,MatrixXd &Klocal, VectorXd &Flocal)
{
  int ii, jj;
  int IU,IV,IW,IP;
  int JU,JV,JW,JP;

  int GP = 2;
  
  VDOUBLE2D x_current(numOfNodeInElm, VDOUBLE1D(3,0e0));
  
  VDOUBLE1D N(numOfNodeInElm, 0e0);
  VDOUBLE2D dNdr(numOfNodeInElm, VDOUBLE1D(3,0e0));
  VDOUBLE2D dNdx(numOfNodeInElm, VDOUBLE1D(3,0e0));

  VDOUBLE2D K(numOfNodeInElm, VDOUBLE1D(numOfNodeInElm, 0e0));
  
  for(int i=0;i<numOfNodeInElm;i++){
    for(int j=0;j<3;j++){
      x_current[i][j] = xFluid[elmFluid[ic]->nodeNumsPrevFluid[i]][j];
    }
  }

  double dxdr[3][3];
  double detJ, weight;
  
  Gauss gauss(GP);
  
  for(int i1=0; i1<GP; i1++){
    for(int i2=0; i2<GP; i2++){
      for(int i3=0; i3<GP; i3++){
        
        ShapeFunction3D::C3D8_N(N,gauss.point[i1],gauss.point[i2],gauss.point[i3]);
        ShapeFunction3D::C3D8_dNdr(dNdr,gauss.point[i1],gauss.point[i2],gauss.point[i3]);
        
        MathFEM::calc_dxdr(dxdr,dNdr,x_current,numOfNodeInElm);
            
        detJ = dxdr[0][0]*dxdr[1][1]*dxdr[2][2]+dxdr[0][1]*dxdr[1][2]*dxdr[2][0]+dxdr[0][2]*dxdr[1][0]*dxdr[2][1]
              -dxdr[0][2]*dxdr[1][1]*dxdr[2][0]-dxdr[0][1]*dxdr[1][0]*dxdr[2][2]-dxdr[0][0]*dxdr[1][2]*dxdr[2][1];
        weight = gauss.weight[i1] * gauss.weight[i2] * gauss.weight[i3];
              
        MathFEM::calc_dNdx(dNdx,dNdr,dxdr,numOfNodeInElm);

        double vel[3]={0e0,0e0,0e0};
        double dvdx[3][3];
        double vdvdx[3]={0e0,0e0,0e0};

        for(int p=0;p<numOfNodeInElm;p++){
          vel[0] += N[p] * uFluid[elmFluid[ic]->nodeNumsPrevFluid[p]];
          vel[1] += N[p] * vFluid[elmFluid[ic]->nodeNumsPrevFluid[p]];
          vel[2] += N[p] * wFluid[elmFluid[ic]->nodeNumsPrevFluid[p]];
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

        VDOUBLE1D tmp(numOfNodeInElm,0e0);

        for(int p=0;p<numOfNodeInElm;p++){
          for(int kk=0;kk<3;kk++){
            tmp[p] += dNdx[p][kk] * vel[kk];
          } 
        }

        double tau = calc_tau2(vel);

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

            // DIFUSION TERM 
            Klocal(IU, JU) += K[ii][jj] / Re * detJ * weight;
            Klocal(IV, JV) += K[ii][jj] / Re * detJ * weight;
            Klocal(IW, JW) += K[ii][jj] / Re * detJ * weight;

            // ADVECTION TERM 
            Klocal(IU, JU) += N[ii] * N[jj] * dvdx[0][0] * detJ * weight;
            Klocal(IU, JV) += N[ii] * N[jj] * dvdx[0][1] * detJ * weight;
            Klocal(IU, JW) += N[ii] * N[jj] * dvdx[0][2] * detJ * weight;
            Klocal(IV, JU) += N[ii] * N[jj] * dvdx[1][0] * detJ * weight;
            Klocal(IV, JV) += N[ii] * N[jj] * dvdx[1][1] * detJ * weight;
            Klocal(IV, JW) += N[ii] * N[jj] * dvdx[1][2] * detJ * weight;
            Klocal(IW, JU) += N[ii] * N[jj] * dvdx[2][0] * detJ * weight;
            Klocal(IW, JV) += N[ii] * N[jj] * dvdx[2][1] * detJ * weight;
            Klocal(IW, JW) += N[ii] * N[jj] * dvdx[2][2] * detJ * weight;
    
            Klocal(IU, JU) += N[ii] * tmp[jj] * detJ * weight;
            Klocal(IV, JV) += N[ii] * tmp[jj] * detJ * weight;
            Klocal(IW, JW) += N[ii] * tmp[jj] * detJ * weight;
            
            // PRESSURE TERM
            Klocal(IU, JP) -= N[jj] * dNdx[ii][0] * detJ * weight;
            Klocal(IV, JP) -= N[jj] * dNdx[ii][1] * detJ * weight;
            Klocal(IW, JP) -= N[jj] * dNdx[ii][2] * detJ * weight;
  
            // CONTINUITY TERM
            Klocal(IP, JU) += N[ii] * dNdx[jj][0] * detJ * weight;
            Klocal(IP, JV) += N[ii] * dNdx[jj][1] * detJ * weight;
            Klocal(IP, JW) += N[ii] * dNdx[jj][2] * detJ * weight;

            
            // SUPG PRESSURE TERM
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
            
            // SUPG ADVECTION TERM
            for(int kk=0;kk<3;kk++){
              Klocal(IU, JU) += tau * N[jj] * dNdx[ii][0] * vel[kk] * dvdx[0][kk] * detJ * weight;
              Klocal(IU, JV) += tau * N[jj] * dNdx[ii][1] * vel[kk] * dvdx[0][kk] * detJ * weight;
              Klocal(IU, JW) += tau * N[jj] * dNdx[ii][2] * vel[kk] * dvdx[0][kk] * detJ * weight;
              Klocal(IV, JU) += tau * N[jj] * dNdx[ii][0] * vel[kk] * dvdx[1][kk] * detJ * weight;
              Klocal(IV, JV) += tau * N[jj] * dNdx[ii][1] * vel[kk] * dvdx[1][kk] * detJ * weight;
              Klocal(IV, JW) += tau * N[jj] * dNdx[ii][2] * vel[kk] * dvdx[1][kk] * detJ * weight;
              Klocal(IW, JU) += tau * N[jj] * dNdx[ii][0] * vel[kk] * dvdx[2][kk] * detJ * weight;
              Klocal(IW, JV) += tau * N[jj] * dNdx[ii][1] * vel[kk] * dvdx[2][kk] * detJ * weight;
              Klocal(IW, JW) += tau * N[jj] * dNdx[ii][2] * vel[kk] * dvdx[2][kk] * detJ * weight;
            }
            Klocal(IU, JU) += tau * tmp[ii] * N[jj] * dvdx[0][0] * detJ * weight;
            Klocal(IU, JV) += tau * tmp[ii] * N[jj] * dvdx[0][1] * detJ * weight;
            Klocal(IU, JW) += tau * tmp[ii] * N[jj] * dvdx[0][2] * detJ * weight;
            Klocal(IV, JU) += tau * tmp[ii] * N[jj] * dvdx[1][0] * detJ * weight;
            Klocal(IV, JV) += tau * tmp[ii] * N[jj] * dvdx[1][1] * detJ * weight;
            Klocal(IV, JW) += tau * tmp[ii] * N[jj] * dvdx[1][2] * detJ * weight;
            Klocal(IW, JU) += tau * tmp[ii] * N[jj] * dvdx[2][0] * detJ * weight;
            Klocal(IW, JV) += tau * tmp[ii] * N[jj] * dvdx[2][1] * detJ * weight;
            Klocal(IW, JW) += tau * tmp[ii] * N[jj] * dvdx[2][2] * detJ * weight;

            Klocal(IU, JU) += tau * tmp[ii] * tmp[jj] * detJ * weight;
            Klocal(IV, JV) += tau * tmp[ii] * tmp[jj] * detJ * weight;
            Klocal(IW, JW) += tau * tmp[ii] * tmp[jj] * detJ * weight;
  
            // PSPG TERM
            Klocal(IP, JP) += tau * K[ii][jj] * detJ * weight;
            for(int k=0;k<3;k++){
              Klocal(IP, JU) +=  tau * (dNdx[ii][k]*N[jj]*dvdx[k][0]+dNdx[ii][0]*vel[k]*dNdx[jj][k]) * detJ * weight; 
              Klocal(IP, JV) +=  tau * (dNdx[ii][k]*N[jj]*dvdx[k][1]+dNdx[ii][1]*vel[k]*dNdx[jj][k]) * detJ * weight; 
              Klocal(IP, JW) +=  tau * (dNdx[ii][k]*N[jj]*dvdx[k][2]+dNdx[ii][2]*vel[k]*dNdx[jj][k]) * detJ * weight; 
            }
            
          } // II loop //

          // DIFUSION TERM
          for(int kk=0;kk<3;kk++){
            Flocal(IU) -= dNdx[ii][kk] * dvdx[0][kk] / Re * detJ * weight;
            Flocal(IV) -= dNdx[ii][kk] * dvdx[1][kk] / Re * detJ * weight;
            Flocal(IW) -= dNdx[ii][kk] * dvdx[2][kk] / Re * detJ * weight;
          }

          // ADVECTION TERM
          Flocal[IU] -= N[ii] * vdvdx[0] * detJ * weight;
          Flocal[IV] -= N[ii] * vdvdx[1] * detJ * weight;
          Flocal[IW] -= N[ii] * vdvdx[2] * detJ * weight;
            
          // PRESSURE TERM 
          Flocal[IU] += pre * dNdx[ii][0] * detJ * weight;
          Flocal[IV] += pre * dNdx[ii][1] * detJ * weight;
          Flocal[IW] += pre * dNdx[ii][2] * detJ * weight;
  
          Flocal[IP] -= N[ii] * div * detJ * weight;

          // SUPG PRESSURE TERM
          Flocal[IU] -= tau * tmp[ii] * dpdx[0] * detJ * weight; 
          Flocal[IV] -= tau * tmp[ii] * dpdx[1] * detJ * weight;
          Flocal[IW] -= tau * tmp[ii] * dpdx[2] * detJ * weight;

          // SUPG ADVECTION TERM
          for(int kk=0;kk<3;kk++){
            Flocal[IU] -= tau * tmp[ii] * vel[kk]*dvdx[0][kk] * detJ * weight;
            Flocal[IV] -= tau * tmp[ii] * vel[kk]*dvdx[1][kk] * detJ * weight;
            Flocal[IW] -= tau * tmp[ii] * vel[kk]*dvdx[2][kk] * detJ * weight;
          }

          // PSPG TERM
          Flocal[IP] -= tau * dNdx[ii][0] * dpdx[0] * detJ * weight;
          Flocal[IP] -= tau * dNdx[ii][1] * dpdx[1] * detJ * weight;
          Flocal[IP] -= tau * dNdx[ii][2] * dpdx[2] * detJ * weight;
          for(int kk=0;kk<3;kk++){
            Flocal[IP] -= tau * dNdx[ii][0] * vel[kk] * dvdx[0][kk] * detJ * weight; 
            Flocal[IP] -= tau * dNdx[ii][1] * vel[kk] * dvdx[1][kk] * detJ * weight;           
            Flocal[IP] -= tau * dNdx[ii][2] * vel[kk] * dvdx[2][kk] * detJ * weight; 
          }
        }
      }
    }
  }
}