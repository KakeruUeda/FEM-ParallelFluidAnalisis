#include "FEM.h"
using namespace std;

void FEM::XFEM_MatAssySTT(const int ic, MatrixXd &Klocal, VectorXd &Flocal)
{
  int ii, jj;
  int IU,IV,IW,IP;
  int JU,JV,JW,JP;
  
  int GP = 2;

  VDOUBLE2D x_current(numOfNodeInElm,VDOUBLE1D(3,0e0));
  VDOUBLE1D sdf_current(numOfNodeInElm,0e0);
  
  VDOUBLE1D NP(numOfNodeInElm,0e0);
  VDOUBLE1D NV(numOfNodeInElm,0e0);
  VDOUBLE2D dNPdr(numOfNodeInElm,VDOUBLE1D(3,0e0));
  VDOUBLE2D dNVdr(numOfNodeInElm,VDOUBLE1D(3,0e0));
  VDOUBLE2D dNVdx(numOfNodeInElm,VDOUBLE1D(3,0e0));
  VDOUBLE2D dNPdx(numOfNodeInElm,VDOUBLE1D(3,0e0));
  
  VDOUBLE2D K(numOfNodeInElm,VDOUBLE1D(numOfNodeInElm,0e0));
  VDOUBLE2D L(numOfNodeInElm,VDOUBLE1D(numOfNodeInElm,0e0));

  double minSDF_node = 1e12;
  double maxSDF_node = -1e12;

  double minSDF = sqrt(dx*dx+dy*dy+dz*dz)/2;
  
  int check = 0;
  for(int i=0;i<numOfNodeInElm;i++){
    sdf_current[i] = sdfFluid[elmFluid[ic]->nodeNumsPrevFluid[i]];
    if(sdf_current[i] <= 0e0){
      check = check +1;
    }
    if(sdf_current[i]>0e0 && sdf_current[i]<minSDF_node) minSDF_node = sdf_current[i];
    if(sdf_current[i]>0e0 && sdf_current[i]>maxSDF_node) maxSDF_node = sdf_current[i];
    for(int j=0;j<3;j++){
      x_current[i][j] = xFluid[elmFluid[ic]->nodeNumsPrevFluid[i]][j];
    }
  }

  VDOUBLE1D b(sub_div);
  LocalRefinement(sub_div,b);

  Gauss gauss(GP);

  double dxdr[3][3];

  //// stabilization parameter ////
  double h = dx/2;
  double tau = h*h/mu/12e0;
  
  double sdf_gauss;
  double g1_local, g2_local, g3_local;
  double detJ, weight;
  int countMinus = 0;

  //#pragma omp parallel for
  for(int n1=0;n1<sub_div;n1++){
    for(int n2=0;n2<sub_div;n2++){
      for(int n3=0;n3<sub_div;n3++){
        for(int i1=0; i1<GP; i1++){
          for(int i2=0; i2<GP; i2++){
            for(int i3=0; i3<GP; i3++){

              g1_local = (gauss.point[i1] + b[n1])/(double)sub_div;
              g2_local = (gauss.point[i2] + b[n2])/(double)sub_div;
              g3_local = (gauss.point[i3] + b[n3])/(double)sub_div;

              ShapeFunction3D::C3D8_N(NP,g1_local,g2_local,g3_local);
              ShapeFunction3D::C3D8_N(NV,g1_local,g2_local,g3_local);
              ShapeFunction3D::C3D8_dNdr(dNPdr,g1_local,g2_local,g3_local);
              ShapeFunction3D::C3D8_dNdr(dNVdr,g1_local,g2_local,g3_local); 
              
              double sdf_gauss = 0e0;
              for(int p=0;p<numOfNodeInElm;p++) sdf_gauss += NP[p] * sdf_current[p];

              if(sdf_gauss<=0){
                countMinus++; 
                //continue;
                for(int p=0;p<numOfNodeInElm;p++){
                  NV[p] = 0;
                  NP[p] = 0;
                }
              }else{
                if(sdf_gauss/minSDF<1.0){
                  double dfdr[3]={0e0,0e0,0e0};
                  for(int i=0;i<3;i++){
                    for(int p=0;p<numOfNodeInElm;p++){
                      dfdr[i] += sdf_current[p]/minSDF * dNPdr[p][i];
                    }
                  }
                  for(int p=0;p<numOfNodeInElm;p++){
                    for(int i=0;i<3;i++) dNVdr[p][i] = dNPdr[p][i] * sdf_gauss/minSDF + NP[p]*dfdr[i];
                    NV[p] = NP[p] * sdf_gauss/minSDF;       
                  }
                }
              }

              
              MathFEM::calc_dxdr(dxdr,dNPdr,x_current,numOfNodeInElm);
            
              detJ = dxdr[0][0]*dxdr[1][1]*dxdr[2][2]+dxdr[0][1]*dxdr[1][2]*dxdr[2][0]+dxdr[0][2]*dxdr[1][0]*dxdr[2][1]
                          -dxdr[0][2]*dxdr[1][1]*dxdr[2][0]-dxdr[0][1]*dxdr[1][0]*dxdr[2][2]-dxdr[0][0]*dxdr[1][2]*dxdr[2][1];
              weight = gauss.weight[i1]*gauss.weight[i2]*gauss.weight[i3]/(sub_div*sub_div*sub_div);
              
              MathFEM::calc_dNdx(dNVdx,dNVdr,dxdr,numOfNodeInElm);
              MathFEM::calc_dNdx(dNPdx,dNPdr,dxdr,numOfNodeInElm);

              if(sdf_gauss<=0){
                for(int p=0;p<numOfNodeInElm;p++){
                  for(int i=0;i<3;i++){
                    dNVdx[p][i] = 0;
                  }
                }
              }
              
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
                  L[ii][jj] = 0e0;
                  
                  for(int k=0;k<3;k++){
                    K[ii][jj] += dNVdx[ii][k] * dNVdx[jj][k];
                    L[ii][jj] += dNPdx[ii][k] * dNPdx[jj][k];
                  }

                  //// disfusion ////
                  Klocal(IU, JU) -= mu * K[ii][jj] * detJ * weight;
                  Klocal(IV, JV) -= mu * K[ii][jj] * detJ * weight;
                  Klocal(IW, JW) -= mu * K[ii][jj] * detJ * weight;

                  //// pressure ////
                  Klocal(IU, JP) += NP[jj] * dNVdx[ii][0] * detJ * weight;
                  Klocal(IV, JP) += NP[jj] * dNVdx[ii][1] * detJ * weight;
                  Klocal(IW, JP) += NP[jj] * dNVdx[ii][2] * detJ * weight;
              
                  //// continuity ////
                  Klocal(IP, JU) += NP[ii] * dNVdx[jj][0] * detJ * weight;
                  Klocal(IP, JV) += NP[ii] * dNVdx[jj][1] * detJ * weight;
                  Klocal(IP, JW) += NP[ii] * dNVdx[jj][2] * detJ * weight;
                  
                  //// PSPG ////
                  Klocal(IP, JP) += tau * L[ii][jj] * detJ * weight;
                }
              }
            }
          }
        }
      }
    }
  }
  
}


void FEM::XFEM_MatAssySTT2(const int ic, MatrixXd &Klocal, VectorXd &Flocal)
{
  int ii, jj;
  int IU, IV, IW, IP;
  int JU, JV, JW, JP;

  VDOUBLE2D x_current(numOfNodeInElm,VDOUBLE1D(3,0e0));
  VDOUBLE1D sdf_current(numOfNodeInElm,0e0);
  
  VDOUBLE1D NP(numOfNodeInElm,0e0);
  VDOUBLE1D NV(numOfNodeInElm,0e0);
  VDOUBLE2D dNPdr(numOfNodeInElm,VDOUBLE1D(3,0e0));
  VDOUBLE2D dNVdr(numOfNodeInElm,VDOUBLE1D(3,0e0));
  VDOUBLE2D dNVdx(numOfNodeInElm,VDOUBLE1D(3,0e0));
  VDOUBLE2D dNPdx(numOfNodeInElm,VDOUBLE1D(3,0e0));
  
  VDOUBLE2D K(numOfNodeInElm,VDOUBLE1D(numOfNodeInElm,0e0));
  VDOUBLE2D L(numOfNodeInElm,VDOUBLE1D(numOfNodeInElm,0e0));

  double minSDF_node = 1e12;
  double minSDF = sqrt(dx*dx+dy*dy+dz*dz)/2e0;
  //double minSDF = dx/2e0;

  for(int i=0;i<numOfNodeInElm;i++){
    for(int j=0;j<3;j++){
      x_current[i][j] = xFluid[elmFluid[ic]->nodeNumsPrevFluid[i]][j];
    }
  }
  for(int i=0; i<numOfNodeInElm; i++){
    sdf_current[i] = sdfFluid[elmFluid[ic]->nodeNumsPrevFluid[i]];
  }

  double dxdr[3][3];

  //// stabilization parameter ////
  double h = dx/2;
  double tau = h*h/mu/12e0;
  
  double sdf_gauss;
  double detJ, weight;
  int countMinus = 0;

  int cc = 0;

  int GP = 2;
  int totalGP = GP * GP * GP;

  for(auto itr = elmFluid[ic]->sub_elm.begin(); itr != elmFluid[ic]->sub_elm.end(); itr++)
  {
    for(int gp=0; gp<totalGP; gp++)
    {
      double g1_local = elmFluid[ic]->sub_elm[cc].sub_gx[gp];
      double g2_local = elmFluid[ic]->sub_elm[cc].sub_gy[gp];
      double g3_local = elmFluid[ic]->sub_elm[cc].sub_gz[gp];
      
      ShapeFunction3D::C3D8_N(NP,g1_local,g2_local,g3_local);
      ShapeFunction3D::C3D8_N(NV,g1_local,g2_local,g3_local);
      ShapeFunction3D::C3D8_dNdr(dNPdr,g1_local,g2_local,g3_local);
      ShapeFunction3D::C3D8_dNdr(dNVdr,g1_local,g2_local,g3_local); 

      double sdf_gauss = 0e0;
      
      for(int p=0;p<numOfNodeInElm;p++){
        sdf_gauss += NP[p] * sdf_current[p];
      }

      if(sdf_gauss<=0){
        //continue;
        countMinus++; 
        for(int p=0;p<numOfNodeInElm;p++){
          NV[p] = 0;
          NP[p] = 0;
        }
      }else{
        if(sdf_gauss/minSDF<1.0){
          double dfdr[3]={0e0,0e0,0e0};
          for(int i=0;i<3;i++){
            for(int p=0;p<numOfNodeInElm;p++){
              dfdr[i] += sdf_current[p]/minSDF * dNPdr[p][i];
            }
          }
          for(int p=0;p<numOfNodeInElm;p++){
            for(int i=0;i<3;i++) dNVdr[p][i] = dNPdr[p][i] * sdf_gauss/minSDF + NP[p]*dfdr[i];
            NV[p] = NP[p] * sdf_gauss/minSDF;            
          }
        }
      }

      MathFEM::calc_dxdr(dxdr,dNPdr,x_current,numOfNodeInElm);
        
      detJ = dxdr[0][0]*dxdr[1][1]*dxdr[2][2]+dxdr[0][1]*dxdr[1][2]*dxdr[2][0]+dxdr[0][2]*dxdr[1][0]*dxdr[2][1]
            -dxdr[0][2]*dxdr[1][1]*dxdr[2][0]-dxdr[0][1]*dxdr[1][0]*dxdr[2][2]-dxdr[0][0]*dxdr[1][2]*dxdr[2][1];
      weight = elmFluid[ic]->sub_elm[cc].sub_weight[gp];
      
      MathFEM::calc_dNdx(dNVdx,dNVdr,dxdr,numOfNodeInElm);
      MathFEM::calc_dNdx(dNPdx,dNPdr,dxdr,numOfNodeInElm);
      
      if(sdf_gauss<=0){
        for(int p=0;p<numOfNodeInElm;p++){
          for(int i=0;i<3;i++){
            dNVdx[p][i] = 0;
          }
        }
      }

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
          L[ii][jj] = 0e0;
          
          for(int k=0;k<3;k++){
            K[ii][jj] += dNVdx[ii][k] * dNVdx[jj][k];
            L[ii][jj] += dNPdx[ii][k] * dNPdx[jj][k];
          }

          //// disfusion ////
          Klocal(IU, JU) -= mu * K[ii][jj] * detJ * weight;
          Klocal(IV, JV) -= mu * K[ii][jj] * detJ * weight;
          Klocal(IW, JW) -= mu * K[ii][jj] * detJ * weight;

          //// pressure ////
          Klocal(IU, JP) += NP[jj] * dNVdx[ii][0] * detJ * weight;
          Klocal(IV, JP) += NP[jj] * dNVdx[ii][1] * detJ * weight;
          Klocal(IW, JP) += NP[jj] * dNVdx[ii][2] * detJ * weight;
      
          //// continuity ////
          Klocal(IP, JU) += NP[ii] * dNVdx[jj][0] * detJ * weight;
          Klocal(IP, JV) += NP[ii] * dNVdx[jj][1] * detJ * weight;
          Klocal(IP, JW) += NP[ii] * dNVdx[jj][2] * detJ * weight;
          
          //// PSPG ////
          Klocal(IP, JP) += tau * L[ii][jj] * detJ * weight;
        }
      }
    }
    cc++;
  }

}




void FEM::SAWADA_XFEM_MatAssySTT(const int ic, MatrixXd &Klocal, VectorXd &Flocal)
{
  int ii, jj;
  int IU,IV,IW,IP;
  int JU,JV,JW,JP;

  VDOUBLE2D x_current(numOfNodeInElm,VDOUBLE1D(3,0e0));
  VDOUBLE1D sdf_current(numOfNodeInElm,0e0);
  
  VDOUBLE1D NP(numOfNodeInElm,0e0);
  VDOUBLE1D NV(numOfNodeInElm,0e0);
  VDOUBLE2D dNPdr(numOfNodeInElm,VDOUBLE1D(3,0e0));
  VDOUBLE2D dNVdr(numOfNodeInElm,VDOUBLE1D(3,0e0));
  VDOUBLE2D dNVdx(numOfNodeInElm,VDOUBLE1D(3,0e0));
  VDOUBLE2D dNPdx(numOfNodeInElm,VDOUBLE1D(3,0e0));
  
  VDOUBLE2D K(numOfNodeInElm,VDOUBLE1D(numOfNodeInElm,0e0));
  VDOUBLE2D L(numOfNodeInElm,VDOUBLE1D(numOfNodeInElm,0e0));

  double minSDF_node = 1e12;
  double minSDF = dx/2;
  
  int check = 0;
  for(int i=0;i<numOfNodeInElm;i++){
    sdf_current[i] = sdfFluid[elmFluid[ic]->nodeNumsPrevFluid[i]];
    if(sdf_current[i] <= 0e0){
      check = check +1;
    }
    if(sdf_current[i]>0e0 && sdf_current[i]<minSDF_node) minSDF_node = sdf_current[i];
    for(int j=0;j<3;j++){
      x_current[i][j] = xFluid[elmFluid[ic]->nodeNumsPrevFluid[i]][j];
    }
  }

  if((minSDF_node < dx/20) && (check == 3)){
    //return;
  }

  bool flag = false;
  if(check==8){
    flag = true;
    cout << "**** ERROR ****" << endl;
  }

  Gauss g1(1),g2(2),g3(3),g5(6),g6(6),g19(19),g30(30);

  double dxdr[3][3];

  //// stabilization parameter ////
  double h = dx/2;
  double tau = h*h/mu/12e0;
  
  double sdf_gauss;
  double detJ, weight;
  int countMinus = 0;
  //#pragma omp parallel for
  for(int i1=0; i1<31; i1++){
    for(int i2=0; i2<31; i2++){
      for(int i3=0; i3<31; i3++){

        ShapeFunction3D::C3D8_N(NP,g30.point[i1],g30.point[i2],g30.point[i3]);
        ShapeFunction3D::C3D8_N(NV,g30.point[i1],g30.point[i2],g30.point[i2]);
        ShapeFunction3D::C3D8_dNdr(dNPdr,g30.point[i1],g30.point[i2],g30.point[i3]);
        ShapeFunction3D::C3D8_dNdr(dNVdr,g30.point[i1],g30.point[i2],g30.point[i3]); 
        
        double sdf_gauss = 0e0;
        for(int p=0;p<numOfNodeInElm;p++) sdf_gauss += NP[p] * sdf_current[p];


        if(sdf_gauss<=0){
          countMinus++; 
          for(int p=0;p<numOfNodeInElm;p++){
            NV[p] = 0;
            NP[p] = 0;
          }
        }else{
          if(sdf_gauss/minSDF<1e0){
            double dfdr[3]={0e0,0e0,0e0};
            for(int i=0;i<3;i++){
              for(int p=0;p<numOfNodeInElm;p++) dfdr[i] += sdf_current[p]/minSDF * dNPdr[p][i];
            }
            for(int p=0;p<numOfNodeInElm;p++){
              for(int i=0;i<3;i++) dNVdr[p][i] = dNPdr[p][i]*sdf_gauss/minSDF + NP[p]*dfdr[i];
              NV[p] = NP[p]*sdf_gauss/minSDF;
            }
          }
        }
        
        MathFEM::calc_dxdr(dxdr,dNPdr,x_current,numOfNodeInElm);
      
        detJ = dxdr[0][0]*dxdr[1][1]*dxdr[2][2]+dxdr[0][1]*dxdr[1][2]*dxdr[2][0]+dxdr[0][2]*dxdr[1][0]*dxdr[2][1]
               -dxdr[0][2]*dxdr[1][1]*dxdr[2][0]-dxdr[0][1]*dxdr[1][0]*dxdr[2][2]-dxdr[0][0]*dxdr[1][2]*dxdr[2][1];
        weight = g30.weight[i1]*g30.weight[i2]*g30.weight[i3];
        
        MathFEM::calc_dNdx(dNVdx,dNVdr,dxdr,numOfNodeInElm);
        MathFEM::calc_dNdx(dNPdx,dNPdr,dxdr,numOfNodeInElm);

        if(sdf_gauss<=0){
          for(int p=0;p<numOfNodeInElm;p++){
            for(int i=0;i<3;i++){
              dNVdx[p][i] = 0;
            }
          }
        }
        
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
            L[ii][jj] = 0e0;
            
            for(int k=0;k<3;k++){
              K[ii][jj] += dNVdx[ii][k]*dNVdx[jj][k];
              L[ii][jj] += dNPdx[ii][k]*dNPdx[jj][k];
            }

            //// disfusion ////
            Klocal(IU, JU) -= mu * K[ii][jj] * detJ * weight;
            Klocal(IV, JV) -= mu * K[ii][jj] * detJ * weight;
            Klocal(IW, JW) -= mu * K[ii][jj] * detJ * weight;

            //// pressure ////
            Klocal(IU, JP) += NP[jj] * dNVdx[ii][0] * detJ * weight;
            Klocal(IV, JP) += NP[jj] * dNVdx[ii][1] * detJ * weight;
            Klocal(IW, JP) += NP[jj] * dNVdx[ii][2] * detJ * weight;
        
            //// continuity ////
            Klocal(IP, JU) += NP[ii] * dNVdx[jj][0] * detJ * weight;
            Klocal(IP, JV) += NP[ii] * dNVdx[jj][1] * detJ * weight;
            Klocal(IP, JW) += NP[ii] * dNVdx[jj][2] * detJ * weight;
            
            //// PSPG ////
            Klocal(IP, JP) += tau * L[ii][jj] * detJ * weight;
          }
        }
      }
    }
  }
}



void FEM::DiffusionInGaussIntegralXFEM(MatrixXd &Klocal, VectorXd &Flocal, VDOUBLE2D &dNPdr, VDOUBLE2D &dNVdr, VDOUBLE2D &x_current, const int numOfNodeInElm,const double weight,const int ic)
{
  int ii, jj;
  int TI,TIp1,TIp2,TIp3;
  int TJ,TJp1,TJp2,TJp3;

  VDOUBLE2D dNVdx(numOfNodeInElm,VDOUBLE1D(3,0));

  double dxdr[3][3];

  MathFEM::calc_dxdr(dxdr,dNPdr,x_current,numOfNodeInElm);
  double detJ=dxdr[0][0]*dxdr[1][1]*dxdr[2][2]+dxdr[0][1]*dxdr[1][2]*dxdr[2][0]+dxdr[0][2]*dxdr[1][0]*dxdr[2][1]
             -dxdr[0][2]*dxdr[1][1]*dxdr[2][0]-dxdr[0][1]*dxdr[1][0]*dxdr[2][2]-dxdr[0][0]*dxdr[1][2]*dxdr[2][1];
  MathFEM::calc_dNdx(dNVdx,dNVdr,dxdr,numOfNodeInElm);


  VDOUBLE2D K(numOfNodeInElm,VDOUBLE1D(numOfNodeInElm,0));


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
      double fact;
      
      for(int k=0;k<3;k++){
        K[ii][jj] += dNVdx[ii][k]*dNVdx[jj][k];
      }
      Klocal(TI,   TJ)   -= mu * K[ii][jj] * detJ * weight;
      Klocal(TIp1, TJp1) -= mu * K[ii][jj] * detJ * weight;
      Klocal(TIp2, TJp2) -= mu * K[ii][jj] * detJ * weight;
    }
  }

}



void FEM::PressureInGaussIntegralXFEM(MatrixXd &Klocal, VectorXd &Flocal, VDOUBLE1D &NP, VDOUBLE2D &dNPdr, VDOUBLE2D &dNVdr, VDOUBLE2D &x_current, const int numOfNodeInElm, const double weight, const int ic)
{ 
  int ii, jj;
  int TI,TIp1,TIp2,TIp3;
  int TJ,TJp1,TJp2,TJp3;

  VDOUBLE2D dNVdx(numOfNodeInElm,VDOUBLE1D(3,0));

  double dxdr[3][3];
  MathFEM::calc_dxdr(dxdr,dNPdr,x_current,numOfNodeInElm);
  double detJ=dxdr[0][0]*dxdr[1][1]*dxdr[2][2]+dxdr[0][1]*dxdr[1][2]*dxdr[2][0]+dxdr[0][2]*dxdr[1][0]*dxdr[2][1]
             -dxdr[0][2]*dxdr[1][1]*dxdr[2][0]-dxdr[0][1]*dxdr[1][0]*dxdr[2][2]-dxdr[0][0]*dxdr[1][2]*dxdr[2][1];
  
  MathFEM::calc_dNdx(dNVdx,dNVdr,dxdr,numOfNodeInElm);

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
      
      // pressure term
      Klocal(TI,   TJp3) += 1e0 * NP[jj] * dNVdx[ii][0] * detJ * weight;
      Klocal(TIp1, TJp3) += 1e0 * NP[jj] * dNVdx[ii][1] * detJ * weight;
      Klocal(TIp2, TJp3) += 1e0 * NP[jj] * dNVdx[ii][2] * detJ * weight;

      // continuity equation
      Klocal(TIp3, TJ)   += 1e0 * NP[ii] * dNVdx[jj][0] * detJ * weight;
      Klocal(TIp3, TJp1) += 1e0 * NP[ii] * dNVdx[jj][1] * detJ * weight;
      Klocal(TIp3, TJp2) += 1e0 * NP[ii] * dNVdx[jj][2] * detJ * weight;
    }
  }
}


void FEM::PSPGInGaussIntegralXFEM(MatrixXd &Klocal, VectorXd &Flocal,VDOUBLE2D &dNPdr,VDOUBLE2D &x_current,const int numOfNodeInElm,const double weight,const int ic)
{  
  int ii, jj;
  int TI,TIp1,TIp2,TIp3;
  int TJ,TJp1,TJp2,TJp3;

  VDOUBLE2D dNPdx(numOfNodeInElm,VDOUBLE1D(3,0));

  double dxdr[3][3];
  MathFEM::calc_dxdr(dxdr,dNPdr,x_current,numOfNodeInElm);
  double detJ=dxdr[0][0]*dxdr[1][1]*dxdr[2][2]+dxdr[0][1]*dxdr[1][2]*dxdr[2][0]+dxdr[0][2]*dxdr[1][0]*dxdr[2][1]
             -dxdr[0][2]*dxdr[1][1]*dxdr[2][0]-dxdr[0][1]*dxdr[1][0]*dxdr[2][2]-dxdr[0][0]*dxdr[1][2]*dxdr[2][1];
  
  MathFEM::calc_dNdx(dNPdx,dNPdr,dxdr,numOfNodeInElm);

  VDOUBLE2D K(numOfNodeInElm,VDOUBLE1D(numOfNodeInElm,0));
  
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
      for(int k=0;k<3;k++){
        K[ii][jj] += dNPdx[ii][k]*dNPdx[jj][k];
      }
      Klocal(TIp3, TJp3) += tau * K[ii][jj] * detJ * weight;
    }
  }
}

void FEM::LocalRefinement(const int n, VDOUBLE1D &b)
{
  int i=1,tmp=0;
  for(int k=0;k<n;k++){
    if(n-i>0){
    b[tmp]   = n-i;
    b[tmp+1] = -(n-i);
    tmp+=2;
    }else if(n-i==0){
      b[tmp]=0;
      tmp++;
    }else{
      break;
    }
    i += 2;
  }
}


