#include "VariationalDataAssimilation.h"

void VarDA::calcCostFunction()
{      
  costFunction.term1 = 0e0;
  costFunction.term2 = 0e0;
  costFunction.term3 = 0e0;

  nx_in_voxel = obs.dx/dx;
  ny_in_voxel = obs.dy/dy;
  nz_in_voxel = obs.dz/dz;

  volume = obs.dx * obs.dy * obs.dz;

  if(estVel == EstimatedVelocity::AVERAGE){
    getAveragedVelocityFromCFD();
  }
  else if(estVel == EstimatedVelocity::CENTER){
    getCenterVelocityFromCFD();
  }
  else{
    PetscPrintf(MPI_COMM_WORLD, "\n estVel not defined \n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // COST FUNCTION TERM1
  for(int k=0; k<obs.nz; k++){
    for(int j=0; j<obs.ny; j++){
      for(int i=0; i<obs.nx; i++){
       
        ue[k][j][i] = obs.u[k][j][i] - cfd.u[k][j][i];
        ve[k][j][i] = obs.v[k][j][i] - cfd.v[k][j][i];
        we[k][j][i] = obs.w[k][j][i] - cfd.w[k][j][i];
        
        double dev = ue[k][j][i] * ue[k][j][i] + ve[k][j][i] * ve[k][j][i] + we[k][j][i] * we[k][j][i];
        costFunction.term1 += 5e-1 * alpha_cf * dev * volume;
      }
    }
  }

   VDOUBLE1D N2D(4, 0e0);
   VDOUBLE2D dNdr2D(4, VDOUBLE1D(2, 0e0));
   VDOUBLE2D x_current2D(4, VDOUBLE1D(2, 0e0));
   Gauss gauss(2);

   // COST FUNCTION TERM2
   for(int ic=0; ic<control_elm_node.size(); ic++){
     for(int i=0; i<4; i++){
       for(int j=0; j<2; j++){
         x_current2D[i][j] = x[control_elm_node[ic][i]][bdface_dir[j]];
       }
     }
     double value = 0e0;
     for(int i1=0; i1<2; i1++){
       for(int i2=0; i2<2; i2++){
         ShapeFunction2D::C2D4_N(N2D, gauss.point[i1], gauss.point[i2]);
         ShapeFunction2D::C2D4_dNdr(dNdr2D, gauss.point[i1], gauss.point[i2]);
         ReguralizationTerm2_inGaussIntegral(N2D, dNdr2D, x_current2D, value, gauss.weight[i1]*gauss.weight[i2], ic);
       }
     }
     costFunction.term2 += 5e-1 * beta1_cf * value;
   }

   // COST FUNCTION TERM3
   for(int ic=0; ic<control_elm_node.size(); ic++){
     for(int i=0; i<4; i++){
       for(int j=0; j<2; j++){
         x_current2D[i][j] = x[control_elm_node[ic][i]][bdface_dir[j]];
       }
     }
     double value = 0e0;
     for(int i1=0; i1<2; i1++){
       for(int i2=0; i2<2; i2++){
         ShapeFunction2D::C2D4_N(N2D, gauss.point[i1], gauss.point[i2]);
         ShapeFunction2D::C2D4_dNdr(dNdr2D, gauss.point[i1], gauss.point[i2]);
         ReguralizationTerm3_inGaussIntegral(N2D, dNdr2D, x_current2D, value, gauss.weight[i1]*gauss.weight[i2], ic);
       }
     }
     costFunction.term3 += 5e-1 * beta2_cf * value;
   }

   costFunction.sum();

   return;
}


void VarDA::getAveragedVelocityFromCFD()
{
  double mic = 1e-10; 
  VDOUBLE4D vel_local_in_voxel(nz_in_voxel+1, VDOUBLE3D(ny_in_voxel+1, VDOUBLE2D(nx_in_voxel+1, VDOUBLE1D(3, 0e0))));
  VDOUBLE4D x_local_in_voxel(nz_in_voxel+1, VDOUBLE3D(ny_in_voxel+1, VDOUBLE2D(nx_in_voxel+1, VDOUBLE1D(3, 0e0))));

  for(int k=0; k<obs.nz; k++){
    for(int j=0; j<obs.ny; j++){
      for(int i=0; i<obs.nx; i++){

        for(int r=0; r<nz_in_voxel+1; r++){
          for(int q=0; q<ny_in_voxel+1; q++){
            for(int p=0; p<nx_in_voxel+1; p++){
              
              double px = i * obs.dx + p * obs.dx / nx_in_voxel;
              double py = j * obs.dy + q * obs.dy / ny_in_voxel; 
              double pz = k * obs.dz + r * obs.dz / nz_in_voxel;

              int ix = px / dx + mic;
              int iy = py / dy + mic;
              int iz = pz / dz + mic;

              if(ix == nx && iy != ny && iz != nz){
                ix = ix - 1;
              }else if(ix != nx && iy == ny && iz != nz){
                iy = iy - 1;
              }else if(ix == nx && iy == ny && iz != nz){
                ix = ix - 1; iy = iy - 1;
              }else if(ix == nx && iy != ny && iz == nz){
                ix = ix - 1; iz = iz - 1;
              }else if(ix != nx && iy == ny && iz == nz){
                iy = iy - 1; iz = iz - 1;
              }else if(ix == nx && iy == ny && iz == nz){
                ix = ix - 1; iy = iy - 1; iz = iz - 1;
              }else if(iz == nz){
                iz = iz - 1;
              }

              double x_trans = px - ((ix * dx) + (dx/2e0));
              double y_trans = py - ((iy * dy) + (dy/2e0));
              double z_trans = pz - ((iz * dz) + (dz/2e0));
              
              x_trans = x_trans / (dx / 2e0);
              y_trans = y_trans / (dy / 2e0);
              z_trans = z_trans / (dz / 2e0);
  
              if(x_trans<-1-mic || x_trans>1+mic){
                PetscPrintf(MPI_COMM_WORLD, "\n x_trans interpolation error found. x_trans = %e \n", x_trans);
                MPI_Abort(MPI_COMM_WORLD, 1);
              }
              if(y_trans<-1-mic || y_trans>1+mic){
                PetscPrintf(MPI_COMM_WORLD, "\n y_trans interpolation error found. y_trans = %e \n", y_trans);
                MPI_Abort(MPI_COMM_WORLD, 1);
              }
              if(z_trans<-1-mic || z_trans>1+mic){
                PetscPrintf(MPI_COMM_WORLD, "\n z_trans interpolation error found. z_trans = %e \n", z_trans);
                MPI_Abort(MPI_COMM_WORLD, 1);
              }

              int elm = ix + (iy * nx) + (iz * nx * ny);

              if(elm > numOfElmGlobal){
                PetscPrintf(MPI_COMM_WORLD, "\n elm error \n");
                MPI_Abort(MPI_COMM_WORLD, 1);
              }
            
              VDOUBLE1D N(numOfNodeInElm, 0e0);
              ShapeFunction3D::C3D8_N(N,x_trans,y_trans,z_trans);

              for(int l=0; l<3; l++){
                vel_local_in_voxel[r][q][p][l] = 0e0;
                x_local_in_voxel[r][q][p][l]   = 0e0;
              }

              for(int e=0; e<numOfNodeInElm; e++)
              { 
                vel_local_in_voxel[r][q][p][0] += N[e] * u[element[elm][e]];
                vel_local_in_voxel[r][q][p][1] += N[e] * v[element[elm][e]];
                vel_local_in_voxel[r][q][p][2] += N[e] * w[element[elm][e]];
                x_local_in_voxel[r][q][p][0]   += N[e] * x[element[elm][e]][0];
                x_local_in_voxel[r][q][p][1]   += N[e] * x[element[elm][e]][1];
                x_local_in_voxel[r][q][p][2]   += N[e] * x[element[elm][e]][2];
              }
            }
          }
        }

       for(int r=0; r<nz_in_voxel; r++){
          for(int q=0; q<ny_in_voxel; q++){
            for(int p=0; p<nx_in_voxel; p++){
              
              VDOUBLE2D vel_current(numOfNodeInElm, VDOUBLE1D(3,0e0));
              VDOUBLE2D x_current(numOfNodeInElm, VDOUBLE1D(3,0e0));
              
              for(int k=0; k<3; k++){
                vel_current[0][k] = vel_local_in_voxel[r][q][p][k];     
                vel_current[1][k] = vel_local_in_voxel[r][q][p+1][k];   
                vel_current[2][k] = vel_local_in_voxel[r][q+1][p+1][k];   
                vel_current[3][k] = vel_local_in_voxel[r][q+1][p][k];     
                vel_current[4][k] = vel_local_in_voxel[r+1][q][p][k];     
                vel_current[5][k] = vel_local_in_voxel[r+1][q][p+1][k];   
                vel_current[6][k] = vel_local_in_voxel[r+1][q+1][p+1][k]; 
                vel_current[7][k] = vel_local_in_voxel[r+1][q+1][p][k];     
              }

              for(int k=0; k<3; k++){
                x_current[0][k] = x_local_in_voxel[r][q][p][k];     
                x_current[1][k] = x_local_in_voxel[r][q][p+1][k];   
                x_current[2][k] = x_local_in_voxel[r][q+1][p+1][k];   
                x_current[3][k] = x_local_in_voxel[r][q+1][p][k];     
                x_current[4][k] = x_local_in_voxel[r+1][q][p][k];     
                x_current[5][k] = x_local_in_voxel[r+1][q][p+1][k];   
                x_current[6][k] = x_local_in_voxel[r+1][q+1][p+1][k]; 
                x_current[7][k] = x_local_in_voxel[r+1][q+1][p][k];   
              }
              
              VDOUBLE1D N(numOfNodeInElm, 0e0);
              VDOUBLE2D dNdr(numOfNodeInElm, VDOUBLE1D(3, 0e0));
              Gauss gauss(2);

              VDOUBLE1D vel_local_ave_in_voxel(3, 0e0);
              
              for(int i1=0; i1<2; i1++){
                for(int i2=0; i2<2; i2++){
                  for(int i3=0; i3<2; i3++){
                    ShapeFunction3D::C3D8_N(N, gauss.point[i1], gauss.point[i2], gauss.point[i3]);
                    ShapeFunction3D::C3D8_dNdr(dNdr, gauss.point[i1], gauss.point[i2], gauss.point[i3]);
                    gaussIntegral(N, dNdr, x_current, vel_current, vel_local_ave_in_voxel, gauss.weight[i1]*gauss.weight[i2]*gauss.weight[i3]);
                  }
                }
              }
              cfd.u[k][j][i] += vel_local_ave_in_voxel[0];
              cfd.v[k][j][i] += vel_local_ave_in_voxel[1];
              cfd.w[k][j][i] += vel_local_ave_in_voxel[2];
            }
          }
        }
        cfd.u[k][j][i] /= volume;
        cfd.v[k][j][i] /= volume;
        cfd.w[k][j][i] /= volume;
      }
    }
  }
}

void VarDA::getCenterVelocityFromCFD()
{
  exit(1);
}


void VarDA::gaussIntegral(VDOUBLE1D &N, VDOUBLE2D &dNdr, VDOUBLE2D &x_current, VDOUBLE2D &u_current, VDOUBLE1D &vel_local_ave_in_voxel, const double weight)
{
  double dxdr[3][3], drdx[3][3];

  MathFEM::calc_dxdr(dxdr, dNdr, x_current, numOfNodeInElm);

  double detJ = dxdr[0][0]*dxdr[1][1]*dxdr[2][2]+dxdr[0][1]*dxdr[1][2]*dxdr[2][0]+dxdr[0][2]*dxdr[1][0]*dxdr[2][1]
               -dxdr[0][2]*dxdr[1][1]*dxdr[2][0]-dxdr[0][1]*dxdr[1][0]*dxdr[2][2]-dxdr[0][0]*dxdr[1][2]*dxdr[2][1];

  for(int k=0; k<3; k++){
     for(int p=0; p<numOfNodeInElm; p++){
      vel_local_ave_in_voxel[k] += N[p] * u_current[p][k] * detJ * weight;
    }
  }
  
  return;
}


void VarDA::ReguralizationTerm2_inGaussIntegral(VDOUBLE1D &N, VDOUBLE2D &dNdr, VDOUBLE2D &x_current, double &value, const double weight, const int ic)
{
  double dxdr[2][2]; 
  VDOUBLE2D dNdx(4, VDOUBLE1D(2, 0e0));

  MathFEM::calc_dxdr2D(dxdr, dNdr, x_current, 4);
  double detJ = dxdr[0][0]*dxdr[1][1] - dxdr[0][1]*dxdr[1][0];
  MathFEM::calc_dNdx2D(dNdx, dNdr, dxdr, 4);

  double u_gp, v_gp, w_gp;

  for(int p=0; p<4; p++){
    u_gp = N[p] * u[control_elm_node[ic][p]];
    v_gp = N[p] * v[control_elm_node[ic][p]];
    w_gp = N[p] * w[control_elm_node[ic][p]];
  }

  value += (u_gp*u_gp +  v_gp*v_gp + w_gp*w_gp) * detJ * weight;

  return;
}


void VarDA::ReguralizationTerm3_inGaussIntegral(VDOUBLE1D &N, VDOUBLE2D &dNdr, VDOUBLE2D &x_current, double &value, const double weight, const int ic)
{
  double dxdr[2][2]; 
  VDOUBLE2D dNdx(4, VDOUBLE1D(2, 0e0));

  MathFEM::calc_dxdr2D(dxdr, dNdr, x_current, 4);
  double detJ = dxdr[0][0]*dxdr[1][1] - dxdr[0][1]*dxdr[1][0];
  MathFEM::calc_dNdx2D(dNdx, dNdr, dxdr, 4);

  double dudx, dudy, dvdx, dvdy, dwdx, dwdy;

  for(int p=0; p<4; p++){
    dudx = dNdx[p][0] * u[control_elm_node[ic][p]];
    dudy = dNdx[p][1] * u[control_elm_node[ic][p]];
    dvdx = dNdx[p][0] * v[control_elm_node[ic][p]];
    dvdy = dNdx[p][1] * v[control_elm_node[ic][p]];
    dwdx = dNdx[p][0] * w[control_elm_node[ic][p]];
    dwdy = dNdx[p][1] * w[control_elm_node[ic][p]];
  }

  value += ((dudx*dudx +  dvdx*dvdx + dwdx*dwdx) + (dudy*dudy +  dvdy*dvdy + dwdy*dwdy)) * detJ * weight;

  return;
}



void VarDA::calcFeedbackForce()
{
  VDOUBLE3D feedbackIntegral(numOfElmGlobal, VDOUBLE2D(numOfNodeInElm, VDOUBLE1D(3, 0e0)));
  VDOUBLE2D x_current(numOfNodeInElm, VDOUBLE1D(3, 0e0));
  
  calcEdgeVale();

  for(int ic=0; ic<numOfElmGlobal; ic++)
  {
    if(phiVOF[ic] < 1e-10) continue;
    for(int i=0; i<numOfNodeInElm; i++){
      for(int k=0; k<3; k++) x_current[i][k] = x[element[ic][i]][k];
    }
    VDOUBLE1D N(numOfNodeInElm);
    VDOUBLE2D dNdr(numOfNodeInElm, VDOUBLE1D(3, 0e0));
    Gauss gauss(2);
    double feedbackGaussPoint[3] = {0e0, 0e0, 0e0};

    for(int i1=0; i1<2; i1++){
      for(int i2=0; i2<2; i2++){
        for(int i3=0; i3<2; i3++){
          ShapeFunction3D::C3D8_N(N, gauss.point[i1], gauss.point[i2], gauss.point[i3]);
          ShapeFunction3D::C3D8_dNdr(dNdr, gauss.point[i1], gauss.point[i2], gauss.point[i3]);
          calcInterpolatedFeedback(N, x_current, feedbackGaussPoint, ic);
          feedback_inGaussIntegral(feedbackIntegral, N, dNdr, feedbackGaussPoint, x_current, gauss.weight[i1]*gauss.weight[i2]*gauss.weight[i3], ic);
        }
      }
    }
  }

  
  for(int ic=0; ic<numOfNodeGlobal; ic++){
    for(int k=0; k<3; k++) feedbackForce[ic][k] = 0e0;
  }
  for(int ic=0; ic<numOfElmGlobal; ic++){
    for(int p=0; p<numOfNodeInElm; p++){
      for(int k=0; k<3; k++){
        feedbackForce[element[ic][p]][k] += feedbackIntegral[ic][p][k];
      }
    }
  }

  return;
}

void VarDA::calcInterpolatedFeedback(VDOUBLE1D &N, VDOUBLE2D &x_current, double (&feedbackGaussPoint)[3], const int ic)
{      
  double mic = 1e-10;      
  double point[3] = {0e0, 0e0, 0e0};
  
  for(int k=0;k<3;k++){
    for(int p=0; p<numOfNodeInElm; p++) point[k] += N[p] * x_current[p][k];
  }
  
  double px = point[0] + 5e-1 * obs.dx;
  double py = point[1] + 5e-1 * obs.dy;
  double pz = point[2] + 5e-1 * obs.dz;

  int ix = px / obs.dx + mic;
  int iy = py / obs.dy + mic;
  int iz = pz / obs.dz + mic;

  double x_trans = px - ((ix * obs.dx) + (obs.dx / 2e0));
  double y_trans = py - ((iy * obs.dy) + (obs.dy / 2e0));
  double z_trans = pz - ((iz * obs.dz) + (obs.dz / 2e0));
  
  x_trans = x_trans / (obs.dx / 2e0);
  y_trans = y_trans / (obs.dy / 2e0);
  z_trans = z_trans / (obs.dz / 2e0);

  if(x_trans<-1-mic || x_trans>1+mic){
    PetscPrintf(MPI_COMM_WORLD, "\n x_trans interpolation error found. x_trans = %e %e \n", x_trans, px);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  if(y_trans<-1-mic || y_trans>1+mic){
    PetscPrintf(MPI_COMM_WORLD, "\n y_trans interpolation error found. y_trans = %e \n", y_trans);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  if(z_trans<-1-mic || z_trans>1+mic){
    PetscPrintf(MPI_COMM_WORLD, "\n z_trans interpolation error found. z_trans = %e \n", z_trans);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  VDOUBLE1D NFeedBack(numOfNodeInElm, 0e0);
  ShapeFunction3D::C3D8_N(NFeedBack, x_trans, y_trans, z_trans);

  feedbackGaussPoint[0] = NFeedBack[0]*ue_edge[iz][iy][ix]+NFeedBack[1]*ue_edge[iz][iy][ix+1]+NFeedBack[2]*ue_edge[iz][iy+1][ix+1]+NFeedBack[3]*ue_edge[iz][iy+1][ix]
                        + NFeedBack[4]*ue_edge[iz+1][iy][ix]+NFeedBack[5]*ue_edge[iz+1][iy][ix+1]+NFeedBack[6]*ue_edge[iz+1][iy+1][ix+1]+NFeedBack[7]*ue_edge[iz+1][iy+1][ix];
  feedbackGaussPoint[1] = NFeedBack[0]*ve_edge[iz][iy][ix]+NFeedBack[1]*ve_edge[iz][iy][ix+1]+NFeedBack[2]*ve_edge[iz][iy+1][ix+1]+NFeedBack[3]*ve_edge[iz][iy+1][ix]
                        + NFeedBack[4]*ve_edge[iz+1][iy][ix]+NFeedBack[5]*ve_edge[iz+1][iy][ix+1]+NFeedBack[6]*ve_edge[iz+1][iy+1][ix+1]+NFeedBack[7]*ve_edge[iz+1][iy+1][ix];
  feedbackGaussPoint[2] = NFeedBack[0]*we_edge[iz][iy][ix]+NFeedBack[1]*we_edge[iz][iy][ix+1]+NFeedBack[2]*we_edge[iz][iy+1][ix+1]+NFeedBack[3]*we_edge[iz][iy+1][ix]
                        + NFeedBack[4]*we_edge[iz+1][iy][ix]+NFeedBack[5]*we_edge[iz+1][iy][ix+1]+NFeedBack[6]*we_edge[iz+1][iy+1][ix+1]+NFeedBack[7]*we_edge[iz+1][iy+1][ix];

  return;
}


void VarDA::calcEdgeVale()
{
  for(int k=0; k<obs.nz+2; k++){
    for(int j=0; j<obs.ny+2; j++){
      for(int i=0; i<obs.nx+2; i++){ 
        if(k == 0){
          if((j > 0 && j < obs.ny+1) && (i > 0 && i < obs.nx+1)){
            ue_edge[k][j][i] = ue[k][j-1][i-1] + (ue[k][j-1][i-1] - ue[k+1][j-1][i-1]);
            ve_edge[k][j][i] = ve[k][j-1][i-1] + (ve[k][j-1][i-1] - ve[k+1][j-1][i-1]);
            we_edge[k][j][i] = we[k][j-1][i-1] + (we[k][j-1][i-1] - we[k+1][j-1][i-1]);
          }
        }else if(k > 0 && k < obs.nz+1){
          if(j == 0 && i != 0 && i != obs.nx+1){
            ue_edge[k][j][i] = ue[k-1][j][i-1] + (ue[k-1][j][i-1] - ue[k-1][j+1][i-1]);
            ve_edge[k][j][i] = ve[k-1][j][i-1] + (ve[k-1][j][i-1] - ve[k-1][j+1][i-1]);
            we_edge[k][j][i] = we[k-1][j][i-1] + (we[k-1][j][i-1] - we[k-1][j+1][i-1]);
          }else if(i == 0 && j != 0 && j != obs.ny+1){
            ue_edge[k][j][i] = ue[k-1][j-1][i] + (ue[k-1][j-1][i] - ue[k-1][j-1][i+1]);
            ve_edge[k][j][i] = ve[k-1][j-1][i] + (ve[k-1][j-1][i] - ve[k-1][j-1][i+1]);
            we_edge[k][j][i] = we[k-1][j-1][i] + (we[k-1][j-1][i] - we[k-1][j-1][i+1]);            
          }else if(i == obs.nx+1 && j != 0 && j != obs.ny+1){
            ue_edge[k][j][i] = ue[k-1][j-1][i-2] + (ue[k-1][j-1][i-2] - ue[k-1][j-1][i-3]);
            ve_edge[k][j][i] = ve[k-1][j-1][i-2] + (ve[k-1][j-1][i-2] - ve[k-1][j-1][i-3]);
            we_edge[k][j][i] = we[k-1][j-1][i-2] + (we[k-1][j-1][i-2] - we[k-1][j-1][i-3]);            
          }else if(j == obs.ny+1 && i != 0 && i != obs.nx+1){
            ue_edge[k][j][i] = ue[k-1][j-2][i-1] + (ue[k-1][j-2][i-1] - ue[k-1][j-3][i-1]);
            ve_edge[k][j][i] = ve[k-1][j-2][i-1] + (ve[k-1][j-2][i-1] - ve[k-1][j-3][i-1]);
            we_edge[k][j][i] = we[k-1][j-2][i-1] + (we[k-1][j-2][i-1] - we[k-1][j-3][i-1]);
          }else if(i != 0 && i != obs.nx+1 && j != 0 && j != obs.ny+1){
            ue_edge[k][j][i] = ue[k-1][j-1][i-1];
            ve_edge[k][j][i] = ve[k-1][j-1][i-1];
            we_edge[k][j][i] = we[k-1][j-1][i-1];
          }
        }else if(k == obs.nz+1){
          if((j > 0 && j < obs.ny+1) && (i > 0 && i < obs.nx+1)){
            ue_edge[k][j][i] = ue[k-2][j-1][i-1] + (ue[k-2][j-1][i-1]- ue[k-3][j-1][i-1]);
            ve_edge[k][j][i] = ve[k-2][j-1][i-1] + (ve[k-2][j-1][i-1]- ve[k-3][j-1][i-1]);
            we_edge[k][j][i] = we[k-2][j-1][i-1] + (we[k-2][j-1][i-1]- we[k-3][j-1][i-1]);
          }
        }

      } 
    }   
  }

  for(int k=0; k<obs.nz+2; k++){
    for(int j=0; j<obs.ny+2; j++){
      for(int i=0; i<obs.nx+2; i++){
        
        if(k != 0 && k != obs.nz+1){
          if(i == 0 && j == 0){
            ue_edge[k][j][i] = ( (ue_edge[k][j][i+1] + (ue_edge[k][j][i+1] - ue_edge[k][j][i+2])) + (ue_edge[k][j+1][i] + (ue_edge[k][j+1][i] - ue_edge[k][j+2][i])) ) / 2e0;
            ve_edge[k][j][i] = ( (ve_edge[k][j][i+1] + (ve_edge[k][j][i+1] - ve_edge[k][j][i+2])) + (ve_edge[k][j+1][i] + (ve_edge[k][j+1][i] - ve_edge[k][j+2][i])) ) / 2e0;
            we_edge[k][j][i] = ( (we_edge[k][j][i+1] + (we_edge[k][j][i+1] - we_edge[k][j][i+2])) + (we_edge[k][j+1][i] + (we_edge[k][j+1][i] - we_edge[k][j+2][i])) ) / 2e0;        
          }else if(i == obs.nx+1 && j == 0){
            ue_edge[k][j][i] = ( (ue_edge[k][j][i-1] + (ue_edge[k][j][i-1] - ue_edge[k][j][i-2])) + (ue_edge[k][j+1][i] + (ue_edge[k][j+1][i] - ue_edge[k][j+2][i])) ) / 2e0;
            ve_edge[k][j][i] = ( (ve_edge[k][j][i-1] + (ve_edge[k][j][i-1] - ve_edge[k][j][i-2])) + (ve_edge[k][j+1][i] + (ve_edge[k][j+1][i] - ve_edge[k][j+2][i])) ) / 2e0;
            we_edge[k][j][i] = ( (we_edge[k][j][i-1] + (we_edge[k][j][i-1] - we_edge[k][j][i-2])) + (we_edge[k][j+1][i] + (we_edge[k][j+1][i] - we_edge[k][j+2][i])) ) / 2e0;              
          }else if(i == 0 && j == obs.ny+1){
            ue_edge[k][j][i] = ( (ue_edge[k][j][i+1] + (ue_edge[k][j][i+1] - ue_edge[k][j][i+2])) + (ue_edge[k][j-1][i] + (ue_edge[k][j-1][i] - ue_edge[k][j-2][i])) ) / 2e0;
            ve_edge[k][j][i] = ( (ve_edge[k][j][i+1] + (ve_edge[k][j][i+1] - ve_edge[k][j][i+2])) + (ve_edge[k][j-1][i] + (ve_edge[k][j-1][i] - ve_edge[k][j-2][i])) ) / 2e0;
            we_edge[k][j][i] = ( (we_edge[k][j][i+1] + (we_edge[k][j][i+1] - we_edge[k][j][i+2])) + (we_edge[k][j-1][i] + (we_edge[k][j-1][i] - we_edge[k][j-2][i])) ) / 2e0;        
          }else if(i == obs.nx+1 && j == obs.ny+1){
            ue_edge[k][j][i] = ( (ue_edge[k][j][i-1] + (ue_edge[k][j][i-1] - ue_edge[k][j][i-2])) + (ue_edge[k][j-1][i] + (ue_edge[k][j-1][i] - ue_edge[k][j-2][i])) ) / 2e0;
            ve_edge[k][j][i] = ( (ve_edge[k][j][i-1] + (ve_edge[k][j][i-1] - ve_edge[k][j][i-2])) + (ve_edge[k][j-1][i] + (ve_edge[k][j-1][i] - ve_edge[k][j-2][i])) ) / 2e0;
            we_edge[k][j][i] = ( (we_edge[k][j][i-1] + (we_edge[k][j][i-1] - we_edge[k][j][i-2])) + (we_edge[k][j-1][i] + (we_edge[k][j-1][i] - we_edge[k][j-2][i])) ) / 2e0;              
          }
        }
      }
    }
  }

  for(int k=0; k<obs.nz+2; k++){
    for(int j=0; j<obs.ny+2; j++){
      for(int i=0; i<obs.nx+2; i++){
  
        if(k == 0){
          if(i == 0 && j == 0){
            ue_edge[k][j][i] = ( (ue_edge[k][j][i+1] + (ue_edge[k][j][i+1] - ue_edge[k][j][i+2])) + (ue_edge[k][j+1][i] + (ue_edge[k][j+1][i] - ue_edge[k][j+2][i])) + (ue_edge[k+1][j][i] + (ue_edge[k+1][j][i] - ue_edge[k+2][j][i])) ) / 3e0;
            ve_edge[k][j][i] = ( (ve_edge[k][j][i+1] + (ve_edge[k][j][i+1] - ve_edge[k][j][i+2])) + (ve_edge[k][j+1][i] + (ve_edge[k][j+1][i] - ve_edge[k][j+2][i])) + (ve_edge[k+1][j][i] + (ve_edge[k+1][j][i] - ve_edge[k+2][j][i])) ) / 3e0;
            we_edge[k][j][i] = ( (we_edge[k][j][i+1] + (we_edge[k][j][i+1] - we_edge[k][j][i+2])) + (we_edge[k][j+1][i] + (we_edge[k][j+1][i] - we_edge[k][j+2][i])) + (we_edge[k+1][j][i] + (we_edge[k+1][j][i] - we_edge[k+2][j][i])) ) / 3e0;        
          }else if(i == obs.nx+1 && j == 0){
            ue_edge[k][j][i] = ( (ue_edge[k][j][i-1] + (ue_edge[k][j][i-1] - ue_edge[k][j][i-2])) + (ue_edge[k][j+1][i] + (ue_edge[k][j+1][i] - ue_edge[k][j+2][i])) + (ue_edge[k+1][j][i] + (ue_edge[k+1][j][i] - ue_edge[k+2][j][i])) ) / 3e0;
            ve_edge[k][j][i] = ( (ve_edge[k][j][i-1] + (ve_edge[k][j][i-1] - ve_edge[k][j][i-2])) + (ve_edge[k][j+1][i] + (ve_edge[k][j+1][i] - ve_edge[k][j+2][i])) + (ve_edge[k+1][j][i] + (ve_edge[k+1][j][i] - ve_edge[k+2][j][i])) ) / 3e0;
            we_edge[k][j][i] = ( (we_edge[k][j][i-1] + (we_edge[k][j][i-1] - we_edge[k][j][i-2])) + (we_edge[k][j+1][i] + (we_edge[k][j+1][i] - we_edge[k][j+2][i])) + (we_edge[k+1][j][i] + (we_edge[k+1][j][i] - we_edge[k+2][j][i])) ) / 3e0;              
          }else if(i == 0 && j == obs.ny+1){
            ue_edge[k][j][i] = ( (ue_edge[k][j][i+1] + (ue_edge[k][j][i+1] - ue_edge[k][j][i+2])) + (ue_edge[k][j-1][i] + (ue_edge[k][j-1][i] - ue_edge[k][j-2][i])) + (ue_edge[k+1][j][i] + (ue_edge[k+1][j][i] - ue_edge[k+2][j][i])) ) / 3e0;
            ve_edge[k][j][i] = ( (ve_edge[k][j][i+1] + (ve_edge[k][j][i+1] - ve_edge[k][j][i+2])) + (ve_edge[k][j-1][i] + (ve_edge[k][j-1][i] - ve_edge[k][j-2][i])) + (ue_edge[k+1][j][i] + (ue_edge[k+1][j][i] - ue_edge[k+2][j][i])) ) / 3e0;
            we_edge[k][j][i] = ( (we_edge[k][j][i+1] + (we_edge[k][j][i+1] - we_edge[k][j][i+2])) + (we_edge[k][j-1][i] + (we_edge[k][j-1][i] - we_edge[k][j-2][i])) + (ue_edge[k+1][j][i] + (ue_edge[k+1][j][i] - ue_edge[k+2][j][i])) ) / 3e0;        
          }else if(i == obs.nx+1 && j == obs.ny+1){
            ue_edge[k][j][i] = ( (ue_edge[k][j][i-1] + (ue_edge[k][j][i-1] - ue_edge[k][j][i-2])) + (ue_edge[k][j-1][i] + (ue_edge[k][j-1][i] - ue_edge[k][j-2][i])) + (ue_edge[k+1][j][i] + (ue_edge[k+1][j][i] - ue_edge[k+2][j][i])) ) / 3e0;
            ve_edge[k][j][i] = ( (ve_edge[k][j][i-1] + (ve_edge[k][j][i-1] - ve_edge[k][j][i-2])) + (ve_edge[k][j-1][i] + (ve_edge[k][j-1][i] - ve_edge[k][j-2][i])) + (ve_edge[k+1][j][i] + (ve_edge[k+1][j][i] - ve_edge[k+2][j][i])) ) / 3e0;
            we_edge[k][j][i] = ( (we_edge[k][j][i-1] + (we_edge[k][j][i-1] - we_edge[k][j][i-2])) + (we_edge[k][j-1][i] + (we_edge[k][j-1][i] - we_edge[k][j-2][i])) + (we_edge[k+1][j][i] + (we_edge[k+1][j][i] - we_edge[k+2][j][i])) ) / 3e0;              
          }
        }else if(k == obs.nz+1){
          if(i == 0 && j == 0){
            ue_edge[k][j][i] = ( (ue_edge[k][j][i+1] + (ue_edge[k][j][i+1] - ue_edge[k][j][i+2])) + (ue_edge[k][j+1][i] + (ue_edge[k][j+1][i] - ue_edge[k][j+2][i])) + (ue_edge[k-1][j][i] + (ue_edge[k-1][j][i] - ue_edge[k-2][j][i])) ) / 3e0;
            ve_edge[k][j][i] = ( (ve_edge[k][j][i+1] + (ve_edge[k][j][i+1] - ve_edge[k][j][i+2])) + (ve_edge[k][j+1][i] + (ve_edge[k][j+1][i] - ve_edge[k][j+2][i])) + (ve_edge[k-1][j][i] + (ve_edge[k-1][j][i] - ve_edge[k-2][j][i])) ) / 3e0;
            we_edge[k][j][i] = ( (we_edge[k][j][i+1] + (we_edge[k][j][i+1] - we_edge[k][j][i+2])) + (we_edge[k][j+1][i] + (we_edge[k][j+1][i] - we_edge[k][j+2][i])) + (we_edge[k-1][j][i] + (we_edge[k-1][j][i] - we_edge[k-2][j][i])) ) / 3e0;        
          }else if(i == obs.nx+1 && j == 0){
            ue_edge[k][j][i] = ( (ue_edge[k][j][i-1] + (ue_edge[k][j][i-1] - ue_edge[k][j][i-2])) + (ue_edge[k][j+1][i] + (ue_edge[k][j+1][i] - ue_edge[k][j+2][i])) + (ue_edge[k-1][j][i] + (ue_edge[k-1][j][i] - ue_edge[k-2][j][i])) ) / 3e0;
            ve_edge[k][j][i] = ( (ve_edge[k][j][i-1] + (ve_edge[k][j][i-1] - ve_edge[k][j][i-2])) + (ve_edge[k][j+1][i] + (ve_edge[k][j+1][i] - ve_edge[k][j+2][i])) + (ue_edge[k-1][j][i] + (ve_edge[k-1][j][i] - ve_edge[k-2][j][i])) ) / 3e0;
            we_edge[k][j][i] = ( (we_edge[k][j][i-1] + (we_edge[k][j][i-1] - we_edge[k][j][i-2])) + (we_edge[k][j+1][i] + (we_edge[k][j+1][i] - we_edge[k][j+2][i])) + (ue_edge[k-1][j][i] + (we_edge[k-1][j][i] - we_edge[k-2][j][i])) ) / 3e0;              
          }else if(i == 0 && j == obs.ny+1){
            ue_edge[k][j][i] = ( (ue_edge[k][j][i+1] + (ue_edge[k][j][i+1] - ue_edge[k][j][i+2])) + (ue_edge[k][j-1][i] + (ue_edge[k][j-1][i] - ue_edge[k][j-2][i])) + (ue_edge[k-1][j][i] + (ue_edge[k-1][j][i] - ue_edge[k-2][j][i])) ) / 3e0;
            ve_edge[k][j][i] = ( (ve_edge[k][j][i+1] + (ve_edge[k][j][i+1] - ve_edge[k][j][i+2])) + (ve_edge[k][j-1][i] + (ve_edge[k][j-1][i] - ve_edge[k][j-2][i])) + (ve_edge[k-1][j][i] + (ve_edge[k-1][j][i] - ve_edge[k-2][j][i])) ) / 3e0;
            we_edge[k][j][i] = ( (we_edge[k][j][i+1] + (we_edge[k][j][i+1] - we_edge[k][j][i+2])) + (we_edge[k][j-1][i] + (we_edge[k][j-1][i] - we_edge[k][j-2][i])) + (we_edge[k-1][j][i] + (we_edge[k-1][j][i] - we_edge[k-2][j][i])) ) / 3e0;        
          }else if(i == obs.nx+1 && j == obs.ny+1){
            ue_edge[k][j][i] = ( (ue_edge[k][j][i-1] + (ue_edge[k][j][i-1] - ue_edge[k][j][i-2])) + (ue_edge[k][j-1][i] + (ue_edge[k][j-1][i] - ue_edge[k][j-2][i])) + (ue_edge[k-1][j][i] + (ue_edge[k-1][j][i] - ue_edge[k-2][j][i])) ) / 3e0;
            ve_edge[k][j][i] = ( (ve_edge[k][j][i-1] + (ve_edge[k][j][i-1] - ve_edge[k][j][i-2])) + (ve_edge[k][j-1][i] + (ve_edge[k][j-1][i] - ve_edge[k][j-2][i])) + (ve_edge[k-1][j][i] + (ve_edge[k-1][j][i] - ve_edge[k-2][j][i])) ) / 3e0;
            we_edge[k][j][i] = ( (we_edge[k][j][i-1] + (we_edge[k][j][i-1] - we_edge[k][j][i-2])) + (we_edge[k][j-1][i] + (we_edge[k][j-1][i] - we_edge[k][j-2][i])) + (we_edge[k-1][j][i] + (we_edge[k-1][j][i] - we_edge[k-2][j][i])) ) / 3e0;              
          }
        }
      }
    }
  }

  return;
}


void VarDA::feedback_inGaussIntegral(VDOUBLE3D &feedbackIntegral, VDOUBLE1D &N, VDOUBLE2D &dNdr, double (&feedbackGaussPoint)[3], VDOUBLE2D &x_current, const double weight, const int ic)
{
  double dxdr[3][3];
  MathFEM::calc_dxdr(dxdr, dNdr, x_current, numOfNodeInElm);
  double detJ = dxdr[0][0]*dxdr[1][1]*dxdr[2][2]+dxdr[0][1]*dxdr[1][2]*dxdr[2][0]+dxdr[0][2]*dxdr[1][0]*dxdr[2][1]
               -dxdr[0][2]*dxdr[1][1]*dxdr[2][0]-dxdr[0][1]*dxdr[1][0]*dxdr[2][2]-dxdr[0][0]*dxdr[1][2]*dxdr[2][1];

  double area = obs.dx * obs.dy * obs.dz;
  for(int p=0; p<numOfNodeInElm; p++){
    feedbackIntegral[ic][p][0] += area * alpha_cf * N[p] * feedbackGaussPoint[0] * detJ * weight;
    feedbackIntegral[ic][p][1] += area * alpha_cf * N[p] * feedbackGaussPoint[1] * detJ * weight;
    feedbackIntegral[ic][p][2] += area * alpha_cf * N[p] * feedbackGaussPoint[2] * detJ * weight;
  }

  return;
}
