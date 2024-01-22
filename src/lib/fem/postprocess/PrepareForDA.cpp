#include "PostProcess.h"

void PostProcess::prepareForDataAssimilation(FEM &fem, const int xNumOfVoxels, const int yNumOfVoxels, const int zNumOfVoxels)
{ 
  nx_obs = xNumOfVoxels;
  ny_obs = yNumOfVoxels;
  nz_obs = zNumOfVoxels;
  
  dx_obs = Lx_obs/nx_obs;
  dy_obs = Ly_obs/ny_obs;
  dz_obs = Lz_obs/nz_obs;
  
  nx_in_voxel = dx_obs/fem.dx;
  ny_in_voxel = dy_obs/fem.dy;
  nz_in_voxel = dz_obs/fem.dz;

  if(nx_in_voxel == 0 || ny_in_voxel  == 0 || nz_in_voxel  == 0){
     PetscPrintf(MPI_COMM_WORLD, "\n Observed voxels size needs to be bigger than fem element size \n");
     MPI_Abort(MPI_COMM_WORLD, 1);
  }

  numOfElmInVoxel = nx_in_voxel * ny_in_voxel * nz_in_voxel;
  numOfNodeInVoxel = (nx_in_voxel+1) * (ny_in_voxel+1) * (nz_in_voxel+1);

  VDOUBLE4D vel_ave_voxel(nz_obs, VDOUBLE3D(ny_obs, VDOUBLE2D(nx_obs, VDOUBLE1D(3, 0e0))));

  makeVoxelAveragedVelocity(fem, vel_ave_voxel);
  writeAveragedVelocityToFile(fem, vel_ave_voxel);

  return;
}


void PostProcess::readDAParam()
{
  string str,base_label,label;
  double tmp[3], tmp2[3];
  int tmp3[3];

  base_label = "/Postprocess";
      
  label = base_label + "/origin";
  if ( !tp_post.getInspectedVector(label, tmp, 3)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  x0 = tmp[0];
  y0 = tmp[1];
  z0 = tmp[2];

  label = base_label + "/length";
  if ( !tp_post.getInspectedVector(label, tmp2, 3)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  Lx_obs = tmp2[0];
  Ly_obs = tmp2[1];
  Lz_obs = tmp2[2];

  label = base_label + "/nx_opt";
  if ( !tp_post.getInspectedVector(label, tmp3, 3)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  nx_opt = tmp3[0];
  ny_opt = tmp3[1];
  nz_opt = tmp3[2];

  int num = 1;
  label = base_label + "/numOfVoxels" + to_string(num);
  while(tp_post.getInspectedVector(label, tmp3, 3))
  { 
    VINT1D tmp3_tmp(3);
    tmp3_tmp[0] = tmp3[0]; tmp3_tmp[1] = tmp3[1]; tmp3_tmp[2] = tmp3[2];
    numOfObsVoxels.push_back(tmp3_tmp);
    
    num++;
    label = base_label + "/numOfVoxels" + to_string(num);
  }

  return;
}


void PostProcess::validateDADomain(FEM &fem, const int xNumOfVoxels, const int yNumOfVoxels, const int zNumOfVoxels)
{
  double mic = 1e-10; 
  
  nx_obs = xNumOfVoxels;
  ny_obs = yNumOfVoxels;
  nz_obs = zNumOfVoxels;
  
  dx_obs = Lx_obs/nx_obs;
  dy_obs = Ly_obs/ny_obs;
  dz_obs = Lz_obs/nz_obs;
  
  nx_in_voxel = dx_obs/fem.dx;
  ny_in_voxel = dy_obs/fem.dy;
  nz_in_voxel = dz_obs/fem.dz;

  if(nx_in_voxel == 0 || ny_in_voxel  == 0 || nz_in_voxel  == 0){
     PetscPrintf(MPI_COMM_WORLD, "\n Observed voxels size needs to be bigger than fem element size \n");
     MPI_Abort(MPI_COMM_WORLD, 1);
  }

  for(int k=0; k<nz_obs; k++){
    for(int j=0; j<ny_obs; j++){
      for(int i=0; i<nx_obs; i++){
        
        for(int r=0; r<nz_in_voxel+1; r++){
          for(int q=0; q<ny_in_voxel+1; q++){
            for(int p=0; p<nx_in_voxel+1; p++){

              double px = x0 + i * dx_obs + p * dx_obs / nx_in_voxel;
              double py = y0 + j * dy_obs + q * dy_obs / ny_in_voxel; 
              double pz = z0 + k * dz_obs + r * dz_obs / nz_in_voxel;

              int ix = px / fem.dx + mic;
              int iy = py / fem.dy + mic;
              int iz = pz / fem.dz + mic;

              if(ix == fem.nx || iy == fem.ny || iz == fem.nz){
                PetscPrintf(MPI_COMM_WORLD, "\n Define the DA domain not to include edges \n");
                MPI_Abort(MPI_COMM_WORLD, 1);
              }
            }
          }
        }
      
      }
    }
  }

  dx_opt = Lx_obs/nx_opt;
  dy_opt = Ly_obs/ny_opt;
  dz_opt = Lz_obs/nz_opt;

  // VALIDATION FOR EXTRACTING PHIVOF :THIS IS A CELL VALUE.
  int tmp = 0;
  int nx_inside = fem.nx - 1;
  int ny_inside = fem.ny - 1;
  int nz_inside = fem.nz - 1;

  for(int k=0; k<nz_opt; k++){
    for(int j=0; j<ny_opt; j++){
      for(int i=0; i<nx_opt; i++){

        double px = x0 + i * dx_opt;
        double py = y0 + j * dy_opt;
        double pz = z0 + k * dz_opt;

        px = px - 5e-1 * fem.dx;
        py = py - 5e-1 * fem.dy; 
        pz = pz - 5e-1 * fem.dz;

        int ix = px / fem.dx + mic;
        int iy = py / fem.dy + mic;
        int iz = pz / fem.dz + mic;

        if(ix >= nx_inside || iy >= ny_inside || iz >= nz_inside){
          PetscPrintf(MPI_COMM_WORLD, "\n Define the DA domain so the extracted phiVOF can be interpolated \n");
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
      
      }
    }
  }



}


void PostProcess::makeVoxelAveragedVelocity(FEM &fem, VDOUBLE4D &vel_ave_voxel)
{
  VDOUBLE4D vel_local_in_voxel(nz_in_voxel+1, VDOUBLE3D(ny_in_voxel+1, VDOUBLE2D(nx_in_voxel+1, VDOUBLE1D(3, 0e0))));
  VDOUBLE4D x_local_in_voxel(nz_in_voxel+1, VDOUBLE3D(ny_in_voxel+1, VDOUBLE2D(nx_in_voxel+1, VDOUBLE1D(3, 0e0))));
  
  for(int k=0; k<nz_obs; k++){
    for(int j=0; j<ny_obs; j++){
      for(int i=0; i<nx_obs; i++){

        int tmp, tmp2 = 0;
        double mic = 1e-10; 
        
        for(int r=0; r<nz_in_voxel+1; r++){
          for(int q=0; q<ny_in_voxel+1; q++){
            for(int p=0; p<nx_in_voxel+1; p++){

              double px = x0 + i * dx_obs + p * dx_obs / nx_in_voxel;
              double py = y0 + j * dy_obs + q * dy_obs / ny_in_voxel; 
              double pz = z0 + k * dz_obs + r * dz_obs / nz_in_voxel;

              int ix = px / fem.dx + mic;
              int iy = py / fem.dy + mic;
              int iz = pz / fem.dz + mic;

              if(ix == fem.nx || iy == fem.ny || iz == fem.nz){
                PetscPrintf(MPI_COMM_WORLD, "\n please define the DA domain as boundaries of the CFD domain are not included \n");
                MPI_Abort(MPI_COMM_WORLD, 1);
              }
  
              double s = px - ((ix * fem.dx) + (fem.dx/2e0));
              double t = py - ((iy * fem.dy) + (fem.dy/2e0));
              double u = pz - ((iz * fem.dz) + (fem.dz/2e0));
              
              s = s / (fem.dx / 2e0);
              t = t / (fem.dy / 2e0);
              u = u / (fem.dz / 2e0);
  
              if(s<-1-mic || s>1+mic){
                PetscPrintf(MPI_COMM_WORLD, "\n s interpolation error found. s = %e \n", s);
                MPI_Abort(MPI_COMM_WORLD, 1);
              }
              if(t<-1-mic || t>1+mic){
                PetscPrintf(MPI_COMM_WORLD, "\n t interpolation error found. t = %e \n", t);
                MPI_Abort(MPI_COMM_WORLD, 1);
              }
              if(u<-1-mic || u>1+mic){
                PetscPrintf(MPI_COMM_WORLD, "\n u interpolation error found. u = %e \n", u);
                MPI_Abort(MPI_COMM_WORLD, 1);
              }

              int elm = ix + (iy * fem.nx) + (iz * fem.nx * fem.ny);

              if(elm > fem.numOfElmGlobal){
                PetscPrintf(MPI_COMM_WORLD, "\n elm error \n");
                MPI_Abort(MPI_COMM_WORLD, 1);
              }
              
              VDOUBLE1D N(fem.numOfNodeInElm, 0e0);
              ShapeFunction3D::C3D8_N(N,s,t,u);

              for(int l=0; l<3; l++){
                vel_local_in_voxel[r][q][p][l] = 0e0;
                x_local_in_voxel[r][q][p][l]   = 0e0;
              }

              for(int e=0; e<fem.numOfNodeInElm; e++)
              { 
                vel_local_in_voxel[r][q][p][0] += N[e] * fem.u[fem.element[elm][e]];
                vel_local_in_voxel[r][q][p][1] += N[e] * fem.v[fem.element[elm][e]];
                vel_local_in_voxel[r][q][p][2] += N[e] * fem.w[fem.element[elm][e]];
                x_local_in_voxel[r][q][p][0]   += N[e] * fem.x[fem.element[elm][e]][0];
                x_local_in_voxel[r][q][p][1]   += N[e] * fem.x[fem.element[elm][e]][1];
                x_local_in_voxel[r][q][p][2]   += N[e] * fem.x[fem.element[elm][e]][2];
              }
            }
          }
        } // LOCAL DIVISION LOOP

        
        for(int r=0; r<nz_in_voxel; r++){
          for(int q=0; q<ny_in_voxel; q++){
            for(int p=0; p<nx_in_voxel; p++){
              
              VDOUBLE2D vel_current(fem.numOfNodeInElm, VDOUBLE1D(3,0e0));
              VDOUBLE2D x_current(fem.numOfNodeInElm, VDOUBLE1D(3,0e0));
              
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
              
              VDOUBLE1D N(fem.numOfNodeInElm, 0e0);
              VDOUBLE2D dNdr(fem.numOfNodeInElm, VDOUBLE1D(3, 0e0));
              Gauss gauss(2);

              VDOUBLE1D vel_local_ave_in_voxel(3, 0e0);

              for(int i1=0; i1<2; i1++){
                for(int i2=0; i2<2; i2++){
                  for(int i3=0; i3<2; i3++){
                    ShapeFunction3D::C3D8_N(N, gauss.point[i1], gauss.point[i2], gauss.point[i3]);
                    ShapeFunction3D::C3D8_dNdr(dNdr, gauss.point[i1], gauss.point[i2], gauss.point[i3]);
                    gaussIntegral(fem, N, dNdr, x_current, vel_current, vel_local_ave_in_voxel, gauss.weight[i1]*gauss.weight[i2]*gauss.weight[i3]);
                  }
                }
              }
              
              for(int l=0; l<3; l++){
                vel_ave_voxel[k][j][i][l] += vel_local_ave_in_voxel[l];
              }
            
            }
          }
        }/// LOCAL ELEMENT LOOP
        
        double volume = dx_obs * dy_obs * dz_obs;
        for(int l=0; l<3; l++){
          vel_ave_voxel[k][j][i][l] /= volume;
        }
        for(int l=0; l<3; l++){
          if(fabs(vel_ave_voxel[k][j][i][l]) < 1e-10) vel_ave_voxel[k][j][i][l] = 0e0;
        }
      
      }
    }
  }/// DA DOMAIN LOOP

  return;
}


void PostProcess::gaussIntegral(FEM &fem, VDOUBLE1D &N, VDOUBLE2D &dNdr, VDOUBLE2D &x_current, VDOUBLE2D &u_current, VDOUBLE1D &vel_local_ave_in_voxel, const double weight)
{
  double dxdr[3][3], drdx[3][3];

  MathFEM::calc_dxdr(dxdr, dNdr, x_current, fem.numOfNodeInElm);

  double detJ = dxdr[0][0]*dxdr[1][1]*dxdr[2][2]+dxdr[0][1]*dxdr[1][2]*dxdr[2][0]+dxdr[0][2]*dxdr[1][0]*dxdr[2][1]
               -dxdr[0][2]*dxdr[1][1]*dxdr[2][0]-dxdr[0][1]*dxdr[1][0]*dxdr[2][2]-dxdr[0][0]*dxdr[1][2]*dxdr[2][1];

  for(int k=0; k<3; k++){
     for(int p=0; p<fem.numOfNodeInElm; p++){
      vel_local_ave_in_voxel[k] += N[p] * u_current[p][k] * detJ * weight;
    }
  }
  
  return;
}


void PostProcess::writeAveragedVelocityToFile(FEM &fem, VDOUBLE4D &vel_ave_voxel) 
{
  if(fem.myId > 0) return;

  string vtiFile = fem.outputDirDA + "/velocity_MRI_" + to_string(nx_obs) + "x" + to_string(ny_obs) + "x" + to_string(nz_obs) + ".vti";
  export_file_post.export_vti_velocity_cell(vtiFile, vel_ave_voxel, nx_obs, ny_obs, nz_obs, dx_obs, dy_obs, dz_obs);
  
  ofstream out_ave(fem.outputDirDA + "/velocity_MRI_" + to_string(nx_obs) + "x" + to_string(ny_obs) + "x" + to_string(nz_obs) + ".dat");
  out_ave << nx_obs << " " << ny_obs << " " << nz_obs << endl;
  out_ave << Lx_obs << " " << Ly_obs << " " << Lz_obs << endl;

  for (int k = 0; k < nz_obs; k++) {
      for (int j = 0; j < ny_obs; j++) {
          for (int i = 0; i < nx_obs; i++) {
              out_ave << i << " " << j << " " << k << " " 
                      << vel_ave_voxel[k][j][i][0] << " " 
                      << vel_ave_voxel[k][j][i][1] << " " 
                      << vel_ave_voxel[k][j][i][2] << endl;
          }
      }
  }
  out_ave.close();

  return;
}

void PostProcess::extractDomain(FEM &fem)
{
  dx_opt = Lx_obs/nx_opt;
  dy_opt = Ly_obs/ny_opt;
  dz_opt = Lz_obs/nz_opt;

  vel_ref_opt.resize(nz_opt+1, VDOUBLE3D(ny_opt+1, VDOUBLE2D(nx_opt+1, VDOUBLE1D(3, 0e0))));
  x_opt.resize(nz_opt+1, VDOUBLE3D(ny_opt+1, VDOUBLE2D(nx_opt+1, VDOUBLE1D(3, 0e0))));
  sdf_opt.resize(nz_opt+1, VDOUBLE2D(ny_opt+1, VDOUBLE1D(nx_opt+1, 0e0)));
  phiVOF_opt.resize(nz_opt, VDOUBLE2D(ny_opt, VDOUBLE1D(nx_opt, 0e0)));
  
  double mic = 1e-10;

  for(int k=0; k<nz_opt+1; k++){
    for(int j=0; j<ny_opt+1; j++){
      for(int i=0; i<nx_opt+1; i++){

        double px = x0 + i * dx_opt;
        double py = y0 + j * dy_opt; 
        double pz = z0 + k * dz_opt;

        int ix = px / fem.dx + mic;
        int iy = py / fem.dy + mic;
        int iz = pz / fem.dz + mic;

        if(ix == fem.nx || iy == fem.ny || iz == fem.nz){
          PetscPrintf(MPI_COMM_WORLD, "\n please define the DA domain as boundaries of the CFD domain are not included \n");
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
  
        double s = px - ((ix * fem.dx) + (fem.dx/2e0));
        double t = py - ((iy * fem.dy) + (fem.dy/2e0));
        double u = pz - ((iz * fem.dz) + (fem.dz/2e0));
        
        s = s / (fem.dx / 2e0);
        t = t / (fem.dy / 2e0);
        u = u / (fem.dz / 2e0);
  
        if(s<-1-mic || s>1+mic){
          PetscPrintf(MPI_COMM_WORLD, "\n s interpolation error found. s = %e \n", s);
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if(t<-1-mic || t>1+mic){
          PetscPrintf(MPI_COMM_WORLD, "\n t interpolation error found. t = %e \n", t);
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if(u<-1-mic || u>1+mic){
          PetscPrintf(MPI_COMM_WORLD, "\n u interpolation error found. u = %e \n", u);
          MPI_Abort(MPI_COMM_WORLD, 1);
        }

        int elm = ix + (iy * fem.nx) + (iz * fem.nx * fem.ny);

        if(elm > fem.numOfElmGlobal){
          PetscPrintf(MPI_COMM_WORLD, "\n elm error \n");
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        VDOUBLE1D N(fem.numOfNodeInElm, 0e0);
        ShapeFunction3D::C3D8_N(N,s,t,u);
        
        sdf_opt[k][j][i] = 1e0;
        for(int l=0; l<3; l++){  
          vel_ref_opt[k][j][i][l] = 0e0;
        }

        for(int e=0; e<fem.numOfNodeInElm; e++)
        { 
          sdf_opt[k][j][i]        += N[e] * fem.sdf[fem.element[elm][e]];
          vel_ref_opt[k][j][i][0] += N[e] * fem.u[fem.element[elm][e]];
          vel_ref_opt[k][j][i][1] += N[e] * fem.v[fem.element[elm][e]];
          vel_ref_opt[k][j][i][2] += N[e] * fem.w[fem.element[elm][e]];
        }
        if(fabs(sdf_opt[k][j][i]) < 1e-10) sdf_opt[k][j][i] = 0e0;
        if(fabs(vel_ref_opt[k][j][i][0]) < 1e-10) vel_ref_opt[k][j][i][0] = 0e0;
        if(fabs(vel_ref_opt[k][j][i][1]) < 1e-10) vel_ref_opt[k][j][i][1] = 0e0;        
        if(fabs(vel_ref_opt[k][j][i][2]) < 1e-10) vel_ref_opt[k][j][i][2] = 0e0;
      }
    }
  }

  int tmp = 0;
  int nx_inside = fem.nx - 1;
  int ny_inside = fem.ny - 1;
  int nz_inside = fem.nz - 1;

  int numOfInsideElmGlobal = (nx_inside)*(ny_inside)*(nz_inside);
  VINT2D element2(numOfInsideElmGlobal, VINT1D(fem.numOfNodeInElm, 0e0));
  
  for(int k=0; k<nz_inside; k++){
    for(int j=0; j<ny_inside; j++){
      for(int i=0; i<nx_inside; i++){
        element2[tmp][0]= i   + j*(nx_inside+1)     + k*(nx_inside+1)*(ny_inside+1);
        element2[tmp][1]= i+1 + j*(nx_inside+1)     + k*(nx_inside+1)*(ny_inside+1);
        element2[tmp][2]= i+1 + (j+1)*(nx_inside+1) + k*(nx_inside+1)*(ny_inside+1);
        element2[tmp][3]= i   + (j+1)*(nx_inside+1) + k*(nx_inside+1)*(ny_inside+1);
        element2[tmp][4]= i   + j*(nx_inside+1)     + (k+1)*(nx_inside+1)*(ny_inside+1);
        element2[tmp][5]= i+1 + j*(nx_inside+1)     + (k+1)*(nx_inside+1)*(ny_inside+1);
        element2[tmp][6]= i+1 + (j+1)*(nx_inside+1) + (k+1)*(nx_inside+1)*(ny_inside+1);
        element2[tmp][7]= i   + (j+1)*(nx_inside+1) + (k+1)*(nx_inside+1)*(ny_inside+1);
        tmp++;
      }
    }
  }

  for(int k=0; k<nz_opt; k++){
    for(int j=0; j<ny_opt; j++){
      for(int i=0; i<nx_opt; i++){

        double px = x0 + i * dx_opt;
        double py = y0 + j * dy_opt;
        double pz = z0 + k * dz_opt;

        px = px - 5e-1 * fem.dx;
        py = py - 5e-1 * fem.dy; 
        pz = pz - 5e-1 * fem.dz;

        int ix = px / fem.dx + mic;
        int iy = py / fem.dy + mic;
        int iz = pz / fem.dz + mic;

        if(ix == nx_inside || iy == ny_inside || iz == nz_inside){
          PetscPrintf(MPI_COMM_WORLD, "\n please define the DA domain as boundaries of the CFD domain are not included \n");
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
  
        double s = px - ((ix * fem.dx) + (fem.dx/2e0));
        double t = py - ((iy * fem.dy) + (fem.dy/2e0));
        double u = pz - ((iz * fem.dz) + (fem.dz/2e0));
        
        s = s / (fem.dx / 2e0);
        t = t / (fem.dy / 2e0);
        u = u / (fem.dz / 2e0);
  
        if(s<-1-mic || s>1+mic){
          PetscPrintf(MPI_COMM_WORLD, "\n s interpolation error found. s = %e \n", s);
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if(t<-1-mic || t>1+mic){
          PetscPrintf(MPI_COMM_WORLD, "\n t interpolation error found. t = %e \n", t);
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if(u<-1-mic || u>1+mic){
          PetscPrintf(MPI_COMM_WORLD, "\n u interpolation error found. u = %e \n", u);
          MPI_Abort(MPI_COMM_WORLD, 1);
        }

        int elm = ix + (iy * nx_inside) + (iz * nx_inside * ny_inside);

        if(elm > fem.numOfElmGlobal){
          PetscPrintf(MPI_COMM_WORLD, "\n elm error \n");
          MPI_Abort(MPI_COMM_WORLD, 1);
        }

        VDOUBLE1D N(fem.numOfNodeInElm, 0e0);
        ShapeFunction3D::C3D8_N(N,s,t,u);
        
        phiVOF_opt[k][j][i] = 0e0;

        for(int e=0; e<fem.numOfNodeInElm; e++){ 
          phiVOF_opt[k][j][i] += N[e] * fem.phiVOF[element2[elm][e]];
        }
        if(fabs(phiVOF_opt[k][j][i]) < 1e-10) phiVOF_opt[k][j][i] = 0e0;
      }
    }
  }
  writeDADomainToFile(fem);

  return;
}

double PostProcess::sdfToChi(FEM &fem, const double &sdf)
{
  double chi;
  //double eps = sqrt(fem.dx*fem.dx + fem.dy*fem.dy + fem.dz*fem.dz)/2e0;
  double eps = fem.dx;
  
  if(sdf > eps){ // FLUID
    chi = 1e0;
  }else if(sdf < -eps) { // SOLID
    chi = 0e0;
  }else{  //interface
    //chi = 1e0 - 5e-1*(1e0 + sdf/eps + sin(PI*sdf/eps)/PI);
    chi = 5e-1*(1e0 + sdf/eps + sin(PI*sdf/eps)/PI);
  }

  return chi;
}

void PostProcess::writeDADomainToFile(FEM &fem)
{
  if(fem.myId > 0) return;

  string vtiFile = fem.outputDirDA + "/sdf_opt_" + to_string(nx_opt) + "x" + to_string(ny_opt) + "x" + to_string(nz_opt) + ".vti";
  export_file_post.export_vti_node_xyz(vtiFile, sdf_opt, nx_opt, ny_opt, nz_opt, dx_opt, dy_opt, dz_opt);
  vtiFile = fem.outputDirDA + "/phiVOF_opt_" + to_string(nx_opt) + "x" + to_string(ny_opt) + "x" + to_string(nz_opt) + ".vti";
  export_file_post.export_vti_elm_xyz(vtiFile, phiVOF_opt, nx_opt, ny_opt, nz_opt, dx_opt, dy_opt, dz_opt);
  vtiFile = fem.outputDirDA + "/vel_ref_opt_" + to_string(nx_opt) + "x" + to_string(ny_opt) + "x" + to_string(nz_opt) + ".vti";
  export_file_post.export_vti_velocity_node(vtiFile, vel_ref_opt, nx_opt, ny_opt, nz_opt, dx_opt, dy_opt, dz_opt);


  ofstream out_sdf(fem.outputDirDA + "/sdf_opt_" + to_string(nx_opt) + "x" + to_string(ny_opt) + "x" + to_string(nz_opt) + ".dat");
  for(int k = 0; k < nz_opt+1; k++){
    for(int j = 0; j < ny_opt+1; j++){
      for(int i = 0; i < nx_opt+1; i++){
        out_sdf << i << " " << j << " " << k << " " << sdf_opt[k][j][i] << endl;
      }
    }
  }
  out_sdf.close();

  ofstream out_phi(fem.outputDirDA + "/phiVOF_opt_" + to_string(nx_opt) + "x" + to_string(ny_opt) + "x" + to_string(nz_opt) + ".dat");
  for(int k = 0; k < nz_opt; k++){
    for(int j = 0; j < ny_opt; j++){
      for(int i = 0; i < nx_opt; i++){
        out_phi << i << " " << j << " " << k << " " << phiVOF_opt[k][j][i] << endl;
      }
    }
  }
  out_phi.close();

  ofstream out_vel(fem.outputDirDA + "/vel_ref_opt_" + to_string(nx_opt) + "x" + to_string(ny_opt) + "x" + to_string(nz_opt) + ".dat");
  for(int k = 0; k < nz_opt; k++){
    for(int j = 0; j < ny_opt; j++){
      for(int i = 0; i < nx_opt; i++){
        out_vel << i << " " << j << " " << k << " " 
                << vel_ref_opt[k][j][i][0] << " " 
                << vel_ref_opt[k][j][i][1] << " " 
                << vel_ref_opt[k][j][i][2] << endl;
      }
    }
  }
  out_vel.close();
  
  
  return;
}