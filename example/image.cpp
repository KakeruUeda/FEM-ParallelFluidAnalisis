#include <cstdio>
#include <iostream>
#include <cmath>

int main(void){

  int nx,ny,nz;
  nx = 16,ny=16,nz=32;
  double Lx,Ly,Lz;
  Lx = 1e0; Ly = 1e0; Lz=2e0;
  double dx,dy,dz;
  dx = Lx/(double)nx;
  dy = Ly/(double)ny;
  dz = Lz/(double)nz;

  double x1 = 31.5e0 * dx;
  double x2 = 96.5e0 * dy;

  double value;

  FILE *fp; 
  fp=fopen("image32x32x32.dat","w");
    for(int k=0; k<nz; k++){
    for(int j=0;j<ny;j++){
      for(int i=0;i<nx;i++){
        value = 1;
        fprintf(fp,"%d %d %d %e\n",i,j,k,value);
      }
    }
  }

  fp=fopen("image32x32x32_around.dat","w");
  for(int k=0; k<nz; k++){
    for(int j=0; j<ny; j++){
      for(int i=0; i<nx; i++){
        value = 1;
        fprintf(fp,"%d %d %d %e\n",i,j,k,value);
      }
    }
  }

  fp=fopen("sdf32x32x32.dat","w");
  for(int k=0; k<nz; k++){
    for(int j=0; j<ny; j++){
      for(int i=0; i<nx; i++){
        value = 0;
        fprintf(fp,"%d %d %d %e\n",i,j,k,value);
      }
    }
  }
  fclose(fp);

  return 0;
}