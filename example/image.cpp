#include <cstdio>
#include <iostream>
#include <cmath>

int main(void){
  
  int nx,ny,nz;
  nx = 4,ny=4,nz=4;
  double Lx,Ly,Lz;
  Lx = 1e0; Ly = 1e0; Lz=1e0;
  double dx,dy,dz;
  dx = Lx/(double)nx;
  dy = Ly/(double)ny;
  dz = Lz/(double)nz;

  double x1 = 31.5e0 * dx;
  double x2 = 96.5e0 * dy;

  double value;
  
  FILE *fp; 
  fp=fopen("image4x4x4_rectangle.dat","w");
    for(int k=0; k<nz; k++){
    for(int j=0;j<ny;j++){
      for(int i=0;i<nx;i++){
        if(k>=nz/4 && k<nz*3/4){
          if(i>=nx/4 && i<nx*3/4){
            value = 1;
          }else{
            value = 0;    
          }
        }else{
          value = 0;
        }
        fprintf(fp,"%d %d %d %e\n",i,j,k,value);
      }
    }
  }

  fp=fopen("image128x128x128_around.dat","w");
  for(int k=0; k<nz; k++){
    for(int j=0; j<ny; j++){
      for(int i=0; i<nx; i++){
        value = 1;
        fprintf(fp,"%d %d %d %e\n",i,j,k,value);
      }
    }
  }

  fp=fopen("sdf128x128x128.dat","w");
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