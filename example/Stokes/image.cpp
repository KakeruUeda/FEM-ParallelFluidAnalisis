#include <cstdio>
#include <iostream>
#include <cmath>

int main(void){

  int nx,ny,nz;
  nx = 32,ny=200,nz=32;
  double value;

  FILE *fp; 
  fp=fopen("image512x512x512.dat","w");
  for(int k=0;k<nz;k++){
    for(int j=0;j<ny;j++){
      for(int i=0;i<nx;i++){
      // if(i>=16 && i<48){
      //   value = 1;
      // }else{
      //   value = 0;
      // }
      // if(i==15 || i==48) value = 0.5;
        value = 1e0;
        fprintf(fp,"%d %d %d %e\n",i,j,k,value);
      }
    }
  }

  fclose(fp);

  return 0;
}