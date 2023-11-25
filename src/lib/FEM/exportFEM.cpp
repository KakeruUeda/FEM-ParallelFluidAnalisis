#include <iostream>
#include <string>
#include <cmath>

#include "FEM.h"

using namespace std;

void FEM::export_vti(const string &file, vector<int> &node, vector<int> &element)
{
  FILE *fp; 
  fp=fopen(file.c_str(),"w");

  if(fp==NULL){
    cout <<file << " open error" << endl;
    exit(1);

  }
  fprintf(fp,"<?xml version=\"1.0\"?>\n");
  fprintf(fp,"<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp,"<ImageData WholeExtent= \"%d %d %d %d %d %d\" Origin= \"%e %e %e\" Spacing= \"%e %e %e\" >\n",0,nx,0,ny,0,nz,0e0,0e0,0e0,dx,dy,dz);  
  fprintf(fp,"<Piece Extent= \"%d %d %d %d %d %d\">\n",0,nx,0,ny,0,nz);  
  fprintf(fp,"<PointData>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"nodeProcId\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<numOfNodeGlobal;i++){
    fprintf(fp,"%d\n",node[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");
  fprintf(fp,"<CellData>\n");
  fprintf(fp,"<DataArray type=\"Float32\" Name=\"elmProcId\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<numOfElmGlobal;i++){
    fprintf(fp,"%d\n",element[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</CellData>\n");

  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</ImageData>\n");
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
}