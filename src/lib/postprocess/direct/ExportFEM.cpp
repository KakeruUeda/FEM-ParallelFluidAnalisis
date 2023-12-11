#include "FEM.h"

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


void FEM::export_vti_domain(const string &file)
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
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"sdf\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<numOfNodeGlobal;i++){
    fprintf(fp,"%e\n",sdf[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");
  fprintf(fp,"<CellData>\n");
  fprintf(fp,"<DataArray type=\"Float32\" Name=\"phi\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<numOfElmGlobal;i++){
    fprintf(fp,"%e\n",phi[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float32\" Name=\"phiEX\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<numOfElmGlobal;i++){
    fprintf(fp,"%e\n",phiEX[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</CellData>\n");
  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</ImageData>\n");
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
}

void FEM::export_vti_result(const std::string &file, vector<double> &u, vector<double> &v, vector<double> &w, vector<double> &p)
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
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNodeGlobal;i++){
    fprintf(fp,"%e %e %e\n",u[i],v[i],w[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"pressure\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<numOfNodeGlobal;i++){
    fprintf(fp,"%e\n",p[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");
  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</ImageData>\n");
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
}


void FEM::export_vti_result_2D(const std::string &file, vector<double> &u, vector<double> &v, vector<double> &p)
{
  FILE *fp; 
  fp=fopen(file.c_str(),"w");

  if(fp==NULL){
    cout <<file << " open error" << endl;
    exit(1);
  }
  fprintf(fp,"<?xml version=\"1.0\"?>\n");
  fprintf(fp,"<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp,"<ImageData WholeExtent= \"%d %d %d %d %d %d\" Origin= \"%e %e %e\" Spacing= \"%e %e %e\" >\n",0,nx,0,ny,0,1,0e0,0e0,Lx/2,dx,dy,dz);  
  fprintf(fp,"<Piece Extent= \"%d %d %d %d %d %d\">\n",0,nx,0,ny,0,1);  
  fprintf(fp,"<PointData>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=nz/2*(nx+1)*(ny+1);i<nz/2*(nx+1)*(ny+1)+(nx+1)*(ny+1);i++){
    fprintf(fp,"%e %e 0e0\n",u[i],v[i]);
  }
  for(int i=nz/2*(nx+1)*(ny+1);i<nz/2*(nx+1)*(ny+1)+(nx+1)*(ny+1);i++){
    fprintf(fp,"%e %e 0e0\n",u[i],v[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"pressure\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=nz/2*(nx+1)*(ny+1);i<nz/2*(nx+1)*(ny+1)+(nx+1)*(ny+1);i++){
    fprintf(fp,"%e\n",p[i]);
  }
  for(int i=nz/2*(nx+1)*(ny+1);i<nz/2*(nx+1)*(ny+1)+(nx+1)*(ny+1);i++){
    fprintf(fp,"%e\n",p[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");
  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</ImageData>\n");
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
}