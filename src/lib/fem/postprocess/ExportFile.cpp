#include "ExportFile.h"

void ExportFile::export_vti_metis(const string &file, VINT1D &node, VINT1D &element, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz)
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
  for(int i=0; i<(nx+1)*(ny+1)*(nz+1); i++){
    fprintf(fp, "%d\n", node[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");
  fprintf(fp,"<CellData>\n");
  fprintf(fp,"<DataArray type=\"Float32\" Name=\"elmProcId\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0; i<nx*ny*nz; i++){
    fprintf(fp, "%d\n", element[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</CellData>\n");

  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</ImageData>\n");
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
}


void ExportFile::export_vti_node(const string &file, VDOUBLE1D &node, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz)
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
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"node\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<(nx+1)*(ny+1)*(nz+1);i++){
    fprintf(fp,"%e\n", node[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");
  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</ImageData>\n");
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
}


void ExportFile::export_vti_elm(const string &file, VDOUBLE1D &element, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz)
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
  fprintf(fp,"<CellData>\n");
  fprintf(fp,"<DataArray type=\"Float32\" Name=\"elm\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<nx*ny*nz;i++){
    fprintf(fp,"%e\n", element[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</CellData>\n");

  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</ImageData>\n");
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
}


void ExportFile::export_vti_node_3dir(const string &file, VDOUBLE2D &node, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz)
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
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"node\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0; i<(nx+1)*(ny+1)*(nz+1); i++){
    fprintf(fp,"%e %e %e\n", node[i][0], node[i][1], node[i][2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");
  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</ImageData>\n");
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
}



void ExportFile::export_vti_node_xyz(const string &file, VDOUBLE3D &node, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz)
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
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"node\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int k=0; k<nz+1; k++){
    for(int j=0; j<ny+1; j++){
      for(int i=0; i<nx+1; i++){
        fprintf(fp,"%e \n", node[k][j][i]);
      }
    }
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");
  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</ImageData>\n");
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
}


void ExportFile::export_vti_elm_xyz_3dir(const string &file, VDOUBLE3D &element1, VDOUBLE3D &element2, VDOUBLE3D &element3, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz)
{
  FILE *fp;
  fp=fopen(file.c_str(),"w");

  if(fp==NULL){
    cout <<file << " open error" << endl;
    exit(1);
  }
  fprintf(fp,"<?xml version=\"1.0\"?>\n");
  fprintf(fp,"<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp,"<ImageData WholeExtent= \"%d %d %d %d %d %d\" Origin= \"%e %e %e\" Spacing= \"%e %e %e\" >\n", 0, nx, 0, ny, 0, nz, 0e0, 0e0, 0e0, dx, dy, dz);  
  fprintf(fp,"<Piece Extent= \"%d %d %d %d %d %d\">\n", 0, nx, 0, ny, 0, nz);
  fprintf(fp,"<CellData>\n");
  fprintf(fp,"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int k=0; k<nz; k++){
    for(int j=0; j<ny; j++){
      for(int i=0; i<nx; i++){
        fprintf(fp,"%e %e %e\n", element1[k][j][i], element2[k][j][i], element3[k][j][i]);
      }
    }
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</CellData>\n");
  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</ImageData>\n");
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
}


void ExportFile::export_vti_elm_xyz(const string &file, VDOUBLE3D &element, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz)
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
  fprintf(fp,"<CellData>\n");
  fprintf(fp,"<DataArray type=\"Float32\" Name=\"elm\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int k=0; k<nz; k++){
    for(int j=0; j<ny; j++){
      for(int i=0; i<nx; i++){
        fprintf(fp,"%e \n", element[k][j][i]);
      }
    }
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</CellData>\n");

  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</ImageData>\n");
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
}


void ExportFile::export_vti_domain(const string &file, VDOUBLE1D &sdf, VDOUBLE1D &phi, VDOUBLE1D &phiEX, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz)
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
  for(int i=0;i<(nx+1)*(ny+1)*(nz+1);i++){
    fprintf(fp,"%e\n", sdf[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");
  fprintf(fp,"<CellData>\n");
  fprintf(fp,"<DataArray type=\"Float32\" Name=\"phi\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<nx*ny*nz;i++){
    fprintf(fp,"%e\n", phi[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float32\" Name=\"phiEX\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<nx*ny*nz;i++){
    fprintf(fp,"%e\n", phiEX[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</CellData>\n");
  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</ImageData>\n");
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
}


void ExportFile::export_vti_result(const string &file, VDOUBLE1D &u, VDOUBLE1D &v, VDOUBLE1D &w, VDOUBLE1D &p, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz)
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
  for(int i=0;i<(nx+1)*(ny+1)*(nz+1);i++){
    fprintf(fp,"%e %e %e\n",u[i],v[i],w[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"pressure\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<(nx+1)*(ny+1)*(nz+1);i++){
    fprintf(fp,"%e\n",p[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");
  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</ImageData>\n");
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
}


void ExportFile::export_vti_result_2D(const string &file, VDOUBLE1D &u, VDOUBLE1D &v, VDOUBLE1D &p, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz)
{
  FILE *fp; 
  fp=fopen(file.c_str(),"w");

  if(fp==NULL){
    cout <<file << " open error" << endl;
    exit(1);
  }

  fprintf(fp,"<?xml version=\"1.0\"?>\n");
  fprintf(fp,"<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp,"<ImageData WholeExtent= \"%d %d %d %d %d %d\" Origin= \"%e %e %e\" Spacing= \"%e %e %e\" >\n",0,nx,0,ny,0,1,0e0,0e0,dx*nx/2,dx,dy,dz);  
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
    fprintf(fp,"%e\n", p[i]);
  }
  for(int i=nz/2*(nx+1)*(ny+1);i<nz/2*(nx+1)*(ny+1)+(nx+1)*(ny+1);i++){
    fprintf(fp,"%e\n", p[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");
  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</ImageData>\n");
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
}


void ExportFile::export_vti_velocity_cell(const string &file, VDOUBLE4D &vel, int nx, const int ny, const int nz, const double dx, const double dy, const double dz)
{
  FILE *fp;
  fp=fopen(file.c_str(),"w");

  if(fp==NULL){
    cout <<file << " open error" << endl;
    exit(1);
  }
  fprintf(fp,"<?xml version=\"1.0\"?>\n");
  fprintf(fp,"<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp,"<ImageData WholeExtent= \"%d %d %d %d %d %d\" Origin= \"%e %e %e\" Spacing= \"%e %e %e\" >\n", 0, nx, 0, ny, 0, nz, 0e0, 0e0, 0e0, dx, dy, dz);  
  fprintf(fp,"<Piece Extent= \"%d %d %d %d %d %d\">\n", 0, nx, 0, ny, 0, nz);
  fprintf(fp,"<CellData>\n");
  fprintf(fp,"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int k=0; k<nz; k++){
    for(int j=0; j<ny; j++){
      for(int i=0; i<nx; i++){
        fprintf(fp,"%e %e %e\n", vel[k][j][i][0], vel[k][j][i][1], vel[k][j][i][2]);
      }
    }
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</CellData>\n");
  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</ImageData>\n");
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
}


void ExportFile::export_vti_velocity_cell(const string &file, VDOUBLE3D &u, VDOUBLE3D &v, VDOUBLE3D &w, int nx, const int ny, const int nz, const double dx, const double dy, const double dz)
{
  FILE *fp;
  fp=fopen(file.c_str(),"w");

  if(fp==NULL){
    cout <<file << " open error" << endl;
    exit(1);
  }
  fprintf(fp,"<?xml version=\"1.0\"?>\n");
  fprintf(fp,"<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp,"<ImageData WholeExtent= \"%d %d %d %d %d %d\" Origin= \"%e %e %e\" Spacing= \"%e %e %e\" >\n", 0, nx, 0, ny, 0, nz, 0e0, 0e0, 0e0, dx, dy, dz);  
  fprintf(fp,"<Piece Extent= \"%d %d %d %d %d %d\">\n", 0, nx, 0, ny, 0, nz);
  fprintf(fp,"<CellData>\n");
  fprintf(fp,"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int k=0; k<nz; k++){
    for(int j=0; j<ny; j++){
      for(int i=0; i<nx; i++){
        fprintf(fp,"%e %e %e\n", u[k][j][i], v[k][j][i], w[k][j][i]);
      }
    }
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</CellData>\n");
  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</ImageData>\n");
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
}



void ExportFile::export_vti_velocity_node(const string &file, VDOUBLE4D &vel, int nx, const int ny, const int nz, const double dx, const double dy, const double dz)
{
  FILE *fp;
  fp=fopen(file.c_str(),"w");

  if(fp==NULL){
    cout <<file << " open error" << endl;
    exit(1);
  }
  fprintf(fp,"<?xml version=\"1.0\"?>\n");
  fprintf(fp,"<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp,"<ImageData WholeExtent= \"%d %d %d %d %d %d\" Origin= \"%e %e %e\" Spacing= \"%e %e %e\" >\n", 0, nx, 0, ny, 0, nz, 0e0, 0e0, 0e0, dx, dy, dz);  
  fprintf(fp,"<Piece Extent= \"%d %d %d %d %d %d\">\n", 0, nx, 0, ny, 0, nz);
  fprintf(fp,"<PointData>\n");
  fprintf(fp,"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int k=0; k<nz+1; k++){
    for(int j=0; j<ny+1; j++){
      for(int i=0; i<nx+1; i++){
        fprintf(fp,"%e %e %e\n", vel[k][j][i][0], vel[k][j][i][1], vel[k][j][i][2]);
      }
    }
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");
  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</ImageData>\n");
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
}



void ExportFile::export_vti_node_2dir_2D(const string &file, VDOUBLE2D &node, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz)
{
  FILE *fp; 
  fp=fopen(file.c_str(),"w");

  if(fp==NULL){
    cout <<file << " open error" << endl;
    exit(1);
  }

  fprintf(fp,"<?xml version=\"1.0\"?>\n");
  fprintf(fp,"<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp,"<ImageData WholeExtent= \"%d %d %d %d %d %d\" Origin= \"%e %e %e\" Spacing= \"%e %e %e\" >\n",0,nx,0,ny,0,1,0e0,0e0,dx*nx/2,dx,dy,dz);  
  fprintf(fp,"<Piece Extent= \"%d %d %d %d %d %d\">\n",0,nx,0,ny,0,1);  
  fprintf(fp,"<PointData>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"node\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=nz/2*(nx+1)*(ny+1);i<nz/2*(nx+1)*(ny+1)+(nx+1)*(ny+1);i++){
    fprintf(fp,"%e %e 0e0\n",node[i][0], node[i][1]);
  }
  for(int i=nz/2*(nx+1)*(ny+1);i<nz/2*(nx+1)*(ny+1)+(nx+1)*(ny+1);i++){
    fprintf(fp,"%e %e 0e0\n", node[i][0], node[i][1]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");
  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</ImageData>\n");
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
}

